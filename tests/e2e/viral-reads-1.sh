#!/usr/bin/env bash
# Step 47 (viral DNA in the unmapped reads) on HG002's BAM from step 02, with
# synthetic viral reads added as unmapped reads. Nothing here says anything
# about HG002's own viruses: every expected call comes from the reads the case
# plants. The case:
#   - installs the viral genomes with setup.sh --viral-data (NCBI E-utilities,
#     checked against VIRAL_REF_SHA256) and finds 24 sequences;
#   - measures the human depth of the BAM itself (samtools coverage, chr1 to
#     chr22, over the bases with reads), and plants reads drawn from the HHV-6B
#     genome at half that depth (one integrated copy per cell), 300 reads of
#     CMV (a low level) and 200 reads of phiX (the lab control), with 0.3%
#     substitutions; then the step must call inherited integrated HHV-6B, CMV
#     viral DNA found but no integration, phiX the lab control, and HHV-6A,
#     whose only reads are the segments it shares with HHV-6B, not found;
#   - plants decoys that must not count as viral evidence: three HBV molecules
#     read four times each (3 counted reads: below threshold, not found), and
#     poly-A, (TTAGGG)n and (AT)n reads (VZV, where minimap2 places the AT
#     reads, keeps no counted read);
#   - removes the unmapped reads from the BAM (samtools view -F 4): the step
#     must say "not assessed" for every genome, not "not found";
#   - plants HHV-6B at a tenth of the human depth instead (a viraemia-like
#     level): the same genome, the same breadth, and the call must not be
#     integration. This is the control that the depth threshold decides;
#   - runs the Nextflow processes (VIRAL_READS, VIRAL_REPORT) on the first
#     sample, whose table and coverages must be the bash step's, byte for byte;
#   - runs generate-report.sh: the text report states the HHV-6 call.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
BASE="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
VREF="${G}/reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa"
dk() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$img" "$@"; }

# --- the viral genomes ----------------------------------------------------------------
echo "+ scripts/setup.sh --viral-data"
"${REPO}/scripts/setup.sh" --viral-data "$G" 2>&1 | tail -n 5
check "setup.sh --viral-data wrote the FASTA" test -s "$VREF"
check_eq "the FASTA is the pinned download" "$(sha256sum "$VREF" 2>/dev/null | cut -d' ' -f1)" "$VIRAL_REF_SHA256"
check_eq "24 viral genomes" "$(grep -c '^>' "$VREF" 2>/dev/null)" 24
check "the BAM of step 02 is there" test -s "${BASE}.bai"

# --- the human depth, measured here -------------------------------------------------------
dk "$SAMTOOLS_IMAGE" samtools coverage --ff UNMAP,SECONDARY,QCFAIL "$BASE" > "${T}/human_coverage.tsv"
HUMAN=$(awk -F'\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ && $5 > 0 {t += $7 * ($3 - $2 + 1); c += $5}
  END {if (c) printf "%.2f", t / c}' "${T}/human_coverage.tsv")
check "a human depth over the autosomes (${HUMAN:-none})" test -n "$HUMAN"
echo "human depth over the autosomal bases with reads: ${HUMAN}"

# decoys FASTA: unmapped SAM records that must not count as viral evidence: three HBV molecules
# read four times each (150 bp, forward, the same bases), and 30 reads each of poly-A, (TTAGGG)n
# and (AT)n.
decoys() {
  awk -v OFS='\t' '
    /^>/ {keep = ($1 == ">NC_003977.2"); next}
    keep {hbv = hbv toupper($0)}
    END {
      for (i = 0; i < 150; i++) q = q "I"
      for (m = 0; m < 3; m++) for (c = 0; c < 4; c++)
        print "dechbv_" m "_" c, 4, "*", 0, 0, "*", "*", 0, 0, substr(hbv, 101 + m * 1100, 150), q
      for (i = 0; i < 25; i++) t = t "TTAGGG"
      for (i = 0; i < 75; i++) at = at "AT"
      for (i = 0; i < 150; i++) a = a "A"
      for (i = 0; i < 30; i++) {
        print "decpa_" i, 4, "*", 0, 0, "*", "*", 0, 0, a, q
        print "dectelo_" i, 4, "*", 0, 0, "*", "*", 0, 0, t, q
        print "decat_" i, 4, "*", 0, 0, "*", "*", 0, 0, at, q
      }
    }' "$1"
}

# plant OUT FRACTION: the step-02 BAM plus unmapped reads from HHV-6B at FRACTION of the human
# depth, 300 CMV reads and 200 phiX reads (150 bp, random strand, 0.3% substitutions, seed fixed),
# and the decoys.
plant() {
  "${REPO}/bin/og" dev fixture viral-spike "$VREF" "$HUMAN" "$2" > "${T}/spike.sam"
  mkdir -p "${G}/$1/aligned" "${G}/$1/vcf"
  decoys "$VREF" > "${T}/decoy.sam"
  # The planted records go after the BAM's own unmapped reads, at the end of a coordinate-sorted file.
  # shellcheck disable=SC2016
  dk "$SAMTOOLS_IMAGE" bash -euo pipefail -c 'samtools view -h "$1" | cat - "$2" "$3" | samtools view -b -o "$4"; samtools index "$4"' \
    _ "$BASE" "${T}/spike.sam" "${T}/decoy.sam" "${G}/$1/aligned/$1_sorted.bam"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/$1/vcf/$1.vcf.gz"
  dk "$BCFTOOLS_IMAGE" bcftools index -t -f "${G}/$1/vcf/$1.vcf.gz"
  wc -l < "${T}/spike.sam"
}
# cell NAME ACCESSION COLUMN: one cell of the sample's table.
cell() {
  awk -F'\t' -v a="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $2 == a {print $c[k]}' \
    "${G}/$1/viral/$1_viral.tsv" 2>/dev/null
}
between() { awk -v x="$1" -v lo="$2" -v hi="$3" 'BEGIN {exit !(x != "" && x + 0 >= lo && x + 0 <= hi)}'; }

# --- HHV-6B at half the human depth -----------------------------------------------------
V="${SAMPLE}vir"
N_PLANTED=$(plant "$V" 0.5)
check_ge "reads planted" "$N_PLANTED" 1000
run_step 47-viral.sh "$V"
check_step_exit 47-viral.sh
cat "${G}/${V}/viral/${V}_viral.tsv" 2>/dev/null
check_eq "a row per viral genome" "$(awk 'END {print NR - 1}' "${G}/${V}/viral/${V}_viral.tsv" 2>/dev/null)" 24
check_eq "the step's human depth is the one measured here" "$(cell "$V" NC_000898.1 human_depth)" "$HUMAN"
N6B=$(grep -c '^spk6b_' "${T}/spike.sam")
# Counted reads only: the reads on HHV-6B's two copies of its direct repeat map at MAPQ 0, about
# a tenth of the genome (on a synthetic run of this plant: 13924 of 16211 reads, breadth 0.8925).
check_ge "HHV-6B: the planted reads counted (of ${N6B})" "$(cell "$V" NC_000898.1 reads)" "$((N6B * 80 / 100))"
check "HHV-6B: most of the genome covered by counted reads ($(cell "$V" NC_000898.1 breadth))" \
  between "$(cell "$V" NC_000898.1 breadth)" 0.80 1
check "HHV-6B: depth at about half the human depth ($(cell "$V" NC_000898.1 relative_depth))" \
  between "$(cell "$V" NC_000898.1 relative_depth)" 0.40 0.60
check_eq "HHV-6B: inherited integrated HHV-6 likely" "$(cell "$V" NC_000898.1 call)" "inherited integrated HHV-6 likely"
check "HHV-6B: one copy per cell, from the mode ($(cell "$V" NC_000898.1 copies_per_cell))" \
  between "$(cell "$V" NC_000898.1 copies_per_cell)" 0.80 1.20
check_eq "HHV-6A, with only the reads it shares with HHV-6B: not found" "$(cell "$V" NC_001664.4 call)" "not found"
check_ge "CMV: the 300 planted reads" "$(cell "$V" NC_006273.2 reads)" 280
check_eq "CMV: viral DNA found" "$(cell "$V" NC_006273.2 call)" "viral DNA found"
check "CMV: far below the human depth ($(cell "$V" NC_006273.2 copies_per_cell) copies per cell)" \
  between "$(cell "$V" NC_006273.2 copies_per_cell)" 0 0.05
check_eq "phiX: the lab control" "$(cell "$V" NC_001422.1 call)" "lab control (phiX)"
check_eq "HBV: 12 reads of 3 molecules are 3 counted reads" "$(cell "$V" NC_003977.2 reads)" 3
check_eq "HBV: below threshold, not viral DNA found" "$(cell "$V" NC_003977.2 call)" "below threshold"
check_eq "VZV: the simple-repeat decoys are not counted" "$(cell "$V" NC_001348.1 reads)" 0
check_eq "VZV: not found" "$(cell "$V" NC_001348.1 call)" "not found"
for f in viral_reads.bam viral_reads.bam.bai viral_coverage.tsv human_coverage.tsv; do
  check "the step wrote ${V}_${f}" test -s "${G}/${V}/viral/${V}_${f}"
done
for f in unmapped.fq viral.sam viral_depth.tsv viral_hits.tsv bam_idxstats.tsv; do
  check "the step removed ${V}_${f}" test ! -e "${G}/${V}/viral/${V}_${f}"
done

# --- the same genome at a tenth of the human depth --------------------------------------------
L10="${SAMPLE}virlow"
plant "$L10" 0.1 >/dev/null
run_step 47-viral.sh "$L10"
check_step_exit 47-viral.sh
check "control: HHV-6B at a tenth of the depth is below 0.30 ($(cell "$L10" NC_000898.1 relative_depth))" \
  between "$(cell "$L10" NC_000898.1 relative_depth)" 0 0.29
check_eq "control: and is viral DNA found, not integration" "$(cell "$L10" NC_000898.1 call)" "viral DNA found"
check "control: over most of the genome still ($(cell "$L10" NC_000898.1 breadth))" between "$(cell "$L10" NC_000898.1 breadth)" 0.5 1

# --- the same BAM without its unmapped reads ----------------------------------------------------
NU="${SAMPLE}virnu"
mkdir -p "${G}/${NU}/aligned"
# shellcheck disable=SC2016
dk "$SAMTOOLS_IMAGE" bash -euo pipefail -c 'samtools view -b -F 4 -o "$2" "$1"; samtools index "$2"' \
  _ "${G}/${V}/aligned/${V}_sorted.bam" "${G}/${NU}/aligned/${NU}_sorted.bam"
run_step 47-viral.sh "$NU"
check_step_exit 47-viral.sh
check_eq "no unmapped read in the BAM: all 24 genomes not assessed" \
  "$(awk -F'\t' 'NR > 1 && $11 == "not assessed"' "${G}/${NU}/viral/${NU}_viral.tsv" 2>/dev/null | wc -l | tr -d ' ')" 24
check_eq "and HHV-6B is not called not found" "$(cell "$NU" NC_000898.1 call)" "not assessed"

# --- the text report -------------------------------------------------------------------
run_step generate-report.sh "$V"
check_step_exit generate-report.sh
check "the text report states the HHV-6 call" grep -q 'HHV-6: inherited integrated HHV-6B likely' "${G}/${V}/${V}_report.txt"

# --- Nextflow against the bash step --------------------------------------------------------
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
SHEET="${T}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index,bam,bam_index\n%s,male,%s,%s,%s,%s\n' "$V" \
  "${G}/${V}/vcf/${V}.vcf.gz" "${G}/${V}/vcf/${V}.vcf.gz.tbi" \
  "${G}/${V}/aligned/${V}_sorted.bam" "${G}/${V}/aligned/${V}_sorted.bam.bai" > "$SHEET"
NF_OUT="${G}/nf-viral"
echo "+ nextflow run main.nf -profile docker --tools viral"
( cd "$T" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "${E2E_WORK}/nf-work-viral" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools viral --viral_ref "$VREF" \
    --sex_check warn --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${T}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
for f in viral viral_coverage human_coverage; do
  a="${G}/${V}/viral/${V}_${f}.tsv" b="${NF_OUT}/${V}/viral/${V}_${f}.tsv"
  check "Nextflow wrote ${V}_${f}.tsv" test -s "$b"
  check "Nextflow's ${V}_${f}.tsv is the bash step's" cmp -s "$a" "$b"
  cmp -s "$a" "$b" || diff "$a" "$b" | head -20
done
check "Nextflow published the viral reads" test -s "${NF_OUT}/${V}/viral/${V}_viral_reads.bam"
for f in viral_depth.tsv viral_hits.tsv bam_idxstats.tsv; do
  check "Nextflow kept ${V}_${f} out of the results" test ! -e "${NF_OUT}/${V}/viral/${V}_${f}"
done
check "control: the comparison sees a difference" bash -c '! cmp -s "$1" "$2"' _ \
  "${G}/${L10}/viral/${L10}_viral.tsv" "${NF_OUT}/${V}/viral/${V}_viral.tsv"

finish
