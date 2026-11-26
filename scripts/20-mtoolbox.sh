#!/usr/bin/env bash
# Mitochondrial analysis — Heteroplasmy detection + variant calling with GATK Mutect2
# Input: Sorted BAM
# Output: Mitochondrial VCF with heteroplasmy fractions
# Runtime: ~15-30 minutes
# Note: Originally planned for MToolBox, but no working Docker image exists.
#       GATK Mutect2 in mitochondrial mode is the standard alternative; the file
#       name is historical.
#
# chrM is circular but GRCh38 writes it as a line, so reads across its origin
# (m.16569 joined to m.1) align poorly. As Broad's mitochondria pipeline does,
# the control region (m.16025-16569 and m.1-575) is called a second time on
# chrM rotated by 8000 bp, and those calls replace the linear ones there
# (og setup chrm-shift, og vcf shift-back; mito/shifted/).
#
# NuMTs (nuclear copies of mitochondrial DNA) put reads on chrM that come from
# the autosomes. GATK's NuMTFilterTool marks an allele as possible_numt when
# its depth is what such copies could explain at the sample's median autosomal
# coverage. That median is read from step 16b's mosdepth output, or taken from
# AUTOSOMAL_COVERAGE; without either the filter runs at depth 0 and marks
# nothing: the step warns, and the annotation table's '# numt' line says so.
#
# Last, each call is annotated with its ClinVar record (clinvar/clinvar_chr.vcf.gz,
# which setup.sh installs) and with gnomAD v3.1's chrM sites (GNOMAD_MITO_URL in
# versions.env, which setup.sh installs under gnomad_mito/), into
# mito/<sample>_chrM_annotated.tsv (og mito annotate, bin/og). A
# missing file leaves its columns '.' and the table says why.
#
# Then the mtDNA copy number: chrM reads per base over autosomal reads per base
# in the BAM index (samtools idxstats), times 2, into mito/<sample>_mtdna_copies.tsv
# (not assessed under 20 reads per kb of autosome, as on a sliced BAM). And large deletions: eKLIPse
# (EKLIPSE_IMAGE) on the chrM reads, into mito/<sample>_chrM_deletions.tsv. The
# report names the variants and deletions of assets/mito_confirmed.tsv among
# the calls.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/mito"

echo "=== Mitochondrial Analysis (GATK Mutect2): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
# Tables of an earlier run go first, so a stage that fails now leaves none behind.
rm -f "${OUTPUT_DIR}/${SAMPLE}_mtdna_copies.tsv" "${OUTPUT_DIR}/${SAMPLE}_chrM_deletions.tsv"


echo "[1/8] Extracting chrM reads..."
run_in --cpus 2 --memory 4g \
  "$SAMTOOLS_IMAGE" \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM \
      > /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam && \
    samtools index /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam
  "

echo "[2/8] Checking sequence dictionary..."
if [ ! -f "$REF_DICT" ]; then
  echo "  Creating sequence dictionary..."
  # The dictionary goes next to the FASTA, so its directory is writable here.
  run_in --rw "$(dirname "$REF_FASTA")" --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk CreateSequenceDictionary \
      -R "${REF_FASTA_C}" \
      -O "$(cpath "$REF_DICT")"
else
  echo "  Sequence dictionary already exists, skipping."
fi

echo "[3/8] Running Mutect2 in mitochondrial mode..."
run_in --cpus 4 --memory 8g \
  "$GATK_IMAGE" \
  gatk Mutect2 \
    -R "${REF_FASTA_C}" \
    -I "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" \
    -L chrM \
    --mitochondria-mode \
    --max-mnp-distance 0 \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz"

# The control region again on chrM rotated by 8000 bp, where the origin
# (m.16569 joined to m.1) lies mid-sequence, so reads across it align in
# full: og setup chrm-shift writes the rotated chrM, GATK's BwaSpark realigns
# the chrM reads to it (local Spark; the GATK image has no bwa binary), and
# Mutect2 calls its control region. og vcf shift-back maps those calls back
# and merges them with the linear calls outside the control region, with a
# '##og_shifted_pass=' header line. The filters below then run on the merged
# calls with the linear pass's stats. Should any part fail, the step warns
# and filters the linear calls alone, as before; the summary then keeps the
# origin_edge flag. The large intermediate files never leave the container.
SHIFT_DIR="${OUTPUT_DIR}/shifted"
SHIFT_C="/genome/${SAMPLE}/mito/shifted"
RAW_C="/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz"
shifted_pass() {
  rm -rf "$SHIFT_DIR" && mkdir -p "$SHIFT_DIR" || return 1
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
    /pgp-bin/og setup chrm-shift --fasta "${REF_FASTA_C}" \
      --out "${SHIFT_C}/chrM_shifted.fasta" --intervals "${SHIFT_C}/control.intervals" || return 1
  # --root: BwaSpark starts Hadoop's login, which needs a passwd entry for
  # the process's uid, and the GATK image has none for yours (a
  # LoginException, "invalid null input: name"; HADOOP_USER_NAME and
  # SPARK_USER do not help). The image has no other aligner. The container
  # gives the directory back to you when it ends, so the next run can
  # remove it. Its large files stay in the container's own /tmp, so even a
  # container killed outright (no EXIT trap) leaves in the directory, which
  # is yours, only plain files you can remove. The Nextflow MITO_VARIANTS
  # runs as root for the same reason.
  # shellcheck disable=SC2016  # expanded inside the container
  run_in --root --cpus 4 --memory 8g "$GATK_IMAGE" bash -c '
    set -euo pipefail
    d=$1
    t=$(mktemp -d)
    trap '\''chown -R "$3" "$d" || true'\'' EXIT
    export SPARK_LOCAL_IP=127.0.0.1 SPARK_LOCAL_HOSTNAME=localhost
    gatk CreateSequenceDictionary -R "$d/chrM_shifted.fasta" -O "$d/chrM_shifted.dict"
    gatk BwaMemIndexImageCreator -I "$d/chrM_shifted.fasta" -O "$d/chrM_shifted.img"
    # The reads Mutect2 counts on the linear chrM, without their alignment.
    gatk PrintReads -I "$2" -L chrM -O "$t/reads.bam" \
      --read-filter NotDuplicateReadFilter --read-filter NotSecondaryAlignmentReadFilter \
      --read-filter NotSupplementaryAlignmentReadFilter --read-filter MateOnSameContigOrNoMappedMateReadFilter
    # Unaligned and grouped by name. Not --SANITIZE: it re-encodes qualities
    # that are all Q31 or more as if they were phred+64.
    gatk RevertSam -I "$t/reads.bam" -O "$t/unaligned.sam" --SORT_ORDER queryname \
      --OUTPUT_BY_READGROUP false --RESTORE_ORIGINAL_QUALITIES false --VALIDATION_STRINGENCY LENIENT \
      --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO
    # BwaSpark aligns records two by two as pairs: keep both mates or neither.
    awk -F"\t" "/^@/ { print; next } \$1 == name { print prev; print; name = \"\"; next } { name = \$1; prev = \$0 }" \
      "$t/unaligned.sam" > "$t/paired.sam"
    # Unpaired reads (single-end, HiFi or ONT) leave no pair, and Mutect2
    # would then see few or none of the reads the linear pass saw, so the
    # control region would lose its calls. Under 90% kept fails the pass.
    n=$(grep -vc "^@" "$t/unaligned.sam" || true)
    k=$(grep -vc "^@" "$t/paired.sam" || true)
    echo "  Shifted-reference pass: ${k} of ${n} chrM reads kept as complete pairs" >&2
    if [ "$k" -eq 0 ] || [ $((k * 10)) -lt $((n * 9)) ]; then
      echo "  The chrM reads are not paired (under 90% kept): the rotated pass needs read pairs" >&2
      exit 1
    fi
    gatk --java-options -Xmx4g BwaSpark -I "$t/paired.sam" -R "$d/chrM_shifted.fasta" \
      --bwa-mem-index-image "$d/chrM_shifted.img" --spark-master "local[4]" -O "$t/aligned.bam"
    gatk SortSam -I "$t/aligned.bam" -O "$t/sorted.bam" --SORT_ORDER coordinate --CREATE_INDEX true
    gatk Mutect2 -R "$d/chrM_shifted.fasta" -I "$t/sorted.bam" -L "$d/control.intervals" \
      --mitochondria-mode --max-mnp-distance 0 -O "$d/shifted_mutect2.vcf.gz"
  ' _ "$SHIFT_C" "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" "$(id -u):$(id -g)" || return 1
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
    /pgp-bin/og vcf shift-back --linear "$RAW_C" --shifted "${SHIFT_C}/shifted_mutect2.vcf.gz" \
      --out "${SHIFT_C}/merged_mutect2.vcf.gz" || return 1
  rm -f "$SHIFT_DIR"/*.img
}
echo "  Shifted-reference pass: the control region again on chrM rotated by 8000 bp..."
FILTER_IN="$RAW_C"
if shifted_pass; then
  FILTER_IN="${SHIFT_C}/merged_mutect2.vcf.gz"
else
  echo "  WARNING: the shifted-reference pass failed (see above): the control region keeps the"
  echo "  linear calls alone, and calls at its ends keep the origin_edge flag."
  rm -f "$SHIFT_DIR"/*.img
fi

echo "[4/8] Filtering variants..."
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk FilterMutectCalls \
    -R "${REF_FASTA_C}" \
    -V "$FILTER_IN" \
    --stats "${RAW_C}.stats" \
    --mitochondria-mode \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz"

echo "[5/8] Marking possible NuMTs..."
# Median autosomal depth: the depth at which half of the chr1-22 bases are
# covered at least that deep, from mosdepth's per-chromosome distribution.
MOSDEPTH_PREFIX="${SAMPLE_DIR}/mosdepth/${SAMPLE}.mosdepth"
NUMT_DEFAULT=
if [ -n "${AUTOSOMAL_COVERAGE:-}" ]; then
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (AUTOSOMAL_COVERAGE)"
elif [ -s "${MOSDEPTH_PREFIX}.summary.txt" ] && [ -s "${MOSDEPTH_PREFIX}.global.dist.txt" ]; then
  AUTOSOMAL_COVERAGE=$(awk -F'\t' '
    FNR == NR { if ($1 ~ /^chr[0-9]+$/) len[$1] = $2; next }
    ($1 in len) { at_least[$2] += len[$1] * $3 }
    END {
      for (c in len) total += len[c]
      best = 0
      if (total > 0) for (k in at_least) if (at_least[k] / total >= 0.5 && k + 0 > best) best = k + 0
      print best
    }' "${MOSDEPTH_PREFIX}.summary.txt" "${MOSDEPTH_PREFIX}.global.dist.txt")
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (from ${MOSDEPTH_PREFIX}.global.dist.txt)"
else
  AUTOSOMAL_COVERAGE=0
  NUMT_DEFAULT=--numt-default
  echo "  WARNING: no mosdepth output (run scripts/16b-mosdepth.sh first) and no AUTOSOMAL_COVERAGE:"
  echo "  the NuMT filter runs with the default depth 0 and marks nothing; the annotation table says so."
fi
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk NuMTFilterTool \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz" \
    --autosomal-coverage "$AUTOSOMAL_COVERAGE" \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"

echo "[6/8] Annotating the calls with ClinVar and gnomAD..."
CLINVAR_VCF="${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz"
GNOMAD_MITO="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
[ -s "${CLINVAR_VCF}.tbi" ] || echo "  WARNING: no ${CLINVAR_VCF}.tbi (setup.sh installs ClinVar): no ClinVar columns."
[ -s "$GNOMAD_MITO" ] || echo "  WARNING: no ${GNOMAD_MITO} (setup.sh installs it): no gnomAD columns."
if [ ! -s "${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz" ]; then
  echo "ERROR: NuMTFilterTool wrote no ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz" >&2
  exit 1
fi
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
  /pgp-bin/og mito annotate \
    --vcf "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" \
    --clinvar "$(cpath "$CLINVAR_VCF")" \
    --gnomad "$(cpath "$GNOMAD_MITO")" \
    ${NUMT_DEFAULT:+"$NUMT_DEFAULT"} \
    --out "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_annotated.tsv"

echo "[7/8] mtDNA copy number (chrM against autosomal reads, from the BAM index)..."
run_in --cpus 1 --memory 1g "$SAMTOOLS_IMAGE" \
  bash -c "samtools idxstats /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam > /genome/${SAMPLE}/mito/${SAMPLE}_idxstats.txt"
run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
  /pgp-bin/og mito copies \
    --idxstats "/genome/${SAMPLE}/mito/${SAMPLE}_idxstats.txt" \
    --out "/genome/${SAMPLE}/mito/${SAMPLE}_mtdna_copies.tsv"

echo "[8/8] Large deletions (eKLIPse, soft-clipped chrM reads)..."
# eKLIPse writes into a new eKLIPse_<id> folder under -out and -tmp on each run.
EK_DIR="${OUTPUT_DIR}/eklipse"
rm -rf "$EK_DIR"
mkdir -p "${EK_DIR}/out" "${EK_DIR}/tmp"
printf '%s\t%s\n' "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" "$SAMPLE" > "${EK_DIR}/input.tsv"
# GRCh38's chrM is the rCRS: the reference is the GenBank file the image ships beside eKLIPse.py.
# shellcheck disable=SC2016  # expanded inside the container
run_in --cpus 2 --memory 4g "$EKLIPSE_IMAGE" bash -c '
  EK=$(command -v eKLIPse.py)
  python3 "$EK" -in "$1/input.tsv" -ref "$(dirname "$EK")/data/NC_012920.1.gb" -out "$1/out" -tmp "$1/tmp" -thread 2
' _ "/genome/${SAMPLE}/mito/eklipse"
# eKLIPse stops with exit 1 on an error (set -e stops here then). Should it
# still leave no table, there is none (the old one went at the start) and the
# report says the deletions were not looked for.
CSV=$(find "${EK_DIR}/out" -name eKLIPse_deletions.csv | head -n 1)
if [ -z "$CSV" ]; then
  echo "  WARNING: eKLIPse exited 0 but wrote no eKLIPse_deletions.csv under ${EK_DIR}/out: no deletion table."
else
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$OG_IMAGE" \
    /pgp-bin/og mito deletions \
      --eklipse "/genome/${SAMPLE}/mito/eklipse/out/${CSV#"${EK_DIR}/out/"}" \
      --out "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_deletions.tsv"
fi

echo "=== Mitochondrial analysis complete ==="
echo "Raw calls: ${OUTPUT_DIR}/${SAMPLE}_chrM_mutect2.vcf.gz (linear chrM)"
if [ "$FILTER_IN" != "$RAW_C" ]; then
  echo "Shifted-reference pass: ${SHIFT_DIR}/shifted_mutect2.vcf.gz (rotated positions), merged: ${SHIFT_DIR}/merged_mutect2.vcf.gz"
fi
echo "Filtered (FilterMutectCalls, then possible_numt): ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz"
echo "ClinVar and gnomAD per allele: ${OUTPUT_DIR}/${SAMPLE}_chrM_annotated.tsv"
echo "mtDNA copy number: ${OUTPUT_DIR}/${SAMPLE}_mtdna_copies.tsv"
echo "Large deletions (eKLIPse): ${OUTPUT_DIR}/${SAMPLE}_chrM_deletions.tsv (eKLIPse's own files and plot: ${EK_DIR}/out)"
echo ""
echo "PASS variants below AF 0.95 (possible heteroplasmy; without the shifted-reference"
echo "pass, calls near the ends of the control region, chrM:1-500 and 16000-16569, are"
echo "less reliable, and length variants in the D-loop poly-C tracts m.303-315 and"
echo "m.16184-16193 are common and usually technical, so the reports count them apart):"
echo "  bcftools query -f '%POS\t%REF\t%ALT\t[%AF]\n' -i 'FILTER=\"PASS\"' ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz | awk '\$4 < 0.95'"
