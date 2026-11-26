#!/usr/bin/env bash
# Compare the copy-number engine (og cn) against frozen results from the
# per-gene tools at 605645249e6762c52a84ac58a28968f269b16536. The input is
# HG002's public 35x BWA-MEM GRCh38 no-ALT reads at HBA, RH, and LPA,
# on the fixture header at full depth and after sampling one eighth.
# The RhD phenotype, genotype, confidence, and evidence must match.
# The LPA table must match after normalizing the interval wording and
# bounds. Its current interval must contain the old interval.
# The HBA regions and result table must both match byte for byte.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR" T="$CASE_TMP"
GOLDEN="$REPO/tests/parity/golden/cn-engine"
SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
FIX="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$1" "${@:2}"; }
retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }

# --- HG002's reads at the three loci, on the fixture's header --------------------------
fetch_reads() {
  rm -f "${T}/win.bam" "${T}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${E2E_WORK}:${E2E_WORK}" -w "$T" "$SAMTOOLS_IMAGE" \
    samtools view -b -o "${T}/win.bam" "$SRC" chr16:95001-235000 chr1:25150000-25450000 chr6:160530000-160670000
}
check "HG002's reads at the three loci download" retry fetch_reads
{ dk "$SAMTOOLS_IMAGE" samtools view -H "$FIX" | grep -v -e '^@RG' -e '^@PG'
  dk "$SAMTOOLS_IMAGE" samtools view -H "${T}/win.bam" | grep '^@RG'
  dk "$SAMTOOLS_IMAGE" samtools view "${T}/win.bam" \
    | awk -F'\t' '$7 == "=" || $7 == "*" || $7 == "chr1" || $7 == "chr6" || $7 == "chr16"'
} > "${T}/cn.sam"
dk "$SAMTOOLS_IMAGE" samtools sort -o "${T}/full.bam" "${T}/cn.sam"
dk "$SAMTOOLS_IMAGE" samtools view -b -s 42.125 -o "${T}/thin.bam" "${T}/full.bam"
rm -f "${T}/cn.sam" "${T}/win.bam"
check_ge "HG002 reads at the three loci" "$(dk "$SAMTOOLS_IMAGE" samtools view -c "${T}/full.bam" 2>/dev/null || echo 0)" 50000

# dkb BIN IMAGE CMD...: dk with a tree's bin/ at its own path (this checkout is outside the work folder).
dkb() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${E2E_WORK}:${E2E_WORK}" -v "$1:$1:ro" -w "$T" "$2" "${@:3}"; }
HBA_BED="$T/hba.bed"
dkb "$REPO/bin" "$PARASCOPY_IMAGE" "$REPO/bin/og" cn bed --locus HBA > "$HBA_BED"
check "the HBA regions match the golden" cmp -s "$GOLDEN/regions/hba.bed" "$HBA_BED"

# hba_block OUT BAM: step 35 uses all MAPQ and MAPQ >= 1 before og cn.
hba_block() {
  local out=$1 bam=$2
  dk "$PARASCOPY_IMAGE" samtools bedcov "$HBA_BED" "$bam" > "$T/hba-all.tsv"
  dk "$PARASCOPY_IMAGE" samtools bedcov -Q 1 "$HBA_BED" "$bam" > "$T/hba-q1.tsv"
  dkb "$REPO/bin" "$PARASCOPY_IMAGE" "$REPO/bin/og" cn call --locus HBA \
    --bedcov "$T/hba-all.tsv" --bedcov-mapq1 "$T/hba-q1.tsv" --table "$out"
}

check_golden() {
  local description=$1 expected=$2 actual=$3
  check "$description" cmp -s "$expected" "$actual"
  cmp -s "$expected" "$actual" || diff -u "$expected" "$actual" | sed -n '1,20p'
}

for V in full thin; do
  N="${SAMPLE}cn${V}"
  rm -rf "${G:?}/${N}"; mkdir -p "${G}/${N}/aligned" "${G}/${N}/vcf"
  cp "${T}/${V}.bam" "${G}/${N}/aligned/${N}_sorted.bam"
  dk "$SAMTOOLS_IMAGE" samtools index "${G}/${N}/aligned/${N}_sorted.bam"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${N}/vcf/${N}.vcf.gz"
  cp "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" "${G}/${N}/vcf/${N}.vcf.gz.tbi"
  NEW="${T}/new-${V}" WANT="$GOLDEN/$V"
  rm -rf "$NEW"; mkdir -p "$NEW"

  run_step 46-lpa.sh "$N"; check_step_exit 46-lpa.sh
  cp "${G}/${N}/lpa/${N}_lpa_kiv2.tsv" "${NEW}/lpa.tsv"
  run_step 41-bloodtype.sh "$N"; check_step_exit 41-bloodtype.sh
  cp "${G}/${N}/bloodtype/${N}_bloodtype.tsv" "${NEW}/bloodtype.tsv"
  cp "${G}/${N}/bloodtype/${N}_bloodtype_evidence.tsv" "${NEW}/evidence.tsv"
  hba_block "${NEW}/hba.tsv" "${G}/${N}/aligned/${N}_sorted.bam"
  for f in lpa bloodtype evidence hba; do check "${V}: this tree wrote its ${f} table" test -s "${NEW}/${f}.tsv"; done
  cat "${NEW}/lpa.tsv" "${NEW}/hba.tsv"; grep '^RhD' "${NEW}/bloodtype.tsv"
  echo "cn-engine: ${V} LPA $(awk -F'\t' '$1 == "kiv2_copies" {print $2}' "${NEW}/lpa.tsv")," \
    "HBA $(awk -F'\t' '$1 == "alpha_globin_genes" {print $2}' "${NEW}/hba.tsv")," \
    "RhD $(awk -F'\t' '$1 == "RhD" {print $3}' "${NEW}/bloodtype.tsv")" >> "$E2E_NOTES"

  # Keep the header and the RhD rows. The other systems have their own tests.
  for f in bloodtype evidence; do
    awk -F'\t' -v OFS='\t' -v f="$f" 'NR == 1 || $1 == "RhD" {if (f == "bloodtype") NF = 4; print}' \
      "$NEW/$f.tsv" > "$NEW/$f.rhd.tsv"
    check_ge "$V: RhD rows in the $f table" "$(($(wc -l < "$NEW/$f.rhd.tsv") - 1))" 1
  done
  # Step 46's message names what its interval covers since the old tree, and
  # the interval adds the error measured against assembled genomes (bin/cn_loci.tsv
  # interval cv:SYS, CHANGELOG). The count, the class and every other row
  # stay the old tree's: both tables lose the interval's bounds (the rows
  # and the message's), are compared byte for byte, and this tree's interval
  # must hold the old one (the term only widens it).
  sed -E -e '/^ci95_(low|high)\t/d' \
    -e '/^message\t/s/ in total \(95% interval [0-9.]+ to [0-9.]+\)/ in total over both copies (95% interval)/' \
    "${WANT}/lpa.tsv" > "${NEW}/lpa.golden.now.tsv"
  sed -E -e '/^ci95_(low|high)\t/d' \
    -e '/^message\t/s/ \(95% interval [0-9.]+ to [0-9.]+, read sampling plus the error measured against assembled genomes\)/ (95% interval)/' \
    "${NEW}/lpa.tsv" > "${NEW}/lpa.now.tsv"
  lpa_ci() { awk -F'\t' '$1 == "ci95_low" {l = $2} $1 == "ci95_high" {h = $2} END {print l, h}' "$1"; }
  read -r ol oh <<< "$(lpa_ci "${WANT}/lpa.tsv")"
  read -r nl nh <<< "$(lpa_ci "${NEW}/lpa.tsv")"
  check "${V}: this tree's interval (${nl} to ${nh}) holds the old tree's (${ol} to ${oh})" \
    awk -v ol="$ol" -v oh="$oh" -v nl="$nl" -v nh="$nh" \
    'BEGIN {if (ol == "." || nl == ".") exit !(ol == nl && oh == nh); exit !(nl + 0 <= ol + 0 && nh + 0 >= oh + 0)}'
  check_golden "$V: normalized LPA matches the golden" "$NEW/lpa.golden.now.tsv" "$NEW/lpa.now.tsv"
  for f in bloodtype.rhd evidence.rhd hba; do
    check_golden "$V: the $f table matches the golden" "$WANT/$f.tsv" "$NEW/$f.tsv"
  done
done
check "control: the comparison sees a difference (full against thin)" bash -c '! cmp -s "$1" "$2"' _ \
  "${T}/new-full/lpa.tsv" "${T}/new-thin/lpa.tsv"

rm -rf "${G:?}/${SAMPLE}cnfull" "${G:?}/${SAMPLE}cnthin" "${T:?}"/*.bam
finish
