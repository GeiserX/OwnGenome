#!/usr/bin/env bash
# 22-survivor-merge.sh — SV consensus: the calls two or more SV callers agree on
# Usage: ./scripts/22-survivor-merge.sh <sample_name>
#
# SURVIVOR merge (SURVIVOR_IMAGE) pairs the calls of the callers below when
# both breakpoints lie within 1,000 bp of each other, the SV type and strands
# agree and the event is at least 50 bp long (`SURVIVOR merge LIST 1000 2 1 1
# 0 50`), and keeps the events at least two callers support. Each record
# carries SUPP (how many callers) and SUPP_VEC (which: one digit per caller,
# in the order of sv_files.txt).
#
# SURVIVOR's 1 kb rule alone pairs two deletions of 100 bp that start 900 bp
# apart and share no base, and on some inputs a deletion with a duplication at
# the same place. og sv consensus-check then reads each merged record's calls
# (FORMAT CO, TY, ST and LN) and keeps it only when they are one event: the
# same type (INV, BND and TRA count as one, SURVIVOR's rule), the same strands
# and a reciprocal overlap of at least 0.5; insertions start within 1 kb and
# the shorter is at least half the longer. A record that fails is split by
# event, each with its POS, END, SVLEN, CIPOS and CIEND from its own calls,
# and a part with fewer than 2 callers is dropped.
#
# Inputs, each used when it exists (two or more are needed): Manta (step 04),
# Delly (19), CNVpytor (18), TIDDIT (04a) and Sniffles2 (04c). GRIDSS (04b) is
# left out: it reports every event as a pair of breakends (SVTYPE=BND), and
# nothing here turns a pair into one DEL, DUP or INV. SURVIVOR 1.0.7 counts BND
# as the same type as INV and TRA, so its breakends would pair with the other
# callers' inversions and translocations and count as votes for an event GRIDSS
# never typed.
#
# A caller votes only with the calls it makes: og sv survivor-prep drops a
# record whose FILTER or sample FT is not PASS, or whose genotype is reference
# only (0/0, 0|0, 0, 0/.), and for CNVpytor a call that fails step 05's quality
# rule (e-val < 0.01, q0 < 0.5). A no-call (./.) and a sites-only record vote.
#
# SUPP counts callers that agree, not independent proof: Manta and Delly both
# read split and paired reads, so they share alignment errors. AnnotSV (step 05)
# annotates the callers' own calls; the reports join SUPP onto them.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image SURVIVOR_IMAGE BCFTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${S}/sv_merged"
IN_DIR="${OUTDIR}/inputs"
OUT="${OUTDIR}/${SAMPLE}_sv_consensus.vcf.gz"
mkdir -p "$IN_DIR"

echo "============================================"
echo "  Step 22: SV consensus (SURVIVOR merge)"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUT}"
echo "============================================"
echo ""

# CNVpytor's own table, when step 18 left no VCF: one DEL or DUP record per
# CNV, with a sample column, so SURVIVOR reads it as it reads the others. The
# e-value (column 5), q0 (9) and pN (10) go into INFO under the names CNVpytor's
# VCF uses, so its quality rule reads both inputs the same.
CNV_TXT="${S}/cnvpytor/${SAMPLE}_cnvs.txt"
CNV_VCF="${S}/cnvpytor/${SAMPLE}_cnvs.vcf.gz"
CNV_MAKE=false
if [ -f "$CNV_TXT" ]; then
  if [ ! -f "$CNV_VCF" ]; then
    CNV_MAKE=true
  else
    # A VCF an earlier version of this step made from the table has no
    # e-value: none of its calls would vote, so it is made again. Step 18's
    # own VCF declares pytorP1 and is kept.
    CNV_HEAD=$(run_in "$BCFTOOLS_IMAGE" bcftools view -h "$(cpath "$CNV_VCF")") || CNV_HEAD=
    case "$CNV_HEAD" in
      *'##INFO=<ID=pytorP1,'*) ;;
      *)
        echo "CNVpytor: ${CNV_VCF#"$S"/} has no e-value (pytorP1); making it again from the table."
        CNV_MAKE=true ;;
    esac
  fi
fi
if [ "$CNV_MAKE" = true ]; then
  echo "CNVpytor: converting ${CNV_TXT} to a VCF..."
  {
    echo "##fileformat=VCFv4.2"
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">"
    echo "##INFO=<ID=pytorP1,Number=1,Type=Float,Description=\"e-value (t-test, whole region), column 5 of the CNVpytor table\">"
    echo "##INFO=<ID=pytorQ0,Number=1,Type=Float,Description=\"Fraction of reads with mapping quality 0, column 9 of the CNVpytor table\">"
    echo "##INFO=<ID=pytorPN,Number=1,Type=Float,Description=\"Fraction of N bases, column 10 of the CNVpytor table\">"
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REF_FASTA}.fai"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$SAMPLE"
    awk -F'\t' '$1 == "deletion" || $1 == "duplication" {
      split($2, a, ":"); split(a[2], b, "-")
      t = ($1 == "duplication") ? "DUP" : "DEL"
      len = (t == "DEL") ? -(b[2] - b[1]) : b[2] - b[1]
      printf "%s\t%s\t.\tN\t<%s>\t.\tPASS\tSVTYPE=%s;END=%s;SVLEN=%d;pytorP1=%s;pytorQ0=%s;pytorPN=%s\tGT\t./.\n", a[1], b[1], t, t, b[2], len, $5, $9, $10
    }' "$CNV_TXT"
  } > "${IN_DIR}/cnvpytor_txt.vcf"
  run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
    'bcftools sort "$1" -Oz -o "$2.tmp" && mv "$2.tmp" "$2" && bcftools index -f -t "$2"' \
    _ "$(cpath "${IN_DIR}/cnvpytor_txt.vcf")" "$(cpath "$CNV_VCF")"
  rm -f "${IN_DIR}/cnvpytor_txt.vcf"
fi

# name<TAB>VCF, in the order SUPP_VEC reports them.
CALLERS=(
  "manta	${S}/manta/results/variants/diploidSV.vcf.gz"
  "delly	${S}/delly/${SAMPLE}_sv.vcf.gz"
  "cnvpytor	${CNV_VCF}"
  "tiddit	${S}/sv_tiddit/${SAMPLE}_sv.vcf.gz"
  "sniffles2	${S}/sv_sniffles/${SAMPLE}_sv.vcf.gz"
)
USED=()
rm -f "${IN_DIR}"/*.vcf "${IN_DIR}/sv_files.txt"
for c in "${CALLERS[@]}"; do
  name=${c%%	*} vcf=${c#*	}
  if [ ! -f "$vcf" ]; then
    echo "  [--] ${name}: no ${vcf#"$S"/}"
    continue
  fi
  # The calls that vote, as plain text, which SURVIVOR reads, with one
  # sample column named after the caller: SURVIVOR names its output columns
  # after the input samples, and three columns all named after the sample
  # would be one name three times. og prints what it left out.
  run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$BCFTOOLS_IMAGE" \
    /pgp-bin/og sv survivor-prep --caller "$name" "$(cpath "$vcf")" "$(cpath "${IN_DIR}/${name}.vcf")"
  n=$(grep -vc '^#' "${IN_DIR}/${name}.vcf" || true)
  echo "  [OK] ${name}: ${n} calls that vote (${vcf#"$S"/})"
  echo "$(cpath "${IN_DIR}/${name}.vcf")" >> "${IN_DIR}/sv_files.txt"
  USED+=("$name")
done
if [ -f "${S}/sv_gridss/${SAMPLE}_gridss.vcf.gz" ]; then
  echo "  [--] gridss: left out of the SV consensus: its breakend pairs (BND) are not turned into DEL, DUP or INV events"
fi
echo ""

if [ "${#USED[@]}" -lt 2 ]; then
  echo "ERROR: the consensus needs at least 2 SV callers; found ${#USED[@]} (${USED[*]:-none})." >&2
  echo "  Run at least two of: steps 04 (Manta), 19 (Delly), 18 (CNVpytor), 04a (TIDDIT), 04c (Sniffles2)." >&2
  exit 1
fi

echo "[1/3] SURVIVOR merge of ${USED[*]}: breakpoints within 1,000 bp, same type and strands, >= 50 bp, 2+ callers..."
MERGED="${OUTDIR}/${SAMPLE}_sv_survivor.vcf"
CHECKED="${OUTDIR}/${SAMPLE}_sv_merged.vcf"
rm -f "$MERGED" "$CHECKED"
run_in --cpus 1 --memory 4g "$SURVIVOR_IMAGE" \
  SURVIVOR merge "$(cpath "${IN_DIR}/sv_files.txt")" 1000 2 1 1 0 50 "$(cpath "$MERGED")"
# SURVIVOR exits 0 when it cannot open an input, so its output is the check.
if ! wrote_vcf "$MERGED"; then
  echo "ERROR: SURVIVOR wrote no VCF (${MERGED})." >&2
  exit 1
fi

echo "[2/3] Checking that each merged record is one event (og sv consensus-check)..."
run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$BCFTOOLS_IMAGE" \
  /pgp-bin/og sv consensus-check "$(cpath "$MERGED")" "$(cpath "$CHECKED")"
rm -f "$MERGED"

echo "[3/3] Sorting, compressing and indexing..."
rm -f "${OUT}.tbi"
run_in --cpus 2 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
  'bcftools sort "$1" -Oz -o "$2.tmp" && mv "$2.tmp" "$2" && bcftools index -f -t "$2"' \
  _ "$(cpath "$CHECKED")" "$(cpath "$OUT")"
rm -f "$CHECKED"

COUNT=$(run_in "$BCFTOOLS_IMAGE" bcftools view -H "$(cpath "$OUT")" | wc -l | tr -d ' ')
echo ""
echo "============================================"
echo "  SV consensus complete: ${SAMPLE}"
echo "  Callers (SUPP_VEC order): ${USED[*]}"
echo "  Consensus SVs (2+ callers): ${COUNT}"
echo "  Output: ${OUT}"
echo "============================================"
echo "SUPP is caller agreement, not validation: Manta and Delly read the same split and"
echo "  paired reads. Annotate the callers' own calls (step 05), not this file."
echo "By support (SUPP_VEC: one digit per caller above):"
run_in "$BCFTOOLS_IMAGE" bcftools query -f '%INFO/SVTYPE\t%INFO/SUPP_VEC\n' "$(cpath "$OUT")" \
  | sort | uniq -c | sort -rn | head -20
