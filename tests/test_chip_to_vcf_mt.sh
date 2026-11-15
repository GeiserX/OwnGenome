#!/usr/bin/env bash
# scripts/chip-to-vcf.sh: MyHeritage input and the mitochondrion.
#
# Runs the script with the real bcftools and Picard images on a made-up
# reference: 3 kb contigs 1, X and Y, and a 16,569-base mitochondrion that is
# the same sequence in the "hg19" FASTA (MT, the rCRS of human_g1k_v37) and in
# the "GRCh38" one (chrM). The chain maps 1, X and Y one to one, and chrM the
# way UCSC's hg19ToHg38 chain does: hg19's chrM is 16,571 bases, so past an
# early gap every position moves two bases (1555 to 1553 here).
#
# Checks:
#   1. tests/fixtures/chip/myheritage.csv (quoted CSV, a column header row,
#      CRLF line endings): the script exits 0, the TSV it hands bcftools has
#      ten rows and no carriage return, and the VCF has ten records, the
#      no-call as a missing genotype, every genotype a plain call, on chr1,
#      chrM, chrX and chrY.
#   2. chrM keeps its rCRS positions: m.73, m.1555 and m.3243 sit at 73, 1555
#      and 3243 with the REF of the rCRS there, and nothing sits at 1553 or
#      3241 (a lift of chrM through the chain would put m.1555 on 1553).
#   3. The same rows as a 23andMe file, its MT sites numbered on hg19's chrM
#      as 23andMe numbers them (1557 and 3245 under this chain, plus 301
#      inside the chain's gap), give the same chrM records: each MT position
#      moves to the rCRS through the chrM chain and the gap site is dropped
#      with its count printed.
#
# Needs docker. Run: bash tests/test_chip_to_vcf_mt.sh
# It replaces tests/test_chip_to_vcf_mt.py.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
SCRIPT="${REPO}/scripts/chip-to-vcf.sh"
FIXTURE="${REPO}/tests/fixtures/chip/myheritage.csv"
MT_LEN=16569
SHORT=3000
# rCRS-like bases the checks rely on: 1553 and 1555 the same base (so the
# shifted record passes Picard's REF check), 3241 and 3243 different ones.
FIXED="73=A 1553=A 1555=A 3241=G 3243=A"
FAILS=0

check() { # check DESC STATUS [DETAIL]: STATUS 0 passes
  if [ "$2" -eq 0 ]; then
    echo "PASS $1"
  else
    echo "FAIL $1  ${3:-}"
    FAILS=$((FAILS + 1))
  fi
}

command -v docker >/dev/null 2>&1 || { echo "FAIL: docker is needed to run this test"; exit 1; }
# shellcheck disable=SC1091
IMG=$(. "${REPO}/versions.env"; printf %s "${BCFTOOLS_IMAGE:-}")
[ -n "$IMG" ] || { echo "FAIL: BCFTOOLS_IMAGE is not set in versions.env"; exit 1; }

# sequence N SEED [FIXED]: N bases from the generator s = (75 s + 74) mod
# 65537, base "ACGT"[s mod 4], with POS=BASE pairs of FIXED set (1-based).
sequence() {
  awk -v n="$1" -v s="$2" -v fixed="${3:-}" 'BEGIN {
    split("A C G T", b, " ")
    for (i = 1; i <= n; i++) { s = (s * 75 + 74) % 65537; seq[i] = b[s % 4 + 1] }
    k = split(fixed, f, " ")
    for (j = 1; j <= k; j++) { split(f[j], pb, "="); seq[pb[1] + 0] = pb[2] }
    for (i = 1; i <= n; i++) printf "%s", seq[i]
    print ""
  }'
}

# write_fasta FASTA NAME SEQ [NAME SEQ]...: 60 bases a line, and its .fai.
write_fasta() {
  local fa=$1 off=0 name seq len lines
  shift
  : > "$fa"
  : > "${fa}.fai"
  while [ $# -gt 0 ]; do
    name=$1 seq=$2
    shift 2
    len=${#seq}
    lines=$(( (len + 59) / 60 ))
    printf '>%s\n' "$name" >> "$fa"
    off=$((off + ${#name} + 2))
    printf '%s\n' "$seq" | fold -w 60 >> "$fa"
    printf '%s\t%d\t%d\t60\t61\n' "$name" "$len" "$off" >> "${fa}.fai"
    off=$((off + len + lines))
  done
}

# setup_genome DIR: the two references, the dict and the chain; the chrM
# sequence is left in DIR/mt.txt.
setup_genome() {
  local g=$1 n1 nx ny mt c
  mkdir -p "$g/reference_hg19" "$g/reference" "$g/liftover"
  n1=$(sequence "$SHORT" 7)
  nx=$(sequence "$SHORT" 8)
  ny=$(sequence "$SHORT" 9)
  mt=$(sequence "$MT_LEN" 11 "$FIXED")
  printf '%s' "$mt" > "$g/mt.txt"
  write_fasta "$g/reference_hg19/human_g1k_v37.fasta" 1 "$n1" X "$nx" Y "$ny" MT "$mt"
  write_fasta "$g/reference/test38.fasta" chr1 "$n1" chrX "$nx" chrY "$ny" chrM "$mt"
  {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    for c in chr1 chrX chrY; do printf '@SQ\tSN:%s\tLN:%d\n' "$c" "$SHORT"; done
    printf '@SQ\tSN:chrM\tLN:%d\n' "$MT_LEN"
  } > "$g/reference/test38.dict"
  {
    for c in chr1 chrX chrY; do
      printf 'chain 3000 %s 3000 + 0 3000 %s 3000 + 0 3000 1\n3000\n\n' "$c" "$c"
    done
    # hg19's chrM has two more bases than the rCRS; the gap after 300 bases
    # moves every later position by two, as UCSC's chain moves m.1555 to 1553.
    printf 'chain 16569 chrM 16571 + 0 16571 chrM 16569 + 0 16569 2\n300\t2\t0\n16269\n\n'
  } | gzip -c > "$g/liftover/hg19ToHg38.over.chain.gz"
}

# run_script DIR FORMAT: runs chip-to-vcf.sh, its output in DIR/log.
run_script() {
  local rc=0
  GENOME_DIR=$1 REF_FASTA="$1/reference/test38.fasta" "$SCRIPT" S1 "$2" > "$1/log" 2>&1 || rc=$?
  sed 's/^/    /' "$1/log"
  return "$rc"
}

# records DIR: CHROM POS REF ALT GT of each record of the VCF, in DIR/records.
records() {
  docker run --rm -v "$1:/g:ro" "$IMG" bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' /g/S1/vcf/S1.vcf.gz > "$1/records"
}

# check_mt NAME DIR
check_mt() {
  local name=$1 g=$2 pos ref got
  for pos in 73 1555 3243; do
    ref=$(cut -c "$pos" "$g/mt.txt")
    got=$(awk -F'\t' -v p="$pos" '$1 == "chrM" && $2 == p { print $3 }' "$g/records")
    st=0
    [ "$got" = "$ref" ] || st=1
    check "${name}: m.${pos} stays at chrM:${pos} with the rCRS base ${ref} as REF" "$st" \
      "$(awk -F'\t' '$1 == "chrM" { printf "%s ", $2 }' "$g/records")"
  done
  st=0
  ! awk -F'\t' '$1 == "chrM" && ($2 == 1553 || $2 == 3241) { found = 1 } END { exit !found }' "$g/records" || st=1
  check "${name}: nothing moved to chrM:1553 or chrM:3241" "$st"
  got=$(awk -F'\t' '$1 == "chrM" && $2 == 1555 { print $4 " " $5 }' "$g/records")
  st=0
  [ "$got" = "G 1" ] || st=1
  check "${name}: m.1555 A>G is the alt call" "$st" "$got"
}

WORK=$(mktemp -d "${TMPDIR:-/tmp}/chip-mt.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# 1 and 2: the MyHeritage fixture
g="$WORK/myheritage"
setup_genome "$g"
mkdir -p "$g/S1/raw"
cp "$FIXTURE" "$g/S1/raw/MyHeritage_raw_dna_data.csv"
st=0
grep -q $'\r$' "$FIXTURE" || st=1
check "the MyHeritage fixture has CRLF line endings" "$st"
echo "--- myheritage: chip-to-vcf.sh S1 auto"
run_script "$g" auto
rc=$?
st=0
grep -q "Auto-detected format: myheritage" "$g/log" || st=1
check "myheritage: auto-detected" "$st"
check "myheritage: chip-to-vcf.sh exits 0 (exit ${rc})" "$rc"
tsv="$g/S1/raw/S1_tsv2vcf.tsv"
st=0
[ -f "$tsv" ] && [ "$(wc -l < "$tsv" | tr -d ' ')" -eq 10 ] && ! grep -q $'\r' "$tsv" || st=1
check "myheritage: the converted TSV has ten rows and no carriage return" "$st" \
  "$(head -c 200 "$tsv" 2>/dev/null | od -c | head -5)"
if [ "$rc" -eq 0 ] && records "$g"; then
  st=0
  [ "$(wc -l < "$g/records" | tr -d ' ')" -eq 10 ] || st=1
  check "myheritage: 10 records" "$st" "$(cat "$g/records")"
  st=0
  ! cut -f5 "$g/records" | grep -qvE '^(\.|[0-9]+)([/|](\.|[0-9]+))?$' || st=1
  check "myheritage: every genotype is a plain call, no carriage return" "$st" \
    "$(cut -f5 "$g/records" | sort -u | od -c | head -5)"
  got=$(awk -F'\t' '$1 == "chr1" && $2 == 404 { print $5 }' "$g/records")
  st=0
  [ "$got" = "./." ] || [ "$got" = "." ] || st=1
  check "myheritage: the no-call row is a missing genotype" "$st" "$got"
  got=$(cut -f1 "$g/records" | sort -u | tr '\n' ' ')
  st=0
  [ "$got" = "chr1 chrM chrX chrY " ] || st=1
  check "myheritage: chromosomes chr1 chrM chrX chrY" "$st" "$got"
  check_mt myheritage "$g"
else
  check "myheritage: the VCF reads back" 1
fi

# 3: the same rows as a 23andMe file
g="$WORK/23andme"
setup_genome "$g"
mkdir -p "$g/S1/raw"
{
  printf '# rsid\tchromosome\tposition\tgenotype\n'
  # rCRS to hg19 chrM numbering under the test chain: m.1555 is 1557 there
  # and m.3243 is 3245 (two bases more past the gap at 301-302).
  tr -d '\r"' < "$FIXTURE" | awk -F',' -v OFS='\t' '
    /^#/ || NF == 0 || $1 == "RSID" { next }
    $2 == "MT" && $3 == 1555 { $3 = 1557 }
    $2 == "MT" && $3 == 3243 { $3 = 3245 }
    { print $1, $2, $3, $4 }'
  # a site inside the chain's gap: no rCRS position, dropped
  printf 'rs9000099\tMT\t301\tA\n'
} > "$g/S1/raw/S1_raw.txt"
echo "--- 23andme: chip-to-vcf.sh S1 23andme"
run_script "$g" 23andme
rc=$?
check "23andme: chip-to-vcf.sh exits 0 (exit ${rc})" "$rc"
st=0
grep -q "3 site(s) moved from hg19 chrM to rCRS positions, 1 in a chain gap dropped" "$g/log" || st=1
check "23andme: the MT move and the gap drop are printed" "$st"
if [ "$rc" -eq 0 ] && records "$g"; then
  check_mt 23andme "$g"
  got=$(awk -F'\t' '$1 == "chrM" { print $2 }' "$g/records" | sort -n | tr '\n' ' ')
  st=0
  [ "$got" = "73 1555 3244 " ] || st=1
  check "23andme: chrM holds 73, 1555 and 3243 only (nothing left at 301, 1557 or 3245)" "$st" "$got"
else
  check "23andme: the VCF reads back" 1
fi

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all checks passed"
