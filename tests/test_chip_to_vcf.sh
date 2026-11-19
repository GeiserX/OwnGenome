#!/usr/bin/env bash
# scripts/chip-to-vcf.sh on the two ten-row vendor fixtures in tests/fixtures/chip/,
# with the real bcftools and Picard images, on a made-up 3 kb "hg19" and
# "GRCh38" (the same sequence under both naming styles) and a one-to-one chain.
#
# Checks, for each fixture:
#   1. the script exits 0 and writes vcf/<sample>.vcf.gz;
#   2. ten records, the no-call row as a missing genotype (./. or .);
#   3. the records sit on chr1, chrX, chrY and chrM (AncestryDNA codes 23-26
#      mapped to X, Y, X and MT, then renamed for GRCh38);
#   4. AncestryDNA only: every genotype is diploid (its two allele columns are
#      joined). The script before this fix read the file as 23andMe, so it
#      took the column header as a row, the first allele as the genotype and
#      the numbers 23-26 as chromosome names;
#   5. every called site has the alleles and genotype of its fixture row
#      (%TGT and %GT), hom-ALT rows included: a converter that writes a
#      homozygous ALT call as 0/0, the plink failure the script exists to
#      avoid, fails here;
#   6. the build check says these few heterozygous calls are too few to check.
# Then a cut AncestryDNA file: the script stops with the line number of the
# damaged row and writes no VCF.
# Then the genome build, which the script checks instead of trusting the
# vendor. Each refused file exits non-zero, names the reason and writes no VCF:
#   - the 23andMe fixture under a header that says build 36;
#   - the AncestryDNA fixture (GRCh37 positions) under a header that says
#     build 38;
#   - the MyHeritage fixture under a header that says GRCh38;
#   - 250 heterozygous rows under a build 37 header where neither allele is
#     the reference base, as positions from another build give.
#   Controls: the 23andMe fixture under a line that names GRCh38 and GRCh37
#   (a file already lifted) converts; the same 250 rows with one allele the
#   reference base, in a 6,000-row file, convert, and the log says the
#   positions fit GRCh37.
# Needs docker.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=../versions.env
. "${REPO}/versions.env"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

command -v docker >/dev/null 2>&1 || { echo "FAIL: docker is needed to run this test"; exit 1; }

# write_fasta FILE NAME...: one 3,000-base contig per NAME, the same sequence
# for each NAME at the same index, plus the .fai.
write_fasta() {
  local out=$1
  shift
  awk -v names="$*" -v fa="$out" -v fai="${out}.fai" 'BEGIN {
    n = split(names, name, " ")
    split("A C G T", base, " ")
    s = 7
    off = 0
    for (c = 1; c <= n; c++) {
      hdr = ">" name[c]
      print hdr > fa
      off += length(hdr) + 1
      seq = ""
      for (i = 1; i <= 3000; i++) { s = (s * 75 + 74) % 65537; seq = seq base[s % 4 + 1] }
      lines = 0
      for (i = 1; i <= 3000; i += 60) { print substr(seq, i, 60) > fa; lines++ }
      printf "%s\t3000\t%d\t60\t61\n", name[c], off > fai
      off += 3000 + lines
    }
  }'
}

# setup_genome DIR: the made-up references, dictionary and chain.
setup_genome() {
  local g=$1 c
  mkdir -p "${g}/reference_hg19" "${g}/reference" "${g}/liftover"
  write_fasta "${g}/reference_hg19/human_g1k_v37.fasta" 1 X Y MT
  write_fasta "${g}/reference/test38.fasta" chr1 chrX chrY chrM
  {
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    for c in chr1 chrX chrY chrM; do printf '@SQ\tSN:%s\tLN:3000\n' "$c"; done
  } > "${g}/reference/test38.dict"
  # After the rename the hg19 side is chr1, chrX, chrY, chrM too.
  for c in chr1 chrX chrY chrM; do
    printf 'chain 3000 %s 3000 + 0 3000 %s 3000 + 0 3000 1\n3000\n\n' "$c" "$c"
  done | gzip -c > "${g}/liftover/hg19ToHg38.over.chain.gz"
}

bcf() { docker run --rm -v "${1}:/g:ro" "$BCFTOOLS_IMAGE" bcftools "${@:2}"; }

# run_fixture NAME FIXTURE FORMAT
run_fixture() {
  local name=$1 fixture=$2 format=$3 g="${WORK}/${1}" rc=0 vcf
  setup_genome "$g"
  mkdir -p "${g}/S1/raw"
  cp "${REPO}/tests/fixtures/chip/${fixture}" "${g}/S1/raw/S1_raw.txt"
  echo "--- ${name}: chip-to-vcf.sh S1 ${format}"
  GENOME_DIR="$g" REF_FASTA="${g}/reference/test38.fasta" \
    "${REPO}/scripts/chip-to-vcf.sh" S1 "$format" > "${WORK}/${name}.log" 2>&1 || rc=$?
  sed 's/^/    /' "${WORK}/${name}.log"
  if [ "$rc" -ne 0 ]; then fail "${name}: chip-to-vcf.sh exited ${rc}"; return; fi
  pass "${name}: chip-to-vcf.sh exited 0"
  vcf="${g}/S1/vcf/S1.vcf.gz"
  if [ ! -s "$vcf" ]; then fail "${name}: no ${vcf}"; return; fi
  bcf "$g" query -f '%CHROM\t%POS\t[%GT]\n' /g/S1/vcf/S1.vcf.gz > "${WORK}/${name}.tsv"
  sed 's/^/    /' "${WORK}/${name}.tsv"
  local n chroms
  n=$(grep -c . "${WORK}/${name}.tsv" || true)
  if [ "$n" -eq 10 ]; then pass "${name}: 10 records"; else fail "${name}: ${n} records, want 10"; fi
  local nocall
  nocall=$(awk -F'\t' '$1 == "chr1" && $2 == 404 {print $3}' "${WORK}/${name}.tsv")
  case "$nocall" in
    ./.|.) pass "${name}: the no-call row is a missing genotype (${nocall})" ;;
    *) fail "${name}: the no-call row at chr1:404 is '${nocall}', want a missing genotype" ;;
  esac
  chroms=$(cut -f1 "${WORK}/${name}.tsv" | sort -u | paste -sd' ' -)
  if [ "$chroms" = "chr1 chrM chrX chrY" ]; then pass "${name}: chromosomes ${chroms}"
  else fail "${name}: chromosomes '${chroms}', want 'chr1 chrM chrX chrY'"; fi
  if [ "$name" = ancestrydna ]; then
    local haploid
    haploid=$(awk -F'\t' '$3 !~ /[\/|]/' "${WORK}/${name}.tsv" | grep -c . || true)
    if [ "$haploid" -eq 0 ]; then pass "${name}: every genotype is diploid"
    else fail "${name}: ${haploid} genotype(s) not diploid"; fi
  fi
  # Per-site truth: the alleles and the genotype of every called fixture row.
  # 23andMe reports one allele on X, Y and MT for this (male) fixture.
  local want got
  if [ "$name" = ancestrydna ]; then
    want="chr1:101 A/G 1/0
chr1:202 C/C 1/1
chr1:303 G/T 1/2
chr1:1010 T/T 1/1
chrX:505 A/A 0/0
chrX:606 C/T 1/2
chrX:808 T/T 1/1
chrY:707 G/G 1/1
chrM:909 A/G 0/1"
  else
    want="chr1:101 A/G 1/0
chr1:202 C/C 1/1
chr1:303 G/T 1/2
chr1:1010 T/T 1/1
chrX:505 A/A 0/0
chrX:606 C/T 1/2
chrX:808 T 1
chrY:707 G 1
chrM:909 A 0"
  fi
  got=$(bcf "$g" query -e 'GT="mis"' -f '%CHROM:%POS [%TGT] [%GT]\n' /g/S1/vcf/S1.vcf.gz)
  if [ "$got" = "$want" ]; then pass "${name}: every called site has its fixture alleles and genotype"
  else
    fail "${name}: the called sites differ from the fixture rows (want, then got)"
    printf '%s\n' "$want" "---" "$got" | sed 's/^/    /'
  fi
  # AncestryDNA's MT row is a diploid A/G: one heterozygous call more.
  local hets=3
  [ "$name" = ancestrydna ] && hets=4
  if grep -q "Build check: ${hets} heterozygous call(s), fewer than 200" "${WORK}/${name}.log"; then
    pass "${name}: the build check says ${hets} heterozygous calls are too few"
  else fail "${name}: no 'too few to check the build' line for ${hets} heterozygous calls"; fi
}

run_fixture ancestrydna ancestrydna.txt auto
run_fixture 23andme 23andme.txt auto

# A cut AncestryDNA file (its last row lost the second allele) stops with the
# line number and writes no VCF. The script before this check skipped the row
# and converted the rest.
cut_g="${WORK}/ancestrydna-cut"
setup_genome "$cut_g"
mkdir -p "${cut_g}/S1/raw"
{ cat "${REPO}/tests/fixtures/chip/ancestrydna.txt"; printf 'rs999\t1\t500\tA\n'; } > "${cut_g}/S1/raw/S1_raw.txt"
cut_line=$(wc -l < "${cut_g}/S1/raw/S1_raw.txt" | tr -d ' ')
echo "--- ancestrydna-cut: chip-to-vcf.sh S1 ancestrydna (line ${cut_line} has four columns)"
cut_rc=0
GENOME_DIR="$cut_g" REF_FASTA="${cut_g}/reference/test38.fasta" \
  "${REPO}/scripts/chip-to-vcf.sh" S1 ancestrydna > "${WORK}/ancestrydna-cut.log" 2>&1 || cut_rc=$?
sed 's/^/    /' "${WORK}/ancestrydna-cut.log"
if [ "$cut_rc" -ne 0 ]; then pass "ancestrydna-cut: chip-to-vcf.sh exited ${cut_rc}"
else fail "ancestrydna-cut: chip-to-vcf.sh exited 0 on a row with four columns"; fi
if grep -q "line ${cut_line}: want five tab-separated columns" "${WORK}/ancestrydna-cut.log"; then
  pass "ancestrydna-cut: the error names line ${cut_line}"
else fail "ancestrydna-cut: no error naming line ${cut_line}"; fi
if [ ! -e "${cut_g}/S1/vcf/S1.vcf.gz" ]; then pass "ancestrydna-cut: no VCF"
else fail "ancestrydna-cut: wrote a VCF from a cut file"; fi

# --- The genome build ---------------------------------------------------------

# convert NAME RAWFILE FORMAT: the raw file (on stdin, by redirection so rc
# stays in this shell) through chip-to-vcf.sh as sample S1; sets rc and leaves the log in WORK/NAME.log and the genome in
# WORK/NAME.
convert() {
  local name=$1 g="${WORK}/${1}"
  setup_genome "$g"
  mkdir -p "${g}/S1/raw"
  cat > "${g}/S1/raw/${2}"
  echo "--- ${name}: chip-to-vcf.sh S1 ${3}"
  rc=0
  GENOME_DIR="$g" REF_FASTA="${g}/reference/test38.fasta" \
    "${REPO}/scripts/chip-to-vcf.sh" S1 "$3" > "${WORK}/${name}.log" 2>&1 || rc=$?
  sed 's/^/    /' "${WORK}/${name}.log"
}

# refused NAME WHY PATTERN: the last convert exited non-zero, its log matches
# PATTERN (grep -E) and it wrote no VCF.
refused() {
  local name=$1 g="${WORK}/${1}"
  if [ "$rc" -ne 0 ]; then pass "${name}: chip-to-vcf.sh exited ${rc} on ${2}"
  else fail "${name}: chip-to-vcf.sh exited 0 on ${2}"; fi
  if grep -qE "$3" "${WORK}/${name}.log"; then pass "${name}: the error says why"
  else fail "${name}: no line matching '${3}'"; fi
  if [ ! -e "${g}/S1/vcf/S1.vcf.gz" ]; then pass "${name}: no VCF"
  else fail "${name}: wrote a VCF"; fi
}

# The fixtures with a header line added after their first line that names
# another build: refused before anything is converted.
labelled() { awk -v l="$1" '{ print } NR == 1 { print l }' "${REPO}/tests/fixtures/chip/${2}"; }

convert 23andme-build36 S1_raw.txt auto \
  < <(labelled '# We are using reference human assembly build 36 (also known as NCBI36).' 23andme.txt)
refused 23andme-build36 "a header naming build 36" 'line 2 names genome build 36: # We are using reference human assembly build 36'

convert ancestrydna-build38 S1_raw.txt auto \
  < <(labelled '#Columns two and three contain the chromosome and basepair position of the SNP using human reference build 38 coordinates.' ancestrydna.txt)
refused ancestrydna-build38 "GRCh37 positions under a header naming build 38" 'line 2 names genome build 38: '

convert myheritage-grch38 MyHeritage_raw_dna_data.csv auto < <(labelled '# Positions are on GRCh38.' myheritage.csv)
refused myheritage-grch38 "a header naming GRCh38" 'line 2 names genome build 38: # Positions are on GRCh38'

# A file already lifted to GRCh37, as the error tells the user to do, may name
# its old build and GRCh37 on one line: it converts.
convert 23andme-lifted S1_raw.txt auto \
  < <(labelled '# Genotypes lifted from GRCh38 to GRCh37 by a liftover tool.' 23andme.txt)
if [ "$rc" -eq 0 ] && [ -s "${WORK}/23andme-lifted/S1/vcf/S1.vcf.gz" ]; then
  pass "23andme-lifted: a header naming GRCh38 and GRCh37 on one line converts"
else fail "23andme-lifted: chip-to-vcf.sh exited ${rc} on a header naming GRCh38 and GRCh37 on one line"; fi

# het_rows MODE: a 23andMe file whose header names build 37, with 250
# heterozygous rows on chr 1 (positions 1001-1250) of the made-up reference.
# MODE miss: neither allele is the reference base, as positions from another
# build give. MODE hit: one allele is the reference base, and every other
# position of chr 1 and chr X is a hom-REF row: 6,000 rows, about 130 KB, more
# than a pipe buffer holds, so a header read that stops early cannot kill the
# script with SIGPIPE on a full-size download.
het_rows() {
  local fa="${WORK}/het-ref/reference_hg19/human_g1k_v37.fasta"
  [ -s "$fa" ] || setup_genome "${WORK}/het-ref"
  printf '# This data file generated by 23andMe at: synthetic test fixture\n'
  printf '# We are using reference human assembly build 37 (also known as Annotation Release 104).\n'
  printf '# rsid\tchromosome\tposition\tgenotype\n'
  awk -v mode="$1" '/^>/ { c = substr($1, 2); next } c == "1" || c == "X" { s[c] = s[c] $0 }
    END {
      for (p = 1001; p <= 1250; p++) {
        r = substr(s["1"], p, 1)
        o = (r == "A" ? "CG" : (r == "C" ? "AT" : (r == "G" ? "AT" : "CG")))
        printf "rs%d\t1\t%d\t%s\n", 8000000 + p, p, (mode == "miss" ? o : r substr(o, 1, 1))
      }
      if (mode == "miss") exit
      for (p = 1; p <= 3000; p++) {
        if (p < 1001 || p > 1250) printf "rs%d\t1\t%d\t%s\n", 8000000 + p, p, substr(s["1"], p, 1) substr(s["1"], p, 1)
        printf "rs%d\tX\t%d\t%s\n", 8100000 + p, p, substr(s["X"], p, 1) substr(s["X"], p, 1)
      }
    }' "$fa"
}

convert het-miss S1_raw.txt 23andme < <(het_rows miss)
refused het-miss "250 heterozygous calls that miss the reference base under a build 37 header" \
  '250 of 250 heterozygous calls have neither allele equal to the GRCh37'

convert het-hit S1_raw.txt 23andme < <(het_rows hit)
if [ "$rc" -eq 0 ]; then pass "het-hit: chip-to-vcf.sh exited 0 on 250 heterozygous calls that carry the reference base"
else fail "het-hit: chip-to-vcf.sh exited ${rc} on 250 heterozygous calls that carry the reference base"; fi
if grep -q 'Build check: 250 of 250 heterozygous calls carry the GRCh37 base' "${WORK}/het-hit.log"; then
  pass "het-hit: the build check says the positions fit GRCh37"
else fail "het-hit: no build check line for 250 of 250"; fi
if [ -s "${WORK}/het-hit/S1/vcf/S1.vcf.gz" ]; then pass "het-hit: wrote the VCF"
else fail "het-hit: no VCF"; fi
if grep -q 'Input: 6000 genotyped positions' "${WORK}/het-hit.log"; then pass "het-hit: read all 6000 rows"
else fail "het-hit: did not read all 6000 rows"; fi

if [ "$FAILS" -gt 0 ]; then
  echo "test_chip_to_vcf: ${FAILS} check(s) failed"
  exit 1
fi
echo "test_chip_to_vcf: all checks passed"
