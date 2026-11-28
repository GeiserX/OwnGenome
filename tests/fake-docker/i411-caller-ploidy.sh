#!/usr/bin/env bash
# The alternative callers call a male's chrX and chrY with ploidy 1 outside
# the pseudoautosomal regions (assets/par_grch38.bed), as step 03 does, on a
# reference with chr1, chrX, chrY and chrM of GRCh38's lengths:
#   03a GATK       `male` cuts the chrX and chrY units at the PARs; the parts
#                  outside run as units of their own with --sample-ploidy 1,
#                  and no other unit gets it. So with INTERVALS (a region that
#                  crosses the end of PAR1 is cut there) and with
#                  SCATTER=false (two processes in turn, each with every
#                  CPU). `female` and no sex pass no ploidy and run one unit
#                  per contig. A male result is called again for a female run
#                  and the other way round; female and no sex share a result.
#                  With no sex, the sex run-all.sh recorded in
#                  logs/run_status.tsv counts, and run-all.sh with
#                  EXTRA_CALLERS=gatk and a male sample reaches 03a that way.
#                  A sex that is not male or female stops the step with exit 2
#                  before any container.
#   03b FreeBayes  the same cut, --ploidy 1 on the haploid units only.
#   03c Strelka2   `male` configures with --ploidy and a VCF of one CN=1
#                  record per stretch of chrX and chrY outside the PARs, its
#                  column named after the BAM's read group sample; a BAM with
#                  two samples in its read groups stops it. `female` passes no
#                  ploidy.
#   03d Octopus    `male` runs Octopus twice, the haploid regions with
#                  --contig-ploidies chrX=1 chrY=1 and the rest with MT=1
#                  chrM=1, and joins the two VCFs; `female` runs it once as
#                  before, with no ploidy and no region file.
#   All four      a cut at the PARs that fails part way (here at the chrX unit)
#                  stops a male run with exit 1 before any call, and leaves no
#                  result to reuse.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" THREADS=4
G=$GENOME_DIR
seed_reference "$G"
REF="${G}/reference/GRCh38_no_alt_analysis_set"
printf 'chr1\t248956422\t6\t60\t61\nchrX\t156040895\t253105752\t60\t61\nchrY\t57227415\t411747563\t60\t61\nchrM\t16569\t469928629\t60\t61\n' > "${REF}.fasta.fai"
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "${REF}.dict"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook
cat > "${CASE_WORK}/ploidy-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
# bgzf HOST_PATH: a finished bgzipped VCF (ends with the BGZF end-of-file block)
bgzf() {
  mkdir -p "$(dirname "$1")"
  { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' | gzip -c
    printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "$1"
}
opt() { sed -n "s/.* $1 \\([^ ]*\\) .*/\\1/p" <<<"$args"; }
# unit TOOL FLAG FILE: one line per call in calls.log: the tool, hap when the
# call carries FLAG (its haploid switch) or dip, and the regions of FILE.
unit() {
  local p=dip
  case "$args" in *" $2 "*) p=hap ;; esac
  printf '%s\t%s\t%s\n' "$1" "$p" "$(tr '\t\n' ':,' < "$(host_path "$3")")" >> "${CASE_WORK}/calls.log"
}
case "$args" in
  *" HaplotypeCaller "*)
    unit gatk "--sample-ploidy 1" "$(opt -L)"
    o=$(host_path "$(opt -O)"); bgzf "$o"; printf 'TBI\001' > "${o}.tbi"
    exit 0 ;;
  *" freebayes "*)
    unit freebayes "--ploidy 1" "$(opt --targets)"
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
    exit 0 ;;
  *" octopus "*)
    if [[ "$args" == *" --regions-file "* ]]; then unit octopus "--contig-ploidies chrX=1 chrY=1" "$(opt --regions-file)"; fi
    o=$(host_path "$(opt -o)"); bgzf "$o"; printf 'TBI\001' > "${o}.tbi"
    exit 0 ;;
  *" samtools view -H "*)
    printf '@HD\tVN:1.6\tSO:coordinate\n@RG\tID:lane1\tSM:NA24385\tPL:ILLUMINA\n'
    [ -z "${FAKE_TWO_SAMPLES:-}" ] || printf '@RG\tID:lane2\tSM:other\tPL:ILLUMINA\n'
    exit 0 ;;
  *" configureStrelkaGermlineWorkflow.py "*)
    d=$(host_path "$(opt --runDir)")
    mkdir -p "$d" && printf '#!/usr/bin/env python\n' > "${d}/runWorkflow.py"
    exit 0 ;;
  *"/runWorkflow.py "*)
    d=$(host_path /genome/sample1/vcf_strelka2/results/variants)
    mkdir -p "$d" && bgzf "${d}/variants.vcf.gz" && printf 'TBI\001' > "${d}/variants.vcf.gz.tbi"
    exit 0 ;;
  *" bcftools concat "*|*" bcftools sort "*)
    bgzf "$(host_path "$(opt -o)")"
    exit 0 ;;
  *" bcftools index "*)
    printf 'TBI\001' > "$(host_path "${*: -1}").tbi"
    exit 0 ;;
  *" bcftools stats "*)
    printf 'SN\t0\tnumber of records:\t42\n'
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/ploidy-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/ploidy-hook"
fresh() { : > "$FAKE_DOCKER_LOG"; : > "${CASE_WORK}/calls.log"; }
# calls TOOL: TOOL's calls of this run, `hap|dip<TAB>regions`, sorted.
calls() { awk -F'\t' -v t="$1" '$1 == t { print $2 "\t" $3 }' "${CASE_WORK}/calls.log" | LC_ALL=C sort; }
count() { grep -cE -- "$1" "$FAKE_DOCKER_LOG" || true; }
same() { [ "$2" = "$3" ] || fail "$1:
  got:
$2
  want:
$3"; }

# The cut of a whole-genome male run, the same for GATK and FreeBayes: PAR1
# and PAR2 stay diploid; 1-10,000, the stretch between the PARs and the end
# past PAR2 of each contig are haploid; chr1 and chrM are untouched.
MALE=$(printf '%s\n' \
  $'dip\tchr1:0:248956422,' \
  $'dip\tchrM:0:16569,' \
  $'dip\tchrX:10000:2781479,chrX:155701382:156030895,' \
  $'dip\tchrY:10000:2781479,chrY:56887902:57217415,' \
  $'hap\tchrX:0:10000,chrX:2781479:155701382,chrX:156030895:156040895,' \
  $'hap\tchrY:0:10000,chrY:2781479:56887902,chrY:57217415:57227415,' | LC_ALL=C sort)
NOSPLIT=$(printf '%s\n' $'dip\tchr1:0:248956422,' $'dip\tchrM:0:16569,' $'dip\tchrX:0:156040895,' $'dip\tchrY:0:57227415,' | LC_ALL=C sort)

# --- 03a GATK -----------------------------------------------------------------------------
GD="${G}/sample1/vcf_gatk"
fresh
run_expect 0 gatk-male "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 male
output_has gatk-male '^Ploidy: male, chrX and chrY haploid outside the PARs \(assets/par_grch38\.bed\)'
output_has gatk-male '6 unit\(s\), 2 at a time'
same "03a male: the units and their ploidy" "$(calls gatk)" "$MALE"
[ "$(count ' HaplotypeCaller .*--sample-ploidy 1 ')" -eq 2 ] || fail "03a male did not pass --sample-ploidy 1 to exactly the two haploid units"
docker_log_has 'bcftools concat -a -D .*scatter/002\.hap\.vcf\.gz .*scatter/003\.hap\.vcf\.gz' "03a did not join the haploid units"
grep -q '^INTERVALS= scatter=true sex=male bam=' "${GD}/sample1.run" || fail "03a male did not record sex=male: $(head -1 "${GD}/sample1.run")"
[ ! -e "${GD}/scatter" ] || fail "03a kept its scatter folder"
# The same call again keeps the VCF.
fresh
run_expect 0 gatk-male-again "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 male
output_has gatk-male-again 'Output already exists'
# A female run calls again, diploid, one unit per contig.
fresh
run_expect 0 gatk-female "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 female
output_has gatk-female '^Ploidy: female, every contig diploid'
output_lacks gatk-female 'Output already exists'
same "03a female: the units" "$(calls gatk)" "$NOSPLIT"
[ "$(count '--sample-ploidy')" -eq 0 ] || fail "03a female passed a ploidy"
grep -q '^INTERVALS= scatter=true bam=' "${GD}/sample1.run" || fail "03a female's record is not the no-sex record"
# No sex calls as female does, so the female result stands.
fresh
run_expect 0 gatk-none "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-none '^Ploidy: no sex given, chrX and chrY are called diploid'
output_has gatk-none 'Output already exists'
# A male run after the female one calls again.
fresh
run_expect 0 gatk-male-after-female "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 male
same "03a male after female" "$(calls gatk)" "$MALE"
# A region that crosses the end of PAR1 is cut there; a region in PAR1 stays
# diploid and one between the PARs is haploid.
fresh
# A region that starts one base before PAR1 keeps that base, haploid.
INTERVALS="chr1:1-100 chrX:10000-20000 chrX:1000001-1200000 chrX:73700001-74000000 chrY:2700001-3000000" \
  run_expect 0 gatk-intervals "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 male
same "03a male INTERVALS" "$(calls gatk)" "$(printf '%s\n' $'dip\tchr1:0:100,' $'dip\tchrX:10000:20000,' $'dip\tchrX:1000000:1200000,' \
  $'dip\tchrY:2700000:2781479,' $'hap\tchrX:9999:10000,' $'hap\tchrX:73700000:74000000,' $'hap\tchrY:2781479:3000000,' | LC_ALL=C sort)"
# SCATTER=false: two processes in turn, each with every CPU.
fresh
SCATTER=false run_expect 0 gatk-one "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 male
output_has gatk-one '2 unit\(s\), 1 at a time'
same "03a male SCATTER=false" "$(calls gatk)" "$(printf '%s\n' \
  $'dip\tchr1:0:248956422,chrX:10000:2781479,chrX:155701382:156030895,chrY:10000:2781479,chrY:56887902:57217415,chrM:0:16569,' \
  $'hap\tchrX:0:10000,chrX:2781479:155701382,chrX:156030895:156040895,chrY:0:10000,chrY:2781479:56887902,chrY:57217415:57227415,' | LC_ALL=C sort)"
[ "$(count ' --cpus 4 .* HaplotypeCaller ')" -eq 2 ] || fail "SCATTER=false did not give each of the two processes every CPU"
# No sex given, but run-all.sh recorded one for the sample: that one counts.
mkdir -p "${G}/sample1/logs"
printf '# run-all.sh\nmeta\tstarted_epoch\t1\nmeta\tdeclared_sex\tmale\n' > "${G}/sample1/logs/run_status.tsv"
rm -rf "$GD"
fresh
run_expect 0 gatk-status "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1
output_has gatk-status '^NOTE: sex male, as run-all.sh was given it \(.*/sample1/logs/run_status\.tsv\)'
same "03a with the sex from run_status.tsv" "$(calls gatk)" "$MALE"
# The argument wins over the recorded sex.
fresh
run_expect 0 gatk-arg-wins "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 female
same "03a female over a recorded male" "$(calls gatk)" "$NOSPLIT"
rm -f "${G}/sample1/logs/run_status.tsv"
# run-all.sh starts 03a with the sample name only, after the pipeline: a male
# run still reaches it as male, through the sex it recorded.
rm -rf "$GD"
fresh
run_rc run-all env SKIP_VALIDATION=true EXTRA_CALLERS=gatk "${SCRIPTS}/run-all.sh" sample1 male
grep -q '^Ploidy: male' "${G}/sample1/logs/03a-gatk-haplotypecaller.log" 2>/dev/null \
  || fail "run-all.sh EXTRA_CALLERS=gatk did not call 03a as male: $(cat "${G}/sample1/logs/03a-gatk-haplotypecaller.log" 2>&1)"
same "03a through run-all.sh, male" "$(calls gatk)" "$MALE"
rm -f "${G}/sample1/logs/run_status.tsv"
# Another word for the sex stops the step before any container.
fresh
run_expect 2 gatk-bad-sex "${SCRIPTS}/03a-gatk-haplotypecaller.sh" sample1 man
output_has gatk-bad-sex "ERROR: sex must be 'male' or 'female' \(or left out\), got 'man'"
[ "$(count '^run ')" -eq 0 ] || fail "03a started a container with a sex it does not know"

# --- 03b FreeBayes ------------------------------------------------------------------------
fresh
run_expect 0 freebayes-male "${SCRIPTS}/03b-freebayes.sh" sample1 male
output_has freebayes-male '^Ploidy: male'
same "03b male: the units and their ploidy" "$(calls freebayes)" "$MALE"
[ "$(count ' freebayes .*--ploidy 1 ')" -eq 2 ] || fail "03b male did not pass --ploidy 1 to exactly the two haploid units"
grep -q '^INTERVALS= scatter=true sex=male bam=' "${G}/sample1/vcf_freebayes/sample1.run" || fail "03b male did not record sex=male"
fresh
run_expect 0 freebayes-female "${SCRIPTS}/03b-freebayes.sh" sample1 female
same "03b female: the units" "$(calls freebayes)" "$NOSPLIT"
[ "$(count ' freebayes .*--ploidy')" -eq 0 ] || fail "03b female passed a ploidy"

# --- 03c Strelka2 -------------------------------------------------------------------------
SD="${G}/sample1/vcf_strelka2"
fresh
run_expect 0 strelka-male "${SCRIPTS}/03c-strelka2-germline.sh" sample1 male
docker_log_has 'configureStrelkaGermlineWorkflow\.py .*--ploidy /genome/sample1/vcf_strelka2/sample1\.ploidy\.vcf\.gz ' \
  "03c male did not configure Strelka2 with its ploidy VCF"
docker_log_has 'samtools view -H /genome/sample1/aligned/sample1_sorted\.bam' "03c did not read the BAM's read groups"
output_has strelka-male 'Ploidy VCF: .*sample1\.ploidy\.vcf\.gz \(6 haploid stretches, sample NA24385\)'
test -s "${SD}/sample1.ploidy.vcf.gz.tbi" || fail "03c did not index the ploidy VCF"
same "03c male: the ploidy VCF's records" "$(grep -v '^##' "${SD}/sample1.ploidy.vcf")" "$(printf '%s\n' \
  $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA24385' \
  $'chrX\t1\t.\tN\t<CNV>\t.\tPASS\tEND=10000\tCN\t1' \
  $'chrX\t2781479\t.\tN\t<CNV>\t.\tPASS\tEND=155701382\tCN\t1' \
  $'chrX\t156030895\t.\tN\t<CNV>\t.\tPASS\tEND=156040895\tCN\t1' \
  $'chrY\t1\t.\tN\t<CNV>\t.\tPASS\tEND=10000\tCN\t1' \
  $'chrY\t2781479\t.\tN\t<CNV>\t.\tPASS\tEND=56887902\tCN\t1' \
  $'chrY\t57217415\t.\tN\t<CNV>\t.\tPASS\tEND=57227415\tCN\t1')"
grep -q '^sex=male bam=' "${SD}/sample1.run" || fail "03c male did not record sex=male"
# A female run is configured again, with no ploidy and no read group lookup.
fresh
run_expect 0 strelka-female "${SCRIPTS}/03c-strelka2-germline.sh" sample1 female
docker_log_has 'configureStrelkaGermlineWorkflow\.py ' "03c female did not configure again after a male run"
[ "$(count '--ploidy')" -eq 0 ] || fail "03c female passed a ploidy"
[ "$(count 'samtools view -H')" -eq 0 ] || fail "03c female read the BAM's read groups"
# Two samples in the BAM's read groups: no name for the column, so it stops.
rm -rf "$SD"
fresh
FAKE_TWO_SAMPLES=1 run_expect 1 strelka-two-samples "${SCRIPTS}/03c-strelka2-germline.sh" sample1 male
output_has strelka-two-samples "ERROR: the BAM's read groups name 2 samples \(@RG SM\), not one"
[ "$(count 'configureStrelkaGermlineWorkflow')" -eq 0 ] || fail "03c configured Strelka2 without a ploidy VCF for a male"

# --- 03d Octopus --------------------------------------------------------------------------
OD="${G}/sample1/vcf_octopus"
octopus_runs() { count '^run image=[^ ]*octopus.* octopus -R '; }
fresh
run_expect 0 octopus-male "${SCRIPTS}/03d-octopus.sh" sample1 male
[ "$(octopus_runs)" -eq 2 ] || fail "03d male did not run Octopus twice"
[ "$(count ' octopus .*--regions-file /genome/sample1/vcf_octopus/ploidy/001\.hap\.regions --contig-ploidies chrX=1 chrY=1( |$)')" -eq 1 ] \
  || fail "03d male did not call the haploid regions with chrX=1 chrY=1"
[ "$(count ' octopus .*--regions-file /genome/sample1/vcf_octopus/ploidy/001\.regions --contig-ploidies MT=1 chrM=1( |$)')" -eq 1 ] \
  || fail "03d male did not call the rest with MT=1 chrM=1"
# Octopus reads chrom:begin-end zero-based, so the regions are the BED's.
same "03d male: the two runs' regions" \
  "$(calls octopus | awk -F'\t' '{ n = split($2, r, ","); for (i = 1; i <= n; i++) if (r[i] != "") print $1, r[i] }' | LC_ALL=C sort)" \
  "$(printf '%s\n' \
  'dip chr1:0-248956422' 'dip chrM:0-16569' 'dip chrX:10000-2781479' 'dip chrX:155701382-156030895' \
  'dip chrY:10000-2781479' 'dip chrY:56887902-57217415' \
  'hap chrX:0-10000' 'hap chrX:2781479-155701382' 'hap chrX:156030895-156040895' \
  'hap chrY:0-10000' 'hap chrY:2781479-56887902' 'hap chrY:57217415-57227415' | LC_ALL=C sort)"
docker_log_has 'bcftools concat -a -D .*/vcf_octopus/ploidy/001\.vcf\.gz /genome/sample1/vcf_octopus/ploidy/001\.hap\.vcf\.gz' "03d male did not join the two runs"
test -s "${OD}/sample1.vcf.gz.tbi" || fail "03d male ended without the index"
grep -q '^INTERVALS= sex=male bam=' "${OD}/sample1.run" || fail "03d male did not record sex=male"
[ ! -e "${OD}/ploidy" ] || fail "03d kept its ploidy folder"
# A female run calls again, once, as before: no region file, no ploidy.
fresh
run_expect 0 octopus-female "${SCRIPTS}/03d-octopus.sh" sample1 female
output_has octopus-female "^Ploidy: female, every contig diploid but chrY and chrM \(Octopus' default\)"
[ "$(octopus_runs)" -eq 1 ] || fail "03d female did not run Octopus once"
[ "$(count ' octopus .*(--regions-file|--contig-ploidies)')" -eq 0 ] || fail "03d female passed a region file or a ploidy"
# Two INTERVALS regions reach Octopus as two words of --regions, not one
# word with a space inside (the fake docker logs each argument %q-quoted).
fresh
INTERVALS="chr1:1-100 chrX:1000001-1200000" run_expect 0 octopus-female-intervals "${SCRIPTS}/03d-octopus.sh" sample1 female
[ "$(count ' octopus .*--regions chr1:1-100 chrX:1000001-1200000( |$)')" -eq 1 ] \
  || fail "03d female did not pass the two INTERVALS regions as two words: $(grep -E ' octopus ' "$FAKE_DOCKER_LOG" | tail -1)"

# --- A failed PAR cut ---------------------------------------------------------------------
# This awk fails when it cuts a unit that holds chrX at the PARs. Whole-genome
# units are cut in order, so GATK and FreeBayes fail after chr1 was cut.
REAL_AWK=$(command -v awk)
mkdir -p "${CASE_WORK}/broken-awk"
cat > "${CASE_WORK}/broken-awk/awk" <<SHIM
#!/usr/bin/env bash
for last in "\$@"; do :; done
case " \$* " in
  *"/par_grch38.bed "*) if grep -q "^chrX" "\$last" 2>/dev/null; then echo "awk: read error" >&2; exit 2; fi ;;
esac
exec "$REAL_AWK" "\$@"
SHIM
chmod +x "${CASE_WORK}/broken-awk/awk"
# broken_cut NAME SCRIPT DIR CALLER: a male run stops, calls nothing, records nothing.
broken_cut() {
  rm -rf "${G}/sample1/$3"
  fresh
  PATH="${CASE_WORK}/broken-awk:${PATH}" run_expect 1 "$1" "${SCRIPTS}/$2" sample1 male
  output_has "$1" "ERROR: could not cut .* at .*par_grch38\\.bed\\."
  [ "$(count "$4")" -eq 0 ] || fail "$2 called after a failed PAR cut"
  [ ! -e "${G}/sample1/$3/sample1.run" ] || fail "$2 recorded a run after a failed PAR cut"
}
broken_cut gatk-broken-cut 03a-gatk-haplotypecaller.sh vcf_gatk " HaplotypeCaller "
broken_cut freebayes-broken-cut 03b-freebayes.sh vcf_freebayes " freebayes "
broken_cut strelka-broken-cut 03c-strelka2-germline.sh vcf_strelka2 "configureStrelkaGermlineWorkflow"
broken_cut octopus-broken-cut 03d-octopus.sh vcf_octopus "^run image=[^ ]*octopus.* octopus -R "

echo "Each alternative caller calls a male's chrX and chrY haploid outside the PARs, and nothing else."
