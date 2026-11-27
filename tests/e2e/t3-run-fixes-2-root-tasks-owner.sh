#!/usr/bin/env bash
# The tasks that run as root under the docker profile (CPSR, CNVPYTOR and
# MITO_VARIANTS) leave no file to root in their work directories, also when
# an attempt is killed, and whatever docker.fixOwnership says. Nextflow's
# docker.fixOwnership gives the files back only after a task that succeeds,
# and with GNU chown's --from: in an image with BusyBox's chown, such as
# CNVpytor's, it fails, and the files stayed root's. Each of the three
# scripts gives its work directory back itself on any exit, to the owner of
# that directory, so a run with docker.fixOwnership = false (no NXF_OWNER)
# leaves none to root either.
#
# The modules run their real scripts under the checkout's docker profile, as
# this runner's uid, in CNVpytor's image (BusyBox chown; CI has no room for the
# PCGR or GATK images, so CPSR and MITO_VARIANTS run in it too), twice: with
# docker.fixOwnership as nextflow.config sets it (true), then false.
# Stand-ins stand for the tools: cpsr is killed on its first attempt, as for
# memory, then passes; cnvpytor (mounted over the image's own) writes the
# files the module expects; gatk writes each file its -O names (MITO_VARIANTS'
# rotated pass then stops at og, which the image lacks, and the linear calls
# go on). A fourth task, ROOT_CONTROL, runs as root in the same image with no
# give-back of its own: its files must stay root's, or the check below could
# not fail here. One staged input of CPSR and one of CNVPYTOR belong to root:
# the give-back changes the staged links only (chown -h), so those files stay
# root's.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
if [ "$(id -u)" -eq 0 ]; then
  NO_CHECKS_REASON="this runner runs as root, so a root-owned file cannot be told from its own"
  finish
fi

ME=$(id -u)
T="${CASE_TMP}/nf"
NP="${T}/project"
rm -rf "$T"
mkdir -p "${NP}/bin" "${T}/res" "${T}/pcgr" "${T}/vep115" "${T}/ref"
cp "${REPO}/bin/cpsr_sample_id" "${NP}/bin/"

cat > "${NP}/bin/cpsr" <<'STANDIN'
#!/bin/bash
id="" dir=""
while [ $# -gt 0 ]; do case "$1" in --sample_id) id=$2; shift ;; --output_dir) dir=$2; shift ;; esac; shift; done
echo partial > "$dir/$id.cpsr.partial.tsv"
[ -e /marks/cpsr-once ] || { touch /marks/cpsr-once; kill -9 $$; }
touch "$dir/$id.cpsr.grch38.html" "$dir/$id.cpsr.grch38.classification.tsv.gz"
STANDIN
cat > "${T}/cnvpytor" <<'STANDIN'
#!/bin/bash
# The module's five cnvpytor calls: -rd, -his, -partition, -call, -view.
root=$2 op=$3
case "$op" in
  -rd) echo pytor > "$root" ;;
  -call) printf 'deletion\tchr1:1-1000\t1000\t0.5\n' ;;
  -view) f=$(sed -n 's/^set print_filename //p'); echo '##fileformat=VCFv4.2' > "$f" ;;
esac
STANDIN
cat > "${T}/gatk" <<'STANDIN'
#!/bin/bash
# Each gatk call of MITO_VARIANTS writes the file its -O names (Mutect2 its
# .stats beside it too) and the uid it ran as.
tool="" out=""
while [ $# -gt 0 ]; do
  case "$1" in
    --java-options) shift ;;
    -O) out=$2; shift ;;
    -*) ;;
    *) [ -n "$tool" ] || tool=$1 ;;
  esac
  shift
done
[ -z "$out" ] || echo x > "$out"
[ "$tool" != Mutect2 ] || echo x > "${out}.stats"
id -u > gatk.uid
STANDIN
chmod +x "${NP}/bin/cpsr" "${T}/cnvpytor" "${T}/gatk"
: > "${T}/res/gc_hg38.pytor"
: > "${T}/res/mask_hg38.pytor"
for f in once.vcf.gz s1.bam m1.bam m1.bam.bai m1.summary.txt m1.dist.txt ref/r.fasta ref/r.fasta.fai ref/r.dict; do : > "${T}/${f}"; done
# One staged input of each of CPSR and CNVPYTOR belongs to root, as a shared
# reference can. The give-back must change the staged links, never what they
# point to, so these two stay root's.
docker run --rm --user 0:0 -v "${T}:/t" "$CNVPYTOR_IMAGE" sh -c ': > /t/once.vcf.gz.tbi && : > /t/s1.bam.bai' >/dev/null 2>&1
check_eq "two staged inputs are root's before the run" "$(stat -c %u "${T}/once.vcf.gz.tbi" "${T}/s1.bam.bai" 2>/dev/null | paste -sd' ' -)" "0 0"

cat > "${NP}/main.nf" <<NF
include { CPSR } from '${REPO}/modules/local/cpsr/main'
include { CNVPYTOR } from '${REPO}/modules/local/cnvpytor/main'
include { MITO_VARIANTS } from '${REPO}/modules/local/mito_variants/main'

process ROOT_CONTROL {
    containerOptions '--network none --user 0:0'
    output:
    path 'o'
    script:
    """
    mkdir d
    echo x > o
    echo y > d/f
    """
}

workflow {
    CPSR(channel.of([[id: 'once'], file('${T}/once.vcf.gz'), file('${T}/once.vcf.gz.tbi')]), file('${T}/pcgr'), file('${T}/vep115'))
    CNVPYTOR(channel.of([[id: 's1'], file('${T}/s1.bam'), file('${T}/s1.bam.bai')]), file('${T}/res'))
    MITO_VARIANTS(channel.of([[id: 'm1'], file('${T}/m1.bam'), file('${T}/m1.bam.bai'),
                                file('${T}/m1.summary.txt'), file('${T}/m1.dist.txt')]),
                  file('${T}/ref/r.fasta'), file('${T}/ref/r.fasta.fai'), file('${T}/ref/r.dict'))
    ROOT_CONTROL()
}
NF

# workdir LEG TASK ATTEMPT: the work directory of that attempt of the task
# named TASK in that leg's run (a plain string, not a regex: in a -v value
# gawk turns "\(" into "(").
workdir() {
  awk -F'\t' -v task="$2" -v a="$3" 'NR > 1 && $2 == task && $3 == a { print $5; exit }' "${T}/$1/trace.txt" 2>/dev/null
}
# foreign DIR: the files under DIR that are not this runner's, one per line; a
# line naming DIR when it is no directory, so a lost work directory fails the check.
foreign() { [ -d "$1" ] || { echo "no directory '$1'"; return; }; find "$1" ! -user "$ME" 2>/dev/null | sed "s|^$1|.|"; }

# leg NAME FIX_OWNERSHIP: one run of the four tasks, its checks named NAME.
leg() {
  local n=$1 L="${T}/$1" rc=0 CPSR1 CPSR2 CNV MITO CTL
  mkdir -p "${L}/marks"
  chmod 777 "${L}/marks"
  # The docker profile as nextflow.config sets it, plus the mounts of the
  # stand-ins, and CNVpytor's image for every task.
  cat > "${L}/run.config" <<NFCONF
docker.runOptions = '-u \${UID}:\${GROUPS} -e HOME=/tmp -v ${T}/cnvpytor:/usr/local/bin/cnvpytor:ro -v ${T}/gatk:/usr/local/bin/gatk:ro -v ${L}/marks:/marks'
docker.fixOwnership = $2
process {
    withName: 'CPSR|MITO_VARIANTS|ROOT_CONTROL' {
        container = '${CNVPYTOR_IMAGE}'
    }
}
executor.cpus = 32
executor.memory = '256.GB'
trace.raw = true
trace.overwrite = true
trace.fields = 'task_id,name,attempt,exit,workdir'
timeline.enabled = false
report.enabled = false
dag.enabled = false
NFCONF
  echo "+ ${n}: nextflow run -profile docker (CPSR, CNVPYTOR, MITO_VARIANTS, ROOT_CONTROL) as uid ${ME}, docker.fixOwnership = $2"
  (cd "$L" && nextflow run "${NP}/main.nf" -ansi-log false -profile docker -c "${REPO}/nextflow.config" -c "${L}/run.config" \
    -with-trace "${L}/trace.txt" --outdir "${L}/out" --max_cpus 4 --max_memory 64.GB -work-dir "${L}/work") \
    > "${L}/nf.log" 2>&1 || rc=$?
  check_eq "${n}: the Nextflow run exits 0" "$rc" 0
  [ "$rc" -eq 0 ] || tail -n 60 "${L}/nf.log"
  cat "${L}/trace.txt" 2>/dev/null

  CPSR1=$(workdir "$n" 'CPSR (once)' 1) CPSR2=$(workdir "$n" 'CPSR (once)' 2)
  CNV=$(workdir "$n" 'CNVPYTOR (s1)' 1) MITO=$(workdir "$n" 'MITO_VARIANTS (m1)' 1) CTL=$(workdir "$n" 'ROOT_CONTROL' 1)
  check "${n}: CPSR's killed attempt and its retry are in the trace" test -d "$CPSR1" -a -d "$CPSR2"
  check "${n}: the killed attempt wrote a file as root before it died" test -e "${CPSR1}/once.cpsr.partial.tsv"
  check_eq "${n}: files left to another user by CPSR's killed attempt" "$(foreign "$CPSR1" | paste -sd' ' -)" ""
  check_eq "${n}: files left to another user by CPSR's retry" "$(foreign "$CPSR2" | paste -sd' ' -)" ""
  check "${n}: CNVPYTOR wrote its calls" test -s "${CNV}/s1_cnvs.txt"
  check_eq "${n}: files left to another user by CNVPYTOR" "$(foreign "$CNV" | paste -sd' ' -)" ""
  check "${n}: MITO_VARIANTS wrote its calls" test -s "${MITO}/m1_chrM_filtered.vcf.gz"
  check_eq "${n}: MITO_VARIANTS ran as root" "$(cat "${MITO}/gatk.uid" 2>/dev/null)" 0
  check_eq "${n}: files left to another user by MITO_VARIANTS" "$(foreign "$MITO" | paste -sd' ' -)" ""
  if [ "$2" = true ]; then
    check "${n}: control: BusyBox's chown refused Nextflow's --from" grep -q "unrecognized option '--from'" "${CTL}/.command.err"
  fi
  check_ge "${n}: control: a root task with no give-back leaves root's files" "$(find "$CTL" ! -user "$ME" 2>/dev/null | grep -c .)" 1
  # The control's files are root's: only a container can delete them.
  [ -z "$CTL" ] || docker run --rm --user 0:0 -v "${CTL}:/c" "$CNVPYTOR_IMAGE" sh -c 'rm -rf /c/d /c/o /c/.command.trace' >/dev/null 2>&1
}

leg fixownership true
leg no-fixownership false
check_eq "the root-owned files the staged links point to stay root's" \
  "$(stat -c %u "${T}/once.vcf.gz.tbi" "${T}/s1.bam.bai" 2>/dev/null | paste -sd' ' -)" "0 0"

rm -rf "$T"
finish
