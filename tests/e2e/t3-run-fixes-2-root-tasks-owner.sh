#!/usr/bin/env bash
# The tasks that run as root under the docker profile (CPSR and CNVPYTOR)
# leave no file to root in their work directories, also when an attempt is
# killed. Nextflow's docker.fixOwnership gives the files back only after a
# task that succeeds, and with GNU chown's --from: in an image with BusyBox's
# chown, such as CNVpytor's, it fails, and the files stayed root's. Each of
# the two scripts now gives its work directory back itself on any exit.
#
# Both modules run their real scripts under the checkout's docker profile, as
# this runner's uid, in CNVpytor's image (BusyBox chown; CI has no room for the
# PCGR image, so CPSR runs in it too). Stand-ins stand for the tools: cpsr is
# killed on its first attempt, as for memory, then passes; cnvpytor (mounted
# over the image's own) writes the files the module expects. A third task,
# ROOT_CONTROL, runs as root in the same image with no give-back of its own:
# its files must stay root's, or the check below could not fail here. One
# staged input of CPSR and one of CNVPYTOR belong to root: the give-back
# changes the staged links only (chown -h), so those files stay root's.
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
mkdir -p "${NP}/bin" "${T}/marks" "${T}/res" "${T}/pcgr" "${T}/vep115"
chmod 777 "${T}/marks"
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
chmod +x "${NP}/bin/cpsr" "${T}/cnvpytor"
: > "${T}/res/gc_hg38.pytor"
: > "${T}/res/mask_hg38.pytor"
for f in once.vcf.gz s1.bam; do : > "${T}/${f}"; done
# One staged input of each root task belongs to root, as a shared reference
# can. The give-back must change the staged links, never what they point to,
# so these two stay root's.
docker run --rm --user 0:0 -v "${T}:/t" "$CNVPYTOR_IMAGE" sh -c ': > /t/once.vcf.gz.tbi && : > /t/s1.bam.bai' >/dev/null 2>&1
check_eq "two staged inputs are root's before the run" "$(stat -c %u "${T}/once.vcf.gz.tbi" "${T}/s1.bam.bai" 2>/dev/null | paste -sd' ' -)" "0 0"

cat > "${NP}/main.nf" <<NF
include { CPSR } from '${REPO}/modules/local/cpsr/main'
include { CNVPYTOR } from '${REPO}/modules/local/cnvpytor/main'

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
    ROOT_CONTROL()
}
NF
# The docker profile as nextflow.config sets it, plus two mounts for the
# stand-ins, and CNVpytor's image for every task.
cat > "${T}/run.config" <<NFCONF
docker.runOptions = '-u \${UID}:\${GROUPS} -e HOME=/tmp -v ${T}/cnvpytor:/usr/local/bin/cnvpytor:ro -v ${T}/marks:/marks'
process {
    withName: 'CPSR|ROOT_CONTROL' {
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

echo "+ nextflow run -profile docker (CPSR, CNVPYTOR, ROOT_CONTROL) as uid ${ME}"
rc=0
(cd "$T" && nextflow run "${NP}/main.nf" -ansi-log false -profile docker -c "${REPO}/nextflow.config" -c "${T}/run.config" \
  -with-trace "${T}/trace.txt" --outdir "${T}/out" --max_cpus 4 --max_memory 64.GB -work-dir "${T}/work") \
  > "${T}/nf.log" 2>&1 || rc=$?
check_eq "the Nextflow run exits 0" "$rc" 0
[ "$rc" -eq 0 ] || tail -n 60 "${T}/nf.log"
cat "${T}/trace.txt" 2>/dev/null

# workdir NAME_REGEX ATTEMPT: the work directory of that attempt of the task.
workdir() {
  awk -F'\t' -v pat="$1" -v a="$2" 'NR > 1 && $2 ~ pat && $3 == a { print $5; exit }' "${T}/trace.txt" 2>/dev/null
}
# foreign DIR: the files under DIR that are not this runner's, one per line.
foreign() { [ -d "$1" ] && find "$1" ! -user "$ME" 2>/dev/null | sed "s|^$1|.|"; }

CPSR1=$(workdir 'CPSR \(once\)' 1) CPSR2=$(workdir 'CPSR \(once\)' 2)
CNV=$(workdir 'CNVPYTOR \(s1\)' 1) CTL=$(workdir 'ROOT_CONTROL' 1)
check "CPSR's killed attempt and its retry are in the trace" test -d "$CPSR1" -a -d "$CPSR2"
check "the killed attempt wrote a file as root before it died" test -e "${CPSR1}/once.cpsr.partial.tsv"
check_eq "files left to another user by CPSR's killed attempt" "$(foreign "$CPSR1" | paste -sd' ' -)" ""
check_eq "files left to another user by CPSR's retry" "$(foreign "$CPSR2" | paste -sd' ' -)" ""
check "CNVPYTOR wrote its calls" test -s "${CNV}/s1_cnvs.txt"
check_eq "files left to another user by CNVPYTOR" "$(foreign "$CNV" | paste -sd' ' -)" ""
check_eq "the root-owned files the staged links point to stay root's" \
  "$(stat -c %u "${T}/once.vcf.gz.tbi" "${T}/s1.bam.bai" 2>/dev/null | paste -sd' ' -)" "0 0"
check "control: BusyBox's chown refused Nextflow's --from" grep -q "unrecognized option '--from'" "${CTL}/.command.err"
check_ge "control: a root task with no give-back leaves root's files" "$(foreign "$CTL" | grep -c .)" 1

# The control's files are root's: only a container can delete them.
[ -z "$CTL" ] || docker run --rm --user 0:0 -v "${CTL}:/c" "$CNVPYTOR_IMAGE" sh -c 'rm -rf /c/d /c/o /c/.command.trace' >/dev/null 2>&1
rm -rf "$T"
finish
