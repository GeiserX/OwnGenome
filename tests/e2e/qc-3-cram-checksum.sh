#!/usr/bin/env bash
# Synthetic read content, including secondary/supplementary reads and tags.
# A quality, sequence or custom-tag change leaves flagstat unchanged. Each
# must stop --delete-bam and restore before a partial output is published.
. "$(dirname "$0")/lib.sh"
export GENOME_DIR="${CASE_TMP}/genome"
export REF_FASTA="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
C=checksum
A="${GENOME_DIR}/${C}/aligned"
rm -rf "$A"
mkdir -p "${GENOME_DIR}/reference" "$A" "${CASE_TMP}/fault"
printf '>chr1\n%s\n' "$(printf 'ACGT%.0s' {1..100})" > "$REF_FASTA"
# Include correct and incorrect MD/NM, missing MD/NM, a custom tag, a reverse
# read, a secondary read, a supplementary read, and placed/unplaced unmapped
# reads. Unmapped MAPQ is normalized by samtools' checksum sanitizer.
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:400\n@RG\tID:rg1\tSM:synthetic\n' > "${A}/in.sam"
for flag in 0 16 256 2048; do
  printf 'read%s\t%s\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:4\tNM:i:0\tXZ:Z:kept\tRG:Z:rg1\n' "$flag" "$flag"
done >> "${A}/in.sam"
printf 'wrong-md\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:0A3\tNM:i:1\nunmapped-placed\t4\tchr1\t9\t30\t*\t*\t0\t0\tACGT\tIIII\nunmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n' >> "${A}/in.sam"
sam faidx reference/GRCh38_no_alt_analysis_set.fasta
sam view -b -o "${C}/aligned/original.bam" "${C}/aligned/in.sam"
sam index "${C}/aligned/original.bam"

# After a real conversion, rewrite its output with one changed field. Both
# commands succeed and all read flags/counts stay the same. For restore the
# same wrapper changes the partial BAM instead.
cat > "${CASE_TMP}/fault/docker" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail
self_dir="$(cd "$(dirname "$0")" && pwd)"
IFS=: read -r -a dirs <<< "$PATH"
keep=(); real=""
for d in "${dirs[@]}"; do
  [ "$d" = "$self_dir" ] && continue
  keep+=("$d")
  if [ -z "$real" ] && [ -x "${d}/docker" ]; then real="${d}/docker"; fi
done
PATH=$(IFS=:; echo "${keep[*]}"); export PATH
[ -n "$real" ] || exit 127
"$real" "$@"
args=" $* "
case "$args" in
  *" samtools view "*" -o "*".part.${FAULT_FORMAT} "*) ;;
  *) exit 0 ;;
esac
out=$(grep -oE "/genome/[^ ]+\.part\.${FAULT_FORMAT}" <<< "$args" | head -n 1)
"$real" run --rm -i -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" "$SAMTOOLS_IMAGE" \
  bash -euo pipefail -c '
    out=$1; fault=$2; ref=/genome/reference/GRCh38_no_alt_analysis_set.fasta
    samtools view -T "$ref" --input-fmt-option decode_md=0 "$out" |
      awk -F "\t" -v OFS="\t" -v fault="$fault" '\''
        fault == "reorder" && NR <= 4 {
          # All four share one coordinate, so reversing them stays sorted.
          t=$NF; $NF=$(NF-1); $(NF-1)=t
          rows[NR]=$0
          if (NR == 4) for (j=4; j>=1; j--) print rows[j]
          next
        }
        NR == 1 && fault == "quality" { $11 = "H" substr($11, 2) }
        NR == 3 && fault == "sequence" { $10 = "T" substr($10, 2) }
        NR == 4 && fault == "tag" { for (i=12; i<=NF; i++) if ($i ~ /^XZ:/) $i="XZ:Z:changed" }
        { print }
      '\'' > /tmp/changed.sam
    samtools view -H "$out" > /tmp/changed-header.sam
    cat /tmp/changed-header.sam /tmp/changed.sam > /tmp/changed-all.sam
    if [[ $out == *.cram ]]; then
      samtools view -C -T "$ref" --output-fmt-option store_md=1 --output-fmt-option store_nm=1 -o "$out" /tmp/changed-all.sam
    else
      samtools view -b -o "$out" /tmp/changed-all.sam
    fi
    echo "fault: $fault in $out" >&2
  ' bash "$out" "$FAULT"
WRAP
chmod +x "${CASE_TMP}/fault/docker"
export SAMTOOLS_IMAGE

seed_bam() {
  cp -f "${A}/original.bam" "${A}/${C}_sorted.bam"
  cp -f "${A}/original.bam.bai" "${A}/${C}_sorted.bam.bai"
}
for fault in quality sequence tag; do
  seed_bam
  FAULT=$fault FAULT_FORMAT=cram PATH="${CASE_TMP}/fault:${PATH}" run_step 34-cram-archive.sh "$C" --delete-bam
  check "${fault}: fault injected" has "^fault: ${fault} " "$(cat "$STEP_LOG")"
  check "${fault}: archive refuses changed content" test "$STEP_RC" -ne 0
  check "${fault}: BAM is unchanged" cmp -s "${A}/original.bam" "${A}/${C}_sorted.bam"
  check "${fault}: BAM index stays" test -s "${A}/${C}_sorted.bam.bai"
  check "${fault}: no partial CRAM" test ! -e "${A}/${C}_sorted.part.cram"
  check "${fault}: no published CRAM" test ! -e "${A}/${C}_sorted.cram"
done
seed_bam
FAULT=reorder FAULT_FORMAT=cram PATH="${CASE_TMP}/fault:${PATH}" run_step 34-cram-archive.sh "$C" --delete-bam
check "reordered read/tag control was injected" has "^fault: reorder " "$(cat "$STEP_LOG")"
check_step_exit 34-cram-archive.sh
check "faithful archive deletes BAM" test ! -e "${A}/${C}_sorted.bam"
check "archive keeps checksum evidence" test -s "${A}/${C}_sorted.cram.checksum"
for fault in quality sequence tag; do
  rm -f "${A}/${C}_sorted.bam" "${A}/${C}_sorted.bam.bai"
  FAULT=$fault FAULT_FORMAT=bam PATH="${CASE_TMP}/fault:${PATH}" run_step 34-cram-archive.sh "$C" --restore
  check "${fault}: restore fault injected" has "^fault: ${fault} " "$(cat "$STEP_LOG")"
  check "${fault}: restore refuses changed content" test "$STEP_RC" -ne 0
  check "${fault}: restore leaves no BAM" test ! -e "${A}/${C}_sorted.bam"
  check "${fault}: restore leaves no partial BAM" test ! -e "${A}/${C}_sorted.part.bam"
done
rm -f "${A}/${C}_sorted.bam" "${A}/${C}_sorted.bam.bai"
run_step 34-cram-archive.sh "$C" --restore
check_step_exit 34-cram-archive.sh
check "faithful restore writes BAM and index" test -s "${A}/${C}_sorted.bam" -a -s "${A}/${C}_sorted.bam.bai"
# No generated MD/NM tags are added and all explicit tags survive, independent
# of tag order. checksum sanitizes only format-level normalization.
sam checksum -PCMc -b 0xfff -f0 -F0 -z all,cigarx -t '*' "${C}/aligned/original.bam" | tail -n +2 | LC_ALL=C sort > "${CASE_TMP}/original.checksum"
sam checksum -PCMc -b 0xfff -f0 -F0 -z all,cigarx -t '*' "${C}/aligned/${C}_sorted.bam" | tail -n +2 | LC_ALL=C sort > "${CASE_TMP}/restored.checksum"
check "round trip preserves read content and explicit tags" cmp -s "${CASE_TMP}/original.checksum" "${CASE_TMP}/restored.checksum"
# A CRAM written without store_md/store_nm (an older step 34, or one made
# elsewhere) drops MD and NM wherever they match the reference. Restoring it
# regenerates them, so read0's NM:i:0 is back in the BAM.
plain_cram() {
  sam view -C -T reference/GRCh38_no_alt_analysis_set.fasta -o "$1" "${C}/aligned/original.bam"
  sam index "$1"
}
nm_of_read0() { sam view "$1" | awk -F '\t' '$1 == "read0" { for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) print $i }'; }
plain_cram "${C}/aligned/${C}_sorted.cram"
check "the plain CRAM stores no NM for read0" test -z "$(sam view -T reference/GRCh38_no_alt_analysis_set.fasta --input-fmt-option decode_md=0 "${C}/aligned/${C}_sorted.cram" | awk -F '\t' '$1 == "read0" && /NM:i:/')"
rm -f "${A}/${C}_sorted.bam" "${A}/${C}_sorted.bam.bai"
run_step 34-cram-archive.sh "$C" --restore
check_step_exit 34-cram-archive.sh
check_eq "restore of a plain CRAM regenerates NM" "$(nm_of_read0 "${C}/aligned/${C}_sorted.bam")" "NM:i:0"
# Run the module's actual script blocks on the same synthetic records. Only
# Nextflow interpolation and version metadata are omitted by this harness.
# Literal Nextflow placeholders are substituted by sed, not by this shell.
# shellcheck disable=SC2016
module_script() {
  awk -v process="$1" '
    $0 == "process " process " {" {found=1}
    found && /^    script:/ {script=1; next}
    script && /^    """/ {body=1; next}
    body && /cat <<-END_VERSIONS/ {exit}
    body {print}
  ' "${REPO}/modules/local/cram_archive/main.nf" |
    sed -e 's/${meta.id}/module/g' -e 's/${task.cpus}/1/g' \
        -e 's|${reference}|reference/GRCh38_no_alt_analysis_set.fasta|g' \
        -e 's|${bam}|checksum/aligned/original.bam|g' \
        -e 's|${cram}|module_sorted.cram|g' \
        -e 's/\\\\/\\/g' -e 's/\\\$/\$/g'
}
# A successful command with no report must never certify a conversion. The
# other samtools verbs still run for real, so flagstat and quickcheck pass.
# The function expands its arguments in the container shell.
# shellcheck disable=SC2016
EMPTY_CHECKSUM='samtools() { if [ "$1" = checksum ]; then return 0; else command samtools "$@"; fi; }'
# A checksum that differs on the CRAM side only: the last report row's
# combined value is replaced. The module must fail on the diff alone.
# shellcheck disable=SC2016
CHANGED_CHECKSUM='samtools() { if [ "$1" = checksum ]; then case "${@: -1}" in *.cram) command samtools "$@" | sed "\$s/[0-9a-f]\{8\}\$/00000000/" ;; *) command samtools "$@" ;; esac; else command samtools "$@"; fi; }'
for process in CRAM_ARCHIVE CRAM_TO_BAM; do
  body=$(module_script "$process")
  rc=0
  in_genome "$SAMTOOLS_IMAGE" bash -euo pipefail -c "$body" || rc=$?
  check_eq "${process}: faithful module script passes" "$rc" 0
  rc=0
  in_genome "$SAMTOOLS_IMAGE" bash -euo pipefail -c "${EMPTY_CHECKSUM}
${body}" || rc=$?
  check "${process}: empty checksum reports fail" test "$rc" -ne 0
  rc=0
  in_genome "$SAMTOOLS_IMAGE" bash -euo pipefail -c "${CHANGED_CHECKSUM}
${body}" > "${CASE_TMP}/changed-${process}.log" 2>&1 || rc=$?
  check "${process}: a changed CRAM checksum fails" test "$rc" -ne 0
  check "${process}: the failure is the content check" has "content .*differs" "$(cat "${CASE_TMP}/changed-${process}.log")"
done
# CRAM_TO_BAM on a CRAM without stored MD/NM regenerates them, like --restore.
plain_cram module_sorted.cram
rc=0
in_genome "$SAMTOOLS_IMAGE" bash -euo pipefail -c "$(module_script CRAM_TO_BAM)" || rc=$?
check_eq "CRAM_TO_BAM: plain CRAM passes" "$rc" 0
check_eq "CRAM_TO_BAM: plain CRAM regenerates NM" "$(nm_of_read0 module_from_cram.bam)" "NM:i:0"
finish
