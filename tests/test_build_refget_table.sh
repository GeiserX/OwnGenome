#!/usr/bin/env bash
# scripts/ci/build-refget-table.sh never cuts the checked-in table down and
# never publishes a row without its MD5. A fake curl on PATH serves a
# made-up chr21 (ACGTACGTNN) and an assembly report, so nothing is
# downloaded:
#   - CHROMS without an OUT, with OUT the checked-in table, or with OUT the
#     table through a symlinked folder, exits 2 before any download, and the
#     table is unchanged; so does an OUT that is a folder;
#   - an md5 tool that fails stops the run, and OUT is not written;
#   - with a working md5 tool, OUT holds chr21's row and no staging file is left.
# The script runs from a copy of its checkout layout, so a regression can only
# cut the copy's table.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "${WORK}/repo/scripts/ci" "${WORK}/repo/assets/ownhealth"
cp "${REPO}/scripts/ci/build-refget-table.sh" "${WORK}/repo/scripts/ci/"
cp "${REPO}/assets/ownhealth/grch38_refget.tsv" "${WORK}/repo/assets/ownhealth/"
SCRIPT=${WORK}/repo/scripts/ci/build-refget-table.sh
TABLE=${WORK}/repo/assets/ownhealth/grch38_refget.tsv
ln -s "${WORK}/repo/assets/ownhealth" "${WORK}/alias"
fail() { echo "FAIL: $*" >&2; exit 1; }

mkdir -p "${WORK}/bin" "${WORK}/badmd5"
cat > "${WORK}/bin/curl" <<'FAKE'
#!/usr/bin/env bash
# The fake download: the URL is the last-but-two argument, the file after -o.
echo "$*" >> "${CURL_LOG:?}"
url="" out=""
while [ $# -gt 0 ]; do case $1 in -o) out=$2; shift 2 ;; http*) url=$1; shift ;; *) shift ;; esac; done
case $url in
  *assembly_report.txt) printf '21\tassembled-molecule\t21\tChromosome\tCM000683.2\t=\tNC_000021.9\tPrimary Assembly\t10\tchr21\n' > "$out" ;;
  *chr21.fna.gz) printf '>NC_000021.9 Homo sapiens chromosome 21\nACGTA\ncgtNN\n' | gzip -c > "$out" ;;
  *) exit 22 ;;
esac
FAKE
for t in md5sum md5; do printf '#!/usr/bin/env bash\nexit 1\n' > "${WORK}/badmd5/$t"; done
chmod +x "${WORK}/bin/curl" "${WORK}/badmd5/md5sum" "${WORK}/badmd5/md5"
export CURL_LOG=${WORK}/curl.log
: > "$CURL_LOG"
before=$(cksum < "$TABLE")

rc=0; CHROMS=21 PATH="${WORK}/bin:$PATH" bash "$SCRIPT" 2>"${WORK}/err" || rc=$?
[ "$rc" -eq 2 ] || fail "CHROMS without OUT exited ${rc}, not 2"
grep -q 'give an OUT other than' "${WORK}/err" || fail "no message: $(cat "${WORK}/err")"
rc=0; CHROMS=21 PATH="${WORK}/bin:$PATH" bash "$SCRIPT" "$TABLE" 2>/dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "CHROMS with OUT the checked-in table exited ${rc}, not 2"
rc=0; CHROMS=21 PATH="${WORK}/bin:$PATH" bash "$SCRIPT" "${WORK}/alias/grch38_refget.tsv" 2>/dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "CHROMS with OUT the table through a symlinked folder exited ${rc}, not 2"
rc=0; CHROMS=21 PATH="${WORK}/bin:$PATH" bash "$SCRIPT" "${WORK}/bin" 2>/dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "an OUT that is a folder exited ${rc}, not 2"
[ ! -s "$CURL_LOG" ] || fail "it downloaded before refusing: $(cat "$CURL_LOG")"
[ "$(cksum < "$TABLE")" = "$before" ] || fail "the checked-in table changed"

rc=0; CHROMS=21 PATH="${WORK}/badmd5:${WORK}/bin:$PATH" bash "$SCRIPT" "${WORK}/bad.tsv" 2>/dev/null || rc=$?
[ "$rc" -ne 0 ] || fail "a failing md5 tool still ended with exit 0"
[ ! -e "${WORK}/bad.tsv" ] || fail "a failing md5 tool still wrote the table: $(tail -1 "${WORK}/bad.tsv")"

CHROMS=21 PATH="${WORK}/bin:$PATH" bash "$SCRIPT" "${WORK}/ok.tsv" >/dev/null
row=$(awk -F'\t' '$1 == "chr21"' "${WORK}/ok.tsv"); [ -n "$row" ] || fail "no chr21 row"
md5=$(printf 'ACGTACGTNN' | { md5sum 2>/dev/null || md5 -q; } | cut -d' ' -f1)
[ "$(cut -f5 <<<"$row")" = "$md5" ] || fail "chr21's MD5 is not the sequence's: ${row}"
[ "$(cut -f3 <<<"$row")" = 10 ] || fail "chr21's length: ${row}"
left=$(find "$WORK" -maxdepth 1 -name 'ok.tsv.*')
[ -z "$left" ] || fail "staging files left: ${left}"
echo "build-refget-table.sh keeps the checked-in table whole and writes no row without its MD5."
