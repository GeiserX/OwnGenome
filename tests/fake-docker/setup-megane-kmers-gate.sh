#!/usr/bin/env bash
# setup.sh --megane-kmers (step 50's k-mer set, built by og in OG_IMAGE):
#   1. The RAM gate on podman. Podman has no top-level .MemTotal in
#      `podman info` (it answers with a template error); its memory is
#      .Host.MemTotal. An engine of 2 GB must stop the build with the RAM
#      warning before og runs, not skip the gate.
#   2. A build killed on the way leaves og's hidden temporary files
#      (.<name>.mk.tmp-<pid>-<n> and .<name>.mi.tmp-<pid>-<n>). The next
#      build removes them before it starts, and leaves other files alone.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
DIR="${GENOME_DIR}/reference/megane"
MK="${DIR}/GRCh38_no_alt_analysis_set.mk"
MI="${DIR}/GRCh38_no_alt_analysis_set.mi"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

# The hook below plays og, so a stand-in for each Linux build bin/.og lacks
# is enough for setup.sh's og check; it goes when the case ends.
OGDIR="${REPO_ROOT}/bin/.og"
TREE=$("${REPO_ROOT}/rust/og.sh" tree-hash) || fail "could not compute the rust/ tree hash"
STANDINS=()
trap 'for f in "${STANDINS[@]}"; do rm -f "$f"; done' EXIT
mkdir -p "$OGDIR"
for a in x86_64 aarch64; do
  [ ! -e "${OGDIR}/og-linux-${a}" ] || continue
  printf 'stand-in og-tree-hash:%s;\n' "$TREE" > "${OGDIR}/og-linux-${a}"; chmod 755 "${OGDIR}/og-linux-${a}"
  STANDINS+=("${OGDIR}/og-linux-${a}")
done
"${REPO_ROOT}/rust/og.sh" check --required > /dev/null || fail "bin/.og does not pass setup.sh's og check"

# --- 1. podman of 2 GB: the gate reads .Host.MemTotal and stops ----------------------
mkdir -p "${CASE_WORK}/engine"
cat > "${CASE_WORK}/engine/podman" <<'ENGINE'
#!/usr/bin/env bash
# podman's answers to `info --format`; every other call goes to the fake docker.
if [ "${1:-}" = info ]; then
  case " $* " in
    *'{{.MemTotal}}'*)
      echo 'Error: template: info:1:2: executing "info" at <.MemTotal>: can'\''t evaluate field MemTotal in type *define.Info' >&2
      exit 125 ;;
    *'{{.Host.MemTotal}}'*) echo 2147483648; exit 0 ;;
  esac
fi
exec docker "$@"
ENGINE
chmod +x "${CASE_WORK}/engine/podman"
: > "$FAKE_DOCKER_LOG"
PATH="${CASE_WORK}/engine:${PATH}" CONTAINER_ENGINE=podman \
  run_expect 1 podman-small "${SCRIPTS}/setup.sh" --megane-kmers "$GENOME_DIR"
output_has podman-small "MEGAnE's k-mer set needs about 6 GB of RAM and podman has 2 GB\."
output_lacks podman-small "Building MEGAnE's k-mer set"
if awk -v img="run image=${OG_IMAGE}" 'index($0, img) { f = 1 } END { exit !f }' "$FAKE_DOCKER_LOG"; then
  fail "og ran on a podman of 2 GB: the RAM gate did not read .Host.MemTotal"
fi
[ ! -e "$MK" ] || fail "the gate stopped the build, yet ${MK} exists"

# --- 2. stale temporary files of a killed build are removed --------------------------
mkdir -p "$DIR"
printf 'partial' > "${DIR}/.GRCh38_no_alt_analysis_set.mk.tmp-1-0"
printf 'partial' > "${DIR}/.GRCh38_no_alt_analysis_set.mk.tmp-1-2"
printf 'partial' > "${DIR}/.GRCh38_no_alt_analysis_set.mi.tmp-1-1"
printf 'keep' > "${DIR}/.other.mk.tmp-1-0"
printf 'keep' > "${DIR}/notes.txt"
use_output_hook
# The hook plays og: it checks the temporaries are gone when og starts and
# writes the two files og writes.
cat > "${CASE_WORK}/hook-og" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
[ "${1:-}" = /pgp-bin/og ] || exit 0
out="" prefix=""
while [ $# -gt 0 ]; do
  case "$1" in
    --out) out=$2; shift 2 ;;
    --prefix) prefix=$2; shift 2 ;;
    *) shift ;;
  esac
done
h=$(host_path "$out")
for f in "${h}/.${prefix}.mk.tmp-"* "${h}/.${prefix}.mi.tmp-"*; do
  if [ -e "$f" ]; then
    echo "fake og: the stale ${f##*/} is still there when og starts" >&2
    exit 9
  fi
done
printf '\0\0\0\0\0\0\0\1' > "${h}/${prefix}.mk"
printf '1\n' > "${h}/${prefix}.mi"
echo "og setup megane-kmers: 1 repeated 32-mers of 2 in ${prefix}.mk (1 passes on 1 threads, 0 s)"
HOOK
chmod +x "${CASE_WORK}/hook-og"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-og"
: > "$FAKE_DOCKER_LOG"
run_expect 0 rebuild "${SCRIPTS}/setup.sh" --megane-kmers "$GENOME_DIR"
output_has rebuild "\[OK\] MEGAnE k-mer set \(step 50, 1 repeated 32-mers\)"
for f in .GRCh38_no_alt_analysis_set.mk.tmp-1-0 .GRCh38_no_alt_analysis_set.mk.tmp-1-2 .GRCh38_no_alt_analysis_set.mi.tmp-1-1; do
  [ ! -e "${DIR}/${f}" ] || fail "setup.sh left the stale ${f} of a killed build"
done
[ -e "${DIR}/.other.mk.tmp-1-0" ] || fail "setup.sh removed .other.mk.tmp-1-0, the temporary file of another set"
[ -e "${DIR}/notes.txt" ] || fail "setup.sh removed notes.txt"
[ -s "$MK" ] && [ -s "$MI" ] || fail "setup.sh wrote no ${MK} and ${MI}"
docker_log_has "--memory 6g" "og did not run under --memory 6g"
