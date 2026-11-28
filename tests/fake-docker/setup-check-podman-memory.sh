#!/usr/bin/env bash
# The check reads the engine's memory from CONTAINER_ENGINE. podman answers
# {{.MemTotal}} with an error and keeps the memory in {{.Host.MemTotal}}:
# og setup check asks the second when the first gives no number. The bash
# checks (OG_SETUP=bash) asked a literal docker, here the fake one's 16 GiB.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"

mkdir -p "${CASE_WORK}/podman-bin"
cat > "${CASE_WORK}/podman-bin/podman" <<'PODMAN'
#!/usr/bin/env bash
# A podman in front of the fake docker: its own memory answers, the rest
# passed on (and logged) as docker calls.
if [ "$*" = "info --format {{.MemTotal}}" ]; then
  echo "podman: template: info:1:2: executing \"info\" at <.MemTotal>: can't evaluate field MemTotal" >&2
  exit 125
fi
if [ "$*" = "info --format {{.Host.MemTotal}}" ]; then echo 8589934592; exit 0; fi
exec docker "$@"
PODMAN
chmod +x "${CASE_WORK}/podman-bin/podman"
export PATH="${CASE_WORK}/podman-bin:${PATH}" CONTAINER_ENGINE=podman

run_expect 0 podman-og "${SCRIPTS}/validate-setup.sh"
output_has podman-og '\[WARN\].*Docker memory: 8\.0 GB — 16 GB\+ recommended'
output_lacks podman-og 'Docker memory: 16\.0 GB'

OG_SETUP=bash run_expect 0 podman-bash "${SCRIPTS}/validate-setup.sh"
output_has podman-bash '\[OK\].*Docker memory: 16\.0 GB'
