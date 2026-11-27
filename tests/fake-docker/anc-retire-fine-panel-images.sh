#!/usr/bin/env bash
# setup.sh --ancestry-fine-panel pulls the images step 26b runs and no other:
# og, Beagle, and bcftools with CA certificates for the install. RFMix and
# ADMIXTURE left step 26b when og took over the painting and the fractions,
# so neither is pulled, named in the install log or pinned in versions.env.
# Every image reads as not pulled yet, so each one the install wants is
# pulled; then rust/og.sh reports og as not built, which stops the install
# right after the pulls, before any download. OG names bin/og for the
# downloads (og setup fetch), so only need_og sees the og.sh that fails.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

COPY="${CASE_WORK}/checkout"
mkdir -p "${COPY}/scripts" "${COPY}/rust"
ln -s "${REPO_ROOT}/scripts/lib" "${COPY}/scripts/lib"
ln -s "${REPO_ROOT}/scripts/ci" "${COPY}/scripts/ci"
ln -s "${REPO_ROOT}/bin" "${COPY}/bin"
ln -s "${REPO_ROOT}/assets" "${COPY}/assets"
cp "${SCRIPTS}/setup.sh" "${COPY}/scripts/setup.sh"
cp "${REPO_ROOT}/versions.env" "${COPY}/versions.env"
printf '#!/usr/bin/env bash\nexit 1\n' > "${COPY}/rust/og.sh"
chmod +x "${COPY}/rust/og.sh"

G="${CASE_WORK}/genome"
mkdir -p "$G"
export FAKE_DOCKER_MISSING_IMAGES='.' FETCH_TRIES=1 FETCH_WAIT=0 OG="${REPO_ROOT}/bin/og"

run_rc fine "${COPY}/scripts/setup.sh" --ancestry-fine-panel "$G"
expect_rc fine 1
output_has fine '^\[WARN\] Could not build og'

# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
for v in OG_IMAGE BCFTOOLS_HTTPS_IMAGE BEAGLE_IMAGE; do
  grep -qF "pull :: pull ${!v} " "$FAKE_DOCKER_LOG" || fail "setup.sh --ancestry-fine-panel did not pull ${v} (${!v})"
done
if hit=$(grep -inE 'rfmix|admixture' "$FAKE_DOCKER_LOG"); then
  fail "setup.sh --ancestry-fine-panel called docker with RFMix or ADMIXTURE: ${hit}"
fi
output_lacks fine '[Rr][Ff][Mm]ix|ADMIXTURE|[Aa]dmixture'
if hit=$(grep -nE '^(RFMIX|ADMIXTURE)_IMAGE=' "${REPO_ROOT}/versions.env"); then
  fail "versions.env still pins ${hit}"
fi
echo "PASS setup.sh --ancestry-fine-panel pulls og, Beagle and bcftools, and neither RFMix nor ADMIXTURE"
