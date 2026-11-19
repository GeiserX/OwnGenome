# shellcheck shell=sh
# shellcheck disable=SC2154  # YLEAF_DATA_VERSION and YLEAF_DATA_SHA256 come from versions.env, set by image-smoke.sh
# YLEAF_IMAGE row. Steps 37 and Y_HAPLOGROUP make Yleaf's pileup in
# SAMTOOLS_IMAGE at the marker positions of the data folder setup.sh
# --yleaf-data installs, and run Yleaf's own command line on it (-rp); the
# e2e case sv-mito-telomere-steps-3-y-haplogroup runs that on the fixture's
# chrY slice. This row fetches that archive (tag v<version>), checks its
# sha256 (YLEAF_DATA_SHA256), and checks what the steps rely on: the marker
# table and the tree in the archive are the ones this Yleaf predicts with,
# and this Yleaf still has the options, the pileup file and the samtools
# flags the steps use.
set -eu
v=$YLEAF_DATA_VERSION
t=$(mktemp -d)
# codeload is where github.com/archive/refs/tags/v<version>.tar.gz (setup.sh's
# URL) redirects; the image's busybox wget stops on github.com's long response
# headers. The sha256 below shows the bytes are the ones setup.sh checks.
wget -q -T 120 -O "$t/a.tgz" "https://codeload.github.com/genid/Yleaf/tar.gz/refs/tags/v${v}"
got=$(sha256sum "$t/a.tgz" | cut -d' ' -f1)
echo "archive sha256 ${got}, versions.env ${YLEAF_DATA_SHA256}"
if [ "$got" = "$YLEAF_DATA_SHA256" ]; then echo ok > /out/sha_ok.txt; else echo "got $got" > /out/sha_ok.txt; fi
d="Yleaf-${v}/yleaf/data"
tar -xzf "$t/a.tgz" -C "$t" "${d}/hg38/new_positions.txt" "${d}/hg_prediction_tables/tree.json"
n=$(grep -c . "$t/${d}/hg38/new_positions.txt")
echo "${n} marker positions, tree $(wc -c < "$t/${d}/hg_prediction_tables/tree.json") bytes"
echo "$n" > /out/positions_count.txt

# The installed Yleaf's own package folder (its data and Yleaf.py).
pkg=$(dirname "$(find / -path '*/site-packages/yleaf/Yleaf.py' 2>/dev/null | head -n 1)")
missing=""
for f in hg38/new_positions.txt hg_prediction_tables/tree.json; do
  cmp -s "$t/${d}/$f" "${pkg}/data/$f" || missing="${missing} data/${f} (differs from the archive)"
done
help=$(Yleaf -h 2>&1 || true)
for o in -rp --reuse_pileup --ref-fasta --no-update-check -rg -tree -force; do
  printf '%s\n' "$help" | grep -q -e " ${o}[ ,]" -e " ${o}\$" || missing="${missing} option ${o}"
done
for s in 'temp_haplogroup_pileup.pu' '-ABQ{quality_thresh}q1' 'Reusing existing pileup file' 'LOG.warning(f"Skipping'; do
  grep -qF -e "$s" "${pkg}/Yleaf.py" || missing="${missing} Yleaf.py: ${s}"
done
echo "missing:${missing:- none}"
if [ -z "$missing" ]; then echo ok > /out/launcher_contract.txt; else echo "missing:${missing}" > /out/launcher_contract.txt; fi
rm -rf "$t"
