#!/bin/sh
# Fetch members of VNtyper's public test data without the whole archive.
#
# VNtyper publishes its test BAMs as one 1.2 GB zip on Zenodo (record 19181821,
# "VNtyper 2.0 test data", the archive_file of tests/test_data_config.json in
# VNtyper v2.0.47). A Zenodo record never changes, so where each member's
# deflated bytes sit in the archive is fixed: MEMBERS below pins it, read once
# from the archive's central directory and the member's local header. This
# script fetches only those bytes with one HTTP range request per member,
# about 10 MB instead of 1.2 GB (the .bai files are not fetched: samtools
# index writes them). It wraps them in a gzip header and a trailer holding the
# CRC-32 and size the archive stores, so gzip inflates them and checks both,
# and then checks the result against the pinned sha256. A wrong offset, a
# changed archive or a cut transfer fails one of those checks.
#
# Usage: vntyper_testdata.sh OUTDIR NAME...   (NAME: positive.bam or negative.bam)
# POSIX sh with curl, gzip and sha256sum: it runs in VNTYPER_IMAGE and on the host.
set -eu

URL="https://zenodo.org/records/19181821/files/data.zip?download=1"

# member NAME: "path data_offset deflated_size size crc32 sha256" of NAME.
member() {
  case $1 in
    # Positive: a MUC1 VNTR frameshift VNtyper calls High_Precision* (its
    # integration test example_dfc3_GRCh38_fast runs the same reads on NCBI
    # contig names). Reads realigned by VNtyper to GRCh38 chr1, UCSC names.
    positive.bam) echo "data/remapped/bwa/hg38/example_dfc3_hg38_bwa.bam 1023038956 6194566 6194149 f8fa50d2 5934330a8ac2144ac0cc6ae8c84d7d9262e90b28127417e4068f7d662b7c64bb" ;;
    # Negative: no variant called (its integration test
    # example_40cf_hg38_subset_fast_gdp_guard), mean VNTR depth about 148x.
    negative.bam) echo "data/example_40cf_hg38_subset.bam 154 3376349 3378669 286d97d5 44aff8e3b2096b7ea3529989660052c61b29c8a26afd7f48de4988f7047d0652" ;;
    *) return 1 ;;
  esac
}

# le32 N: N as four little-endian bytes.
le32() {
  printf '%b' "$(printf '\\0%03o\\0%03o\\0%03o\\0%03o' $(($1 & 255)) $(($1 >> 8 & 255)) $(($1 >> 16 & 255)) $(($1 >> 24 & 255)))"
}

[ $# -ge 2 ] || { echo "usage: vntyper_testdata.sh OUTDIR NAME... (NAME of positive.bam, negative.bam)" >&2; exit 2; }
out=$1
shift
for name in "$@"; do
  member "$name" > /dev/null || { echo "vntyper_testdata.sh: unknown member ${name} (positive.bam, negative.bam)" >&2; exit 2; }
done
mkdir -p "$out"
for name in "$@"; do
  read -r path off csize size crc want <<EOF
$(member "$name")
EOF
  part="${out}/${name}.part"
  curl -fsSL --retry 4 --retry-all-errors --retry-delay 15 -r "${off}-$((off + csize - 1))" -o "${part}.raw" "$URL"
  got=$(wc -c < "${part}.raw" | tr -d ' ')
  [ "$got" -eq "$csize" ] || { echo "${path}: got ${got} bytes, want ${csize} (the server ignored the range?)" >&2; exit 1; }
  # gzip member: magic, deflate, no flags, no time, no extra flags, unknown OS;
  # then the raw deflate stream; then CRC-32 and size, both little-endian.
  { printf '\037\213\010\000\000\000\000\000\000\377'; cat "${part}.raw"; le32 "$((0x$crc))"; le32 "$size"; } > "${part}.gz"
  gzip -dc "${part}.gz" > "$part" || { echo "${path}: size or CRC-32 differs from the archive's directory" >&2; exit 1; }
  got=$(sha256sum "$part" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "${path}: sha256 ${got}, pinned ${want}" >&2; exit 1; }
  rm -f "${part}.raw" "${part}.gz"
  mv "$part" "${out}/${name}"
  printf '%s\t%s\t%s\t%s\n' "$name" "$path" "$size" "$got"
done
