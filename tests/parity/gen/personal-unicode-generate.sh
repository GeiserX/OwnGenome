#!/usr/bin/env bash
# Rebuild the guard from SHA-256-pinned official UCD files. Network is only
# used here, never in the generator binary or the personal-data scanner.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
work=${1:?usage: personal-unicode-generate.sh SCRATCH_DIR}
mkdir -p "$work/ucd"
while IFS=$'\t' read -r hash file; do
  mkdir -p "$work/ucd/${file%/*}"
  curl --fail --silent --show-error --retry 2 "https://www.unicode.org/Public/${file%/*}/ucd/${file##*/}" -o "$work/ucd/$file"
  printf '%s  %s\n' "$hash" "$work/ucd/$file" | shasum -a 256 --check
 done < "$root/tests/parity/gen/personal-unicode-sources.tsv"
rustc --edition=2024 -D warnings "$root/tests/parity/gen/personal-unicode-generator.rs" -o "$work/generator"
"$work/generator" "$work/ucd" "$work/table.rs"
awk '/^\/\/ BEGIN GENERATED UNICODE GUARD$/ { copy=1; next }
     /^\/\/ END GENERATED UNICODE GUARD$/ { copy=0 }
     copy { print }' "$root/rust/crates/og-dev/src/personal_data.rs" > "$work/checked-in.rs"
cmp "$work/table.rs" "$work/checked-in.rs"
echo 'Unicode guard regenerates byte for byte.'
