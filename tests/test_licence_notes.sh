#!/usr/bin/env bash
# Check each licence notice where users meet it: setup's usage output, the docs
# for each resource, setup's install functions (before their first download),
# the chip and ORA scripts, and the source of og paraphase's help text.
# None of the calls below reaches Docker or the network: setup without a data
# directory prints usage, and the two scripts stop on their missing inputs.
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

read -r -d '' NOTICES <<'ROWS' || true
annotsv|AnnotSV's annotation bundle includes third-party data with their own terms.|00-reference-setup.md 05-annotsv.md
pcgr|PCGR/CPSR reference data includes third-party databases with their own terms.|00-reference-setup.md 17-cpsr.md
exomiser|Exomiser's hg38 and phenotype bundles include third-party data with their own terms.|00-reference-setup.md 48-exomiser.md
revel|The REVEL authors make the scores freely available for non-commercial use.|00-reference-setup.md 30-vcfanno.md hardware-requirements.md
ucsc|UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence.|00-reference-setup.md 43-archaic.md 44-ancient.md chip-data-guide.md
ora|Illumina's DRAGEN ORA Decompression Software requires acceptance of its end-user licence.|01-ora-to-fastq.md
paraphase|The PacBio Software License Agreement permits Paraphase use only with data generated on a PacBio instrument or provided by PacBio.|35-paralogs.md
ROWS
notice() { awk -F'|' -v k="$1" '$1 == k { print $2 }' <<< "$NOTICES"; }

rc=0
output=$(env -u GENOME_DIR bash "$repo/scripts/setup.sh" 2>&1) || rc=$?
if [ "$rc" -ne 1 ] || ! grep -q '^Usage:' <<< "$output"; then
  printf 'Expected setup usage (exit 1), got exit %s:\n%s\n' "$rc" "$output" >&2
  exit 1
fi

rows=0
while IFS='|' read -r _ notice pages; do
  rows=$((rows + 1))
  if ! grep -Fxq -- "  $notice" <<< "$output"; then
    printf 'Missing setup notice: %s\n' "$notice" >&2
    exit 1
  fi
  read -r -a docs <<< "$pages"
  for page in licences.md "${docs[@]}"; do
    if ! grep -Fq -- "$notice" "$repo/docs/$page"; then
      printf 'Missing notice in docs/%s: %s\n' "$page" "$notice" >&2
      exit 1
    fi
  done
done <<< "$NOTICES"
if [ "$rows" -ne 7 ]; then
  printf 'Expected 7 notices, found %s\n' "$rows" >&2
  exit 1
fi

# Each install function prints its notice before its first download.
installs=0
while IFS='|' read -r fn key var page; do
  installs=$((installs + 1))
  if ! grep -Fxq -- "${var}=\"$(notice "$key")\"" "$repo/scripts/setup.sh"; then
    printf 'setup.sh does not set %s to the %s notice\n' "$var" "$key" >&2
    exit 1
  fi
  if ! awk -v fn="$fn" -v want="echo \"  \${${var}} See docs/${page}.\"" '
    $1 == fn "()" { on = 1; next }
    on && /^}/ { exit }
    on && !seen && index($0, want) { seen = NR }
    on && !first && $0 !~ /^[[:space:]]*#/ && $0 ~ /(^|[^_[:alnum:]])fetch / { first = NR }
    END { exit !(seen && first && seen < first) }
  ' "$repo/scripts/setup.sh"; then
    printf 'setup.sh %s does not print the %s notice before its first fetch\n' "$fn" "$key" >&2
    exit 1
  fi
done <<'INSTALLS'
install_annotsv_bundle|annotsv|ANNOTSV_TERMS|05-annotsv.md
install_exomiser_data|exomiser|EXOMISER_TERMS|48-exomiser.md
install_ancient_panel|ucsc|UCSC_CHAIN_TERMS|44-ancient.md
install_archaic|ucsc|UCSC_CHAIN_TERMS|43-archaic.md
INSTALLS
if [ "$installs" -ne 4 ]; then
  printf 'Expected 4 install notices, found %s\n' "$installs" >&2
  exit 1
fi

# The chip converter and the ORA step print their notice before they check their inputs.
expect_run() {
  local name=$1 rc=$2 out=$3 line=$4
  if [ "$rc" -ne 1 ] || ! grep -Fxq -- "$line" <<< "$out"; then
    printf '%s (exit %s) did not print: %s\n%s\n' "$name" "$rc" "$line" "$out" >&2
    exit 1
  fi
}
rc=0
out=$(GENOME_DIR="$tmp" bash "$repo/scripts/chip-to-vcf.sh" S1 23andme 2>&1) || rc=$?
expect_run chip-to-vcf.sh "$rc" "$out" "$(notice ucsc) See docs/chip-data-guide.md."
rc=0
out=$(GENOME_DIR="$tmp" ORAD=/nonexistent bash "$repo/scripts/01-ora-to-fastq.sh" S1 /nonexistent /nonexistent.ora 2>&1) || rc=$?
expect_run 01-ora-to-fastq.sh "$rc" "$out" "$(notice ora) See docs/01-ora-to-fastq.md."

# og paraphase prints this notice on --help; the crate's unit tests check both help paths.
if ! grep -Fq -- "\"$(notice paraphase) See docs/35-paralogs.md.\"" "$repo/rust/crates/og-steps/src/paraphase.rs"; then
  printf 'og paraphase help does not carry the notice: %s\n' "$(notice paraphase)" >&2
  exit 1
fi

echo "Licence notices match setup output, install order, step output and documentation."
