#!/usr/bin/env bash
# Records the oracle of rust/crates/og-setup/tests/megane_kmers.rs: the .fai
# of each vector FASTA (samtools faidx, then twice.fa.fai by hand) and the
# .mk and .mi MEGAnE's own build_kmerset writes for it, in the pinned
# MEGANE_IMAGE (amd64; emulation is fine for inputs this small).
# Run from anywhere: tests/rust/vectors/megane-kmers/record.sh
set -euo pipefail
cd "$(dirname "$0")"
# shellcheck source=../../../../versions.env
. ../../../../versions.env
dk() { docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "${PWD}:/v" -w /v "$@"; }
dk --platform linux/amd64 "$SAMTOOLS_IMAGE" sh -c 'samtools faidx rules.fa && samtools faidx random.fa && samtools faidx twice.fa'
# twice.fa.fai lists y twice, before and after x, so its k-mers count twice.
awk -F'\t' '$1 == "y" { y = $0 } { l[NR] = $0 } END { print y; for (i = 1; i <= NR; i++) print l[i] }' twice.fa.fai > twice.fai.tmp
mv twice.fai.tmp twice.fa.fai
for n in rules random twice; do
  rm -rf out && mkdir out
  dk --platform linux/amd64 "$MEGANE_IMAGE" build_kmerset -fa "/v/${n}.fa" -prefix "$n" -outdir /v/out
  mv "out/${n}.mk" "out/${n}.mi" .
done
rm -rf out
shasum -a 256 ./*.fa ./*.fai ./*.mk ./*.mi
