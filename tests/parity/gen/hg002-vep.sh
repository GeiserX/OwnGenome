#!/usr/bin/env bash
# hg002-vep.sh: writes hg002_vep.vcf in the current directory: the public
# HG002 VEP slice og-vcf's tests read (rust/crates/og-vcf/tests/data, VEP
# --everything output with every CSQ field), decompressed, for the
# og vcf csq-rank case hg002-vep.
set -euo pipefail
gzip -dc "$(dirname "$0")/../../../rust/crates/og-vcf/tests/data/hg002_vep.vcf.gz" > hg002_vep.vcf
