#!/usr/bin/env bash
set -euo pipefail
"${OG:?set OG to the og binary}" dev fixture parity-consumer synthetic.g.vcf.gz synthetic-sites.tsv.gz
