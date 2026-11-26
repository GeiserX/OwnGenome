#!/usr/bin/env bash
set -euo pipefail
"${OG:?set OG to the og binary}" dev fixture parity-clinvar synthetic.fa synthetic.fa.fai synthetic.vcf
