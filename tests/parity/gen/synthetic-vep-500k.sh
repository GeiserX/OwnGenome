#!/usr/bin/env bash
set -euo pipefail
"${OG:?set OG to the og binary}" dev fixture parity-vep-500k vep-500k.vcf
