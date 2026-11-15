#!/usr/bin/env bash
# Step 07 (PharmCAT) writes a report.json that parses and calls at least one gene
# (the fixture covers CYP2C19, CYP2C9 and CYP2D6).
. "$(dirname "$0")/lib.sh"

run_step 07-pharmacogenomics.sh "$SAMPLE"
check_step_exit 07-pharmacogenomics.sh

JSON="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.report.json"
check "report.json exists" test -s "$JSON"
check_pharmcat_called "genes called in report.json" "$JSON"
check "HTML report exists" nonempty "${SAMPLE}/vcf/${SAMPLE}.report.html"

finish
