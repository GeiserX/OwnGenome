#!/usr/bin/env bash
# Step 09 (ExpansionHunter) and step 09b (Stranger with its GRCh38 catalog).
# Stranger writes its output through a temporary name and annotates every
# ExpansionHunter locus it finds in the catalog.
. "$(dirname "$0")/lib.sh"

# The fixture reference holds 14 chromosomes, and ExpansionHunter stops on a
# catalog locus whose contig the reference lacks ("Invalid contig name 14"),
# in its reference region or in its off-target regions. The case keeps the bundled
# GRCh38 catalog's loci on the fixture's contigs and passes them as EH_CATALOG.
CAT="${GENOME_DIR}/reference/eh_catalog_fixture.json"
docker run --rm "$EXPANSIONHUNTER_IMAGE" cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json \
  > "${CASE_TMP}/catalog.json"
# The bundled catalog names contigs without "chr" (X, 14); ExpansionHunter
# matches them to chrX, chr14. It stops on the first one the reference lacks.
# shellcheck disable=SC2016  # jq program
jq --rawfile fai "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai" '
  ($fai | split("\n") | map(split("\t")[0] | select(. != null and . != "")) | map({(.): true}) | add) as $c
  | [.[] | select([(.ReferenceRegion, .OfftargetRegions) | select(. != null) | if type == "array" then .[] else . end
                   | split(":")[0] | $c[.] or $c["chr" + .]] | all)]' "${CASE_TMP}/catalog.json" > "$CAT"
N_CAT=$(jq length "$CAT" 2>/dev/null || echo 0)
echo "catalog loci: $(jq length "${CASE_TMP}/catalog.json"), on the fixture's contigs: ${N_CAT}"
check_ge "catalog loci on the fixture's contigs" "$N_CAT" 5

EH_CATALOG="$CAT" run_step 09-expansion-hunter.sh "$SAMPLE" male
check_step_exit 09-expansion-hunter.sh
check "the log names EH_CATALOG" has 'Variant catalog: .*eh_catalog_fixture\.json \(EH_CATALOG\)' "$(cat "$STEP_LOG")"
EH="${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
N_EH=$(grep -vc '^#' "${GENOME_DIR}/${EH}" 2>/dev/null || true)
check_ge "ExpansionHunter records" "$N_EH" 1

run_step 09b-stranger.sh "$SAMPLE"
check_step_exit 09b-stranger.sh
check "the log names the GRCh38 catalog" has 'variant_catalog_grch38' "$(cat "$STEP_LOG")"
OUT="${GENOME_DIR}/${SAMPLE}/expansion_hunter/${SAMPLE}_eh_stranger.vcf"
check "no temporary file is left" test ! -e "${OUT}.tmp"
check "the output declares STR_STATUS" grep -q '^##INFO=<ID=STR_STATUS' "$OUT"
check_eq "records in = records out" "$(grep -vc '^#' "$OUT" 2>/dev/null || true)" "$N_EH"
check_ge "records with the catalog's pathologic threshold" "$(grep -v '^#' "$OUT" 2>/dev/null | grep -c 'STR_PATHOLOGIC_MIN=' || true)" 3

finish
