#!/usr/bin/env bash
# Step 30 and Nextflow stage the files differently. Compare the track keys
# and all annotation options after normalizing the two file path templates.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

extract_annotations() {
  awk '
    BEGIN {
      basename["cadd_snv"] = "whole_genome_SNVs.tsv.gz"
      basename["cadd_indel"] = "gnomad.genomes.r4.0.indel.tsv.gz"
      basename["revel"] = "revel_grch38.tsv.gz"
      basename["alphamissense"] = "AlphaMissense_hg38.tsv.gz"
      quote = sprintf("%c", 39)
      for (track in basename) path[track] = "/genome/annotations/" basename[track]
      path["spliceai_snv"] = "/genome/annotations/" quote "\"$(basename \"$SPLICEAI_SNV\")\"" quote
      path["spliceai_indel"] = "/genome/annotations/" quote "\"$(basename \"$SPLICEAI_INDEL\")\"" quote
    }
    /^[[:space:]]*if \([a-z_]+\) \{$/ { track = $2; gsub(/[()]/, "", track) }
    /^[[:space:]]*[a-z_]+\)$/ { track = $1; sub(/\)$/, "", track) }
    /\[\[annotation\]\]$/ { print "[[annotation]]"; block = 1; count++; next }
    block && ($0 == "\"\"\"" || $0 == quote) { block = 0; next }
    block && /^file=/ {
      if ($0 == "file=\"${" track "}\"" || $0 == "file=\"" path[track] "\"") print "file=" track
      else print
      next
    }
    block && NF { print }
    END { if (!count || block) exit 1 }
  ' "$1"
}

extract_annotations "${1:-$ROOT/scripts/30-vcfanno.sh}" > "$TMP/bash.toml"
extract_annotations "${2:-$ROOT/modules/local/vcfanno/main.nf}" > "$TMP/nextflow.toml"
diff -u "$TMP/bash.toml" "$TMP/nextflow.toml"
echo "Step 30 and Nextflow configure the same vcfanno tracks and operations."
