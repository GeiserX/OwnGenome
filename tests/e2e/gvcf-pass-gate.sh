#!/usr/bin/env bash
# Synthetic controls for the actual bcftools view commands in both callers.
# No HG002 outputs or downloaded panels are needed. Stream the fixture into
# Docker so this case also works where Docker cannot mount the checkout.
set -euo pipefail
. "$(dirname "$0")/lib.sh"

mkdir -p "${CASE_TMP}/input"
INPUT="${CASE_TMP}/input"
COMMANDS="${INPUT}/commands.tsv"
: > "$COMMANDS"
for step in 26b-ancestry-fine 44-ancient 25-prs; do
  case "$step" in
    26b-*) module=ancestry_fine ;;
    44-*) module=ancient ;;
    25-*) module='prs' ;;
  esac
  for source in "scripts/${step}.sh" "modules/local/${module}/main.nf"; do
    # Keep each view's options verbatim. Only replace its panel path, remove
    # its filename (read stdin instead), and decode Groovy's escaped quotes.
    awk -v source="$source" '
      /bcftools view .*-[T] / {
        line = $0
        sub(/^.*bcftools view /, "bcftools view ", line)
        if (line !~ /-O[uv]/) next
        mode = line ~ /-Ou/ ? "-Ou" : "-Ov"
        sub(/-T [^ ]+/, "-T sites.tsv", line)
        sub(/-O[uv].*$/, mode, line)
        if (source ~ /^modules\//) gsub(/\\"/, "\"", line)
        kind = line ~ /--trim-alt-alleles/ ? "gvcf" : "vcf"
        print source "\t" kind "\t" line
      }
    ' "${REPO}/${source}" >> "$COMMANDS"
  done
done
check_eq "all ten extraction branches are exercised" "$(wc -l < "$COMMANDS" | tr -d ' ')" 10

cat > "${INPUT}/ref.fa" <<'EOF'
>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF
printf 'chr1\t100\t6\t100\t101\n' > "${INPUT}/ref.fa.fai"
for pos in 10 20 30 40 50 60 70 80 90; do
  printf 'chr1\t%s\n' "$pos"
done > "${INPUT}/sites.tsv"
cat > "${INPUT}/sample.g.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FILTER=<ID=LowQual,Description="Failed variant">
##FILTER=<ID=RefCall,Description="Filtered reference call">
##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SYNTHETIC
chr1	10	.	A	G,<*>	40	PASS	.	GT:GQ:DP	0/1:40:30
chr1	20	.	A	G,<*>	1	LowQual	.	GT:GQ:DP	1/1:1:2
chr1	30	.	A	<*>	0	RefCall	.	GT:GQ:DP	0/0:40:30
chr1	35	.	A	<*>	.	.	END=45	GT:GQ:DP	0/0:40:30
chr1	50	.	A	G,<*>	.	.	.	GT:GQ:DP	0/1:20:10
chr1	55	.	A	<*>	.	.	END=65	GT:GQ:DP	./.:0:0
chr1	70	.	A	G,<*>	40	PASS	.	GT:GQ:DP	0/1:19:1
chr1	75	.	A	<*>	.	.	END=85	GT:GQ:DP	0/0:0:0
chr1	90	.	A	G,<*>	40	PASS	.	GT:DP	0/1:30
EOF
cat > "${INPUT}/run.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
bcftools view -Oz -o sample.g.vcf.gz sample.g.vcf
bcftools index -t sample.g.vcf.gz
bcftools convert --gvcf2vcf -f ref.fa -R sites.tsv -Ob -o expanded.bcf sample.g.vcf.gz
bcftools view -e 'INFO/END>0' -Ob -o variants.bcf sample.g.vcf.gz
while IFS=$'\t' read -r source kind command; do
  input=expanded.bcf
  [ "$kind" = gvcf ] || input=variants.bcf
  printf 'RESULT\t%s\t%s\t' "$source" "$kind"
  bash -euo pipefail -c "$command" < "$input" \
    | bcftools query -f '%POS,' | tr -d '\n'
  printf '\n'
done < commands.tsv
EOF

LOG="${CASE_TMP}/extract.log"
if ! tar -C "$INPUT" -cf - . \
  | docker run --rm -i --network none --entrypoint bash "$BCFTOOLS_IMAGE" \
      -euo pipefail -c 'mkdir /w; cd /w; tar -xf -; bash run.sh' > "$LOG" 2>&1; then
  cat "$LOG"
  fail "the real bcftools extraction completed"
else
  cat "$LOG"
  while IFS=$'\t' read -r source kind _command; do
    case "$source:$kind" in
      *ancient*) expected=10,40,50, ;;
      *:gvcf) expected=10,40,50,70,80,90, ;;
      *:vcf) expected=10,50,70,90, ;;
    esac
    got=$(awk -F'\t' -v src="$source" -v k="$kind" '$1 == "RESULT" && $2 == src && $3 == k {print $4}' "$LOG")
    check_eq "$source ($kind) keeps exactly the eligible panel sites" "$got" "$expected"
  done < "$COMMANDS"
fi
finish
