#!/usr/bin/env bash
# mei-filtered.sh: writes the MEGAnE VCFs of the mei-report cases that put a
# call MEGAnE flagged inside the FKTN window of tests/parity/inputs/mei/sites.tsv,
# in the current directory. scripts/ci/oracle-diff.sh runs it.
#   filtered.vcf       an Alu outside the catalog, and the FKTN SVA flagged LC
#   filtered-pass.vcf  the FKTN SVA flagged LC at 105639400, then a PASS one at
#                      105639660 in the same window
# og reads FILTER: a flagged call alone makes the site "call filtered out (LC)",
# and a PASS call in the window wins over a flagged one. The oracle,
# bin/mei_report.py, read any call in the window as seen and named the first,
# so both cases have a golden s.tsv (tests/parity/golden/mei-report/).
set -euo pipefail
export LC_ALL=C

vcf_head() {
  printf '##fileformat=VCFv4.1\n'
  printf '##source=MEI search version "v1.0.1 2022/04/12"\n'
  printf '##contig=<ID=chr20,length=64444167>\n'
  printf '##ALT=<ID=INS:ME,Description="Insertion of mobile element">\n'
  printf '##FILTER=<ID=LC,Description="Low confidence">\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
  printf '##INFO=<ID=MEPRED,Number=1,Type=String,Description="ME prediction status">\n'
  printf '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
}

# sva POS FILTER: an SVA_F call of one copy that MEGAnE identified.
sva() {
  printf 'chr9\t%s\t.\tN\t<INS:ME>\t.\t%s\tSVTYPE=Retroposon/SVA;MEPRED=PASS;HOMLEN=9;MEI_rpos=%s;MEI=SVA_F,1/1500,+/+;CN_conf=H\tCN\t1\n' \
    "$1" "$2" "$(($1 + 8))"
}

{
  vcf_head
  printf 'chr20\t5040478\t.\tN\t<INS:ME>\t.\tPASS\tSVTYPE=SINE/Alu;MEPRED=PASS;HOMLEN=9;MEI_rpos=5040486;MEI=AluYa5,0/288,+/+;CN_conf=H\tCN\t1\n'
  sva 105639660 LC
} > filtered.vcf

{
  vcf_head
  sva 105639400 LC
  sva 105639660 PASS
} > filtered-pass.vcf
