#!/usr/bin/env bash
# summary-collect.sh: writes the sample folders of the summary-collect
# parity cases in the current directory, the same every time. Each folder
# is one case of og summary collect --only a against bin/collect_summary.py:
# the bash and the Nextflow names of each section's file, a stale file of
# each kind of step result and one written in the very second the run
# started, failed steps with and without a start time and with an empty
# one, unreadable files of each section, a corrupt gzip stream, the
# karyotypes, CRLF and lone-CR line ends, the three places a declared sex
# comes from, symlinked folders, and a VCF of 20,000 records.
# assets/ is a copy of the repository's callability catalog.
# scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)

put() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" > "$1"; }
putgz() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" | gzip -n -c > "$1"; }
at() { touch -t "$1" "${@:2}"; }   # [[CC]YY]MMDDhhmm[.SS], in UTC

T0=1577836800        # 2020-01-01 00:00:00 UTC, the start of the run
T0_TOUCH=202001010000.00
OLD=201906011230.00  # before the run
NEW=202103040506.07  # after it

mkdir -p assets
cp "${here}/../../../assets/callability_genes.tsv" assets/

COV_HEAD='chrom\tlength\tbases\tmean\tmin\tmax\n'
PED_HEAD='#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\n'
QC_MALE='key\tvalue\nsample\tS\ninferred_sex\tmale\nsex_check\tok\nx_depth_ratio\t0.51\ny_depth_ratio\t0.49\nx_sites\t400\nfreemix\t0.0010\n'
MOSAIC='mosaic_call\tnot_detected\nmosaic_depth_ratio\t0.98\nmosaic_depth_ratio_expected\t1\nmosaic_loss_depth\t0.02\nmosaic_loss_depth_low\t-0.01\nmosaic_loss_depth_high\t0.05\nmosaic_loss_af\tNA\nmosaic_floor\t0.1\nmosaic_x_het_sites\t0\n'
READS='key\tvalue\ninput_reads\t1000\nprimary_records\t1000\nmissing_reads\t0\nread_count_check\tok\n'

# The callability table of a whole catalog: every gene on chr1-22 and X
# (chrY left out, so its genes are not on this reference), a spread of
# shares around the 90% line, one gene of length 0.
callability() {
  awk -F'\t' 'BEGIN { OFS = "\t"; print "#chrom", "start", "end", "region", "10X" }
    /^#/ || $1 == "gene" || $4 == "." || $3 == "chrY" { next }
    { n++; len = $5 - $4; f = (n % 7) / 6; if (n == 5) { $5 = $4; len = 0 }
      print $3, $4, $5, $1, int(len * f) }' assets/callability_genes.tsv | gzip -n -c
}

# A VCF of N records: SNPs, indels, MNPs, symbolic and spanning alleles,
# non-PASS records, written as two gzip members, as BGZF is.
vcf() {
  local n=$1
  {
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
    awk -v n="$n" 'BEGIN { OFS = "\t"
      split("G|AT|A,<NON_REF>|<DEL>|*,C|.|GT,C|CG|T,TA", alts, "|")
      for (i = 1; i <= n; i++) {
        a = alts[i % 9 + 1]; r = (i % 9 == 7) ? "AT" : "A"
        f = (i % 5 == 0) ? "LowQual" : ((i % 11 == 0) ? "." : "PASS")
        print "chr" (i % 22 + 1), i, ".", r, a, 50, f, ".", "GT", "0/1"
      } }'
  } > vcf.tmp
  head -n $((n / 2)) vcf.tmp | gzip -n -c
  tail -n +$((n / 2 + 1)) vcf.tmp | gzip -n -c
  rm vcf.tmp
}

# --- full: the bash layout, every section, a run that started before them ----
put full/mosdepth/S.mosdepth.summary.txt "${COV_HEAD}chr1\t248956422\t7468692660\t30.00\t0\t255\ntotal\t3099734149\t96866692156\t31.25\t0\t255\n"
callability > full/mosdepth/S_callability.bed.gz
put full/indexcov/indexcov-indexcov.ped "${PED_HEAD}S\tS\t-9\t-9\t1\t-9\t1.01\t0.98\n"
put full/qc/S_sample_qc.tsv "${QC_MALE}declared_sex\tmale\n${MOSAIC}"
put full/qc/S_read_count.tsv "$READS"
mkdir -p full/vcf && vcf 20000 > full/vcf/S.vcf.gz
put full/logs/run_status.tsv "# run-all.sh\nmeta\tstarted_epoch\t${T0}\nmeta\tdeclared_sex\tmale\nstep\t16\tok\nstep\t16b\tok\nstep\t33\tok\n"
put full/run_manifest.tsv "# section\tkey\tvalue\nrun\tdeclared_sex\tmale\ndata\tclinvar_file_date\t2026-10-01\ndata\tcpsr_clinvar_release\t20250901\ndata\thla_database\tIPD-IMGT/HLA 3.60\textra\nshort\trow\n"

# --- nf: the Nextflow names, two total rows, two indexcov peds ---------------
put nf/coverage/S.mosdepth.summary.txt "${COV_HEAD}total\t1\t1\t30.5\t0\t9\ntotal\t1\t1\tnan\t0\t9\n"
callability > nf/coverage/S_callability.bed.gz
put nf/indexcov/B-indexcov.ped "${PED_HEAD}B\tB\t-9\t-9\t2\t-9\t2.0\t0.0\n"
put nf/indexcov/A-indexcov.ped "${PED_HEAD}A\tA\t-9\t-9\t2\t-9\t2.1\t0.02\n"
put nf/indexcov/.Z-indexcov.ped "${PED_HEAD}Z\tZ\t-9\t-9\t1\t-9\t1\t1\n"
put nf/qc/S_sample_qc.tsv "$QC_MALE"

# --- hidden: only a hidden ped, which a wildcard never matches ----------------
put hidden/indexcov/.x-indexcov.ped "${PED_HEAD}S\tS\t-9\t-9\t1\t-9\t1\t1\n"

# --- stale: files older than the run whose step was not ok in it ------------
put stale/mosdepth/S.mosdepth.summary.txt "${COV_HEAD}total\t1\t1\t29.9\t0\t9\n"
callability > stale/mosdepth/S_callability.bed.gz
put stale/indexcov/indexcov-indexcov.ped "${PED_HEAD}S\tS\t-9\t-9\t2\t-9\t2.02\t0.01\n"
put stale/qc/S_sample_qc.tsv "${QC_MALE}${MOSAIC}"
put stale/qc/S_read_count.tsv "$READS"
putgz stale/vcf/S.vcf.gz '#CHROM\nchr1\t1\t.\tA\tG\t1\tPASS\n'
put stale/logs/run_status.tsv "meta\tstarted_epoch\t${T0}\nstep\t16b\tskipped (opt-in)\nstep\t16\tfailed (exit 1)\nstep\t02\tskipped\nstep\t03\tskipped\nstep\t33\t\n"

# --- failed: no files, steps that failed, a start time as text only ----------
put failed/logs/run_status.tsv "meta\tstarted_utc\t2020-02-02 02:02 UTC\nstep\t16b\tfailed (exit 1)\nstep\t16\tfailed\nstep\t33\tfailed (exit 2)\nstep\t02\tfailed\nstep\t03\tok\n"
put failed-nodate/logs/run_status.tsv "meta\tstarted_epoch\tsoon\nstep\t33\tfailed (exit 1)\nstep\t16b\tfailed\n"
put failed-emptydate/logs/run_status.tsv "meta\tstarted_utc\t\nstep\t33\tfailed (exit 1)\n"

# --- unreadable files, one fault per section and folder ---------------------
put bad-a/mosdepth/S.mosdepth.summary.txt "${COV_HEAD}total\t1\t1\tabc\t0\t9\n"
putgz bad-a/mosdepth/S_callability.bed.gz 'chrom\tstart\tend\tregion\t10X\nchr1\t0\t10\tA\t5\n'
put bad-a/indexcov/indexcov-indexcov.ped '#family_id\tsample_id\nS\tS\n'
put bad-a/qc/S_sample_qc.tsv 'key\tvalue\nsex_check\tok\nmosaic_call\tmaybe\n'
put bad-a/qc/S_read_count.tsv 'key\tvalue\ninput_reads\t1.5\n'
mkdir -p bad-a/vcf && vcf 2000 > full.tmp && head -c 3000 full.tmp > bad-a/vcf/S.vcf.gz && rm full.tmp
# a byte flipped inside the deflate data, past the 10-byte gzip header
mkdir -p badgz/vcf && vcf 2000 > badgz.tmp
printf '\377' | dd of=badgz.tmp bs=1 seek=40 conv=notrunc 2>/dev/null
mv badgz.tmp badgz/vcf/S.vcf.gz

putgz bad-b/mosdepth/S_callability.bed.gz '#chrom\tstart\tend\tregion\t10X\nchr1\t0\t10\tA\n'
put bad-b/indexcov/indexcov-indexcov.ped "$PED_HEAD"
put bad-b/qc/S_sample_qc.tsv 'key\tvalue\ninferred_sex\tfemale\nshort\nkaryotype\twritten\nx_depth_ratio\t1.0\ny_depth_ratio\t0\nx_sites\tinf\nmosaic_call\tgain\nmosaic_floor\tinf\nmosaic_loss_af\tnan\nmosaic_depth_ratio\t5\nmosaic_x_het_sites\t12a\nmosaic_extra\n'
put bad-b/qc/S_read_count.tsv 'key\tvalue\ninput_reads\t10\nprimary_records\t10\nmissing_reads\t0\nread_count_check\tunknown\n'

putgz bad-c/mosdepth/S_callability.bed.gz '#chrom\tstart\tend\tregion\t10X\nchr1\t0\t10\tNOT_A_GENE\t5\n'
put bad-c/qc/S_read_count.tsv 'key\tvalue\ninput_reads\t10\nprimary_records\t-\nmissing_reads\t0\n'
put bad-c/mosdepth/S.mosdepth.summary.txt "${COV_HEAD}chr1\t1\t1\t30\t0\t9\n"

putgz bad-d/mosdepth/S_callability.bed.gz '#chrom\tstart\tend\tregion\t10X\nchr1\tx\t10\tA1BG\t5\n'
put bad-d/qc/S_read_count.tsv 'key\tvalue\ninput_reads\t007\nprimary_records\t+7\nmissing_reads\t-0\nread_count_check\tmismatch\nnote\n'
put bad-d/mosdepth/S.mosdepth.summary.txt "${COV_HEAD}total\t1\t1\n"

putgz bad-e/mosdepth/S_callability.bed.gz '#chrom\tstart\tend\tregion\t10X\n\n'

# --- karyotypes from indexcov's copy numbers and somalier's depth ratios ------
k=0
for xy in "1 0" "2.2 0.8" "2.5 0" "1.0 1.0" "3 0" "6 0" "inf 0" "nan 1" "1.29 2.31" "-0.2 0"; do
  read -r x y <<< "$xy"
  k=$((k + 1))
  put "kary-${k}/indexcov/indexcov-indexcov.ped" "${PED_HEAD}S\tS\t-9\t-9\t1\t-9\t${x}\t${y}\n"
  put "kary-${k}/qc/S_sample_qc.tsv" "key\tvalue\ninferred_sex\tmale\nx_depth_ratio\t$(awk -v v="$x" 'BEGIN { printf "%.3f", v / 2 }')\ny_depth_ratio\t$(awk -v v="$y" 'BEGIN { printf "%.3f", v / 2 }')\nx_sites\t$([ "$k" -eq 3 ] && echo 10 || echo 11)\n"
done

# --- line ends: CRLF, lone CR, blank lines -----------------------------------
put crlf/indexcov/indexcov-indexcov.ped "#family_id\tsample_id\tsex\tCNchrX\tCNchrY\r\n\r\nS\tS\t2\t2.0\t0.0\r\n"
put crlf/qc/S_sample_qc.tsv 'key\tvalue\r\ninferred_sex\tfemale\r\nx_sites\t5\r\n'
put crlf/logs/run_status.tsv "meta\tstarted_epoch\t${T0}\rstep\t16\tskipped\rmeta\tdeclared_sex\tfemale\r"
put crlf/run_manifest.tsv 'data\tclinvar_file_date\t2025-01-01\r\n'
putgz crlf/vcf/S.vcf.gz '#h\r\nchr1\t1\t.\tA\tG\t1\tPASS\r\nchr1\t2\t.\tA\tAT\t1\tPASS\rchr1\t3\t.\tA\t\xc3\xa9\t1\tPASS\nchr1\t4\t.\t\xff\tG\t1\t.\n'

# --- the declared sex: the run status, the manifest, the step 33 table ------
put declared-status/logs/run_status.tsv 'meta\tdeclared_sex\t\nmeta\tdeclared_sex\tfemale\n'
put declared-status/run_manifest.tsv 'run\tdeclared_sex\tmale\n'
put declared-manifest/run_manifest.tsv 'run\tdeclared_sex\t\nrun\tdeclared_sex\tmale\n'
put declared-manifest/qc/S_sample_qc.tsv 'key\tvalue\ninferred_sex\tmale\ndeclared_sex\tfemale\n'
put declared-qc/qc/S_sample_qc.tsv 'key\tvalue\ninferred_sex\tmale\ndeclared_sex\tfemale\n'

# --- symlinks, as the HTML_REPORT task stages its inputs ---------------------
put store/qc/S_read_count.tsv "$READS"
putgz store/S.vcf.gz '#h\nchr1\t1\t.\tA\tG\t1\tPASS\n'
mkdir -p links/vcf
ln -s ../store/qc links/qc
ln -s ../../store/S.vcf.gz links/vcf/S.vcf.gz
ln -s nowhere links/vcf/broken.vcf.gz

# --- not gzip: a plain VCF named .vcf.gz ---------------------------------------
put notgz/vcf/S.vcf.gz '##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t1\tPASS\n'

mkdir -p empty

# Every file after the run started, then the ones that must be older, and
# the one written in the very second the run started.
find . -type f -exec touch -t "$NEW" {} +
at "$OLD" stale/mosdepth/S_callability.bed.gz stale/indexcov/indexcov-indexcov.ped \
  stale/qc/S_sample_qc.tsv stale/qc/S_read_count.tsv stale/vcf/S.vcf.gz
at "$T0_TOUCH" stale/mosdepth/S.mosdepth.summary.txt
