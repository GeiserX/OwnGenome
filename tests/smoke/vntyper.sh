# shellcheck shell=sh
# VNTYPER_IMAGE row: what step 52 runs, on VNtyper's own public test data
# (tests/smoke/vntyper_testdata.sh fetches two BAMs out of its Zenodo archive):
# a sample with a MUC1 VNTR frameshift and one without. The image's samtools
# indexes them, then VNtyper runs with the step's flags (and no IGV view in its HTML)
# from /opt/vntyper,
# the folder its config.json paths are relative to.
set -e
sh /smoke/vntyper_testdata.sh /out/data positive.bam negative.bam
# Hashing stdin must also work when the destination contains a newline.
newline_out='/out/new
line'
sh /smoke/vntyper_testdata.sh "$newline_out" negative.bam
cmp /out/data/negative.bam "$newline_out/negative.bam"
rm -rf "$newline_out"
for s in positive negative; do
  samtools index "/out/data/${s}.bam"
  (cd /opt/vntyper && vntyper pipeline --bam "/out/data/${s}.bam" -o "/out/${s}" --sample-name "$s" \
     --threads 2 --fast-mode --reference-assembly hg38 --report-igv off) > "/out/${s}.log" 2>&1 \
    || { tail -n 40 "/out/${s}.log"; exit 1; }
  echo "== ${s}: kestrel/kestrel_result.tsv"; cat "/out/${s}/kestrel/kestrel_result.tsv"
  echo "== ${s}: coverage/coverage_summary.tsv"; cat "/out/${s}/coverage/coverage_summary.tsv"
done
# A BAM with the same header and no read over MUC1, as a genome slice without
# chr1:155 Mb gives: VNtyper stops with exit 1 and names the empty read
# layout, which og vntyper report turns into "no read over MUC1".
samtools view -b -o /out/data/empty.bam /out/data/negative.bam chr1:1-1000 && samtools index /out/data/empty.bam
rc=0
(cd /opt/vntyper && vntyper pipeline --bam /out/data/empty.bam -o /out/empty --sample-name empty \
   --threads 2 --fast-mode --reference-assembly hg38 --report-igv off) > /out/empty.log 2>&1 || rc=$?
echo "$rc" > empty_rc.txt
echo "== empty: exit ${rc}"; grep -m1 "FASTQ layout" /out/empty.log || true
vntyper --version > version.txt 2>&1
# calls.tsv: sample, Kestrel's confidence, the variant's name, the mean depth
# over the VNTR and VNtyper's coverage verdict. Each file's first row under its
# header, by column name, as csv.DictReader read it: Kestrel's "#" banner lines
# are skipped, and so are blank lines under the header.
: > /out/calls.tsv
for s in positive negative; do
  awk -F'\t' -v s="$s" '
    function get(f, k) { return ((f, k) in col) ? val[f, col[f, k]] : "" }
    { f = (FILENAME == ARGV[1]) ? 1 : 2 }
    (f == 1 && /^#/) || ($0 == "" && (f in hdr)) { next }
    !(f in hdr) { for (i = 1; i <= NF; i++) col[f, $i] = i; hdr[f] = 1; next }
    seen[f]++ { next }
    { for (i = 1; i <= NF; i++) val[f, i] = $i }
    END {
      if (!seen[1] || !seen[2]) { print s ": no result row in " (seen[1] ? ARGV[2] : ARGV[1]) > "/dev/stderr"; exit 1 }
      n = get(1, "Nomenclature"); if (n == "") n = get(1, "Variant")
      print s "\t" get(1, "Confidence") "\t" n "\t" get(2, "mean") "\t" get(2, "coverage_qc")
    }' "/out/${s}/kestrel/kestrel_result.tsv" "/out/${s}/coverage/coverage_summary.tsv" >> /out/calls.tsv
done
cat /out/calls.tsv
