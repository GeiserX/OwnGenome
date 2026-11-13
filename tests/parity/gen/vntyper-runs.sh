#!/usr/bin/env bash
# The VNtyper run folders of the vntyper-report cases, one folder per case
# (kestrel/kestrel_result.tsv after its ## banner, coverage/coverage_summary.tsv).
set -euo pipefail
python3 - <<'PY'
import os
runs = {
    'pos': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tHigh_Precision*\tNot flagged\t3\t59dupC\tCCCCACC\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'low': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tLow_Precision\tNot flagged\t3\t59dupC\tCCCCACC\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'flag': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tHigh_Precision*\tFalse_Positive_4bp_Insertion\t3\t59dupC\tCCCCACC\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'holes': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tHigh_Precision*\tNot flagged\t3\t59dupC\tCCCCACC\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t62.00\tFAIL\n'),
    'nounc': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tHigh_Precision*\tNot flagged\t3\t59dupC\tCCCCACC\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\tNA\tPASS\n'),
    'neg': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n147.94\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'red': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n70.00\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tREDUCED\n'),
    'thin': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n20.00\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tFAIL\n'),
    'na': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\nNA\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tNOT_EVALUATED\n'),
    'empty': (None, 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'hdr': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n147.94\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'negmean': ('##VNtyper Kestrel result\n##source=test\nMotif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\nNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n-5.00\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
    'crlf': ('##VNtyper Kestrel result\n##source=test\nMotifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\r\n2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\tHigh_Precision*\tNot flagged\t3\t59dupC\tCCCCACC\r\n', 'mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n'),
}
for name, (kestrel, cov) in runs.items():
    os.makedirs(f'{name}/kestrel')
    os.makedirs(f'{name}/coverage')
    if kestrel is not None:
        with open(f'{name}/kestrel/kestrel_result.tsv', 'w', newline='') as f:
            f.write(kestrel)
    with open(f'{name}/coverage/coverage_summary.tsv', 'w') as f:
        f.write(cov)
os.makedirs('nok')
os.makedirs('bad/kestrel')
os.makedirs('bad/coverage')
with open('bad/kestrel/kestrel_result.tsv', 'wb') as f:
    f.write(b'\xff\xfe\x00bad')
with open('bad/coverage/coverage_summary.tsv', 'w') as f:
    f.write('mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n1075.32\t1000.00\t300.00\t95\t1942\t4500\t0\t0.00\tPASS\n')
with open('empty.log', 'w') as f:
    f.write("ValueError: FASTQ layout 'empty' cannot be consumed without dropping reads.\n")
PY
