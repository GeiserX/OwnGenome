#!/usr/bin/env python3
"""lpa_kiv2.py: the LPA KIV-2 repeat count from short-read depth (step 46).

Lp(a) is a mostly inherited cardiovascular risk factor, and the number of
kringle IV type 2 (KIV-2) repeats in LPA sets the size of apo(a): fewer
repeats, a smaller isoform and, on average, a higher Lp(a) level. The KIV-2
repeat is a 5.5 kb unit with two exons; a Danish population study counted 6
to 99 on both chromosomes together (Kamstrup et al. 2009, PMID 19509380).
GRCh38 holds six units, and their exons are near-identical, so an aligner
that reports every read places each KIV-2 read on one of the six (with
mapping quality 0). The depth over those exons, against the depth over
single-copy LPA exons, counts the units the sample carries: the approach of
the published short-read estimators (Zekavat et al. 2018, PMID 29973585;
Behera et al. 2024, PMID 39449055, which also gives the GRCh38 region and the
six reference units).

The regions and the method are the LPA_KIV2 row of bin/cn_loci.tsv, run by
bin/copy_number.py (step 46 and the LPA processes call that engine); this
file lays out the step's table and keeps its own command line.

  lpa_kiv2.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  lpa_kiv2.py estimate --bedcov BEDCOV.tsv --bedcov-mapq1 BEDCOV_Q1.tsv --out LPA.tsv
      read `samtools bedcov -j` over those regions (chrom, start, end, name,
      summed depth; all MAPQ, duplicates left out) and `samtools bedcov -j -Q 1`
      over the same regions, and write LPA.tsv (metric<TAB>value). Exit 0
      whatever the status; 2 on unreadable input. Two empty tables mean the
      BAM has no chr6: status 'not assessed'.
  lpa_kiv2.py null-allele --pileup PILEUP.tsv --table LPA.tsv --out NULL.tsv
      read `samtools mpileup` (no reference, every MAPQ) at the seven copies
      of the KIV-2 4925G>A site (`traits_report.py --marker-regions
      LPA_KIV2_4925`) and write its reading to NULL.tsv (metric<TAB>value),
      with the count N of LPA.tsv (written by `estimate` or `copy_number.py
      call --locus LPA_KIV2`, left as it is). Exit 0 whatever the reading (an
      unreadable pileup reads 'not assessed'); 2 when NULL.tsv cannot be
      written.

The regions are exons of the MANE Select transcript ENST00000316300.10,
numbered from its first exon:
  kiv2_e03 to kiv2_e16  the six KIV-2 units of GRCh38 (exons 4 to 15), plus
                        exon 3 (the last exon of KIV-1) and exon 16 (the first
                        of KIV-3), which are near-identical to KIV-2 exons and
                        take some of its reads. Seven 160 bp and seven 182 bp
                        exons, 2,394 bp: six KIV-2 units plus one single-copy
                        pair per chromosome.
  flank_e19 to flank_e33  single-copy kringle exons downstream (the second
                        exon of KIV-4 to KV), under 90% identical to any KIV-2
                        exon. Exons 17 and 18 (about 95% identical) are in
                        neither group, nor are exons 1 and 2 or the protease
                        exons, which resemble plasminogen (PLG).

With D the mean depth over the flank exons (two copies) and M the mean depth
over the KIV-2 block, the block holds N + 2 copies of a 342 bp pair spread
over seven reference pairs, so M = (N + 2) x D / 14 and
    N = 14 x M / D - 2      (the total over both chromosomes)
The 95% interval comes from the spread of the per-exon depths in each group:
relative standard error sqrt(cv_kiv2^2 / 14 + cv_flank^2 / 15), applied to
N + 2 on the log scale. It covers sampling noise only, not a systematic bias
of an aligner or a library.

status is 'estimated', 'too little depth' or 'not assessed'; the message says
why.

The KIV-2 4925G>A splice variant (Coassin et al. 2017, PMID 28444229) sits at
the donor splice site of the second KIV-2 exon: carriers had a much lower
Lp(a) than non-carriers with the same isoform size, most of all with small
isoforms, and it is usually in one or two of a person's KIV-2 repeats. GRCh38
holds the site in each of its six units and in KIV-1, so null-allele adds up
the reads at all seven (bin/traits_report.py MARKERS) and reads the share
with the variant as a fraction of copies: status is 'present', 'not seen',
'uncertain', 'too little depth' or 'not assessed', and with a count, copies
is the fraction x (N + 2), the copies that carry it.
A short-read fraction cannot say which chromosome or which repeat carries
it, so it cannot say whether it sits on the small or the large isoform. Research estimate, not a diagnostic test: the actionable test is a blood
Lp(a) measurement. Standard library only: it runs in the Python image and
under tests/test_lpa_kiv2.py.
"""
import argparse
import os
import sys

import copy_number as cn

LOCUS = "LPA_KIV2"
_ROW = []
ADVICE = ("A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 "
          "European Atherosclerosis Society consensus advises measuring it at least once in every adult.")


def row():
    """The LPA_KIV2 row of the catalog, read once."""
    if not _ROW:
        _ROW.append(cn.load_catalog()[LOCUS])
    return _ROW[0]


def __getattr__(name):
    # The regions and thresholds come from the catalog row, read on first use.
    r = row()
    if name == "REGIONS":   # in catalog order: the KIV-2 block, then the single-copy exons
        return [(r["chrom"], s, e, n) for n, s, e, _ in r["_targets"] + r["_controls"]]
    vals = {"MIN_FLANK_DEPTH": r["_min_control"], "MIN_MAPQ0_SHARE": r["_mapq0_min"],
            "LOW": r["_call"][1][0], "HIGH": r["_call"][1][1], "READINGS": r["_readings"]}
    if name in vals:
        return vals[name]
    raise AttributeError(name)


def bed_lines():
    return ["\t".join(map(str, r)) for r in cn.regions(row())]


def read_bedcov(path):
    """{name: mean depth} from `samtools bedcov` output over the regions,
    which must all be there and nothing else; {} for an empty table (no chr6)."""
    return cn.read_bedcov(path, row())


def table(b):
    """([(metric, value)], summary line) of the engine's LPA_KIV2 block."""
    r = row()
    out = {"status": "not assessed", "kiv2_copies": ".", "ci95_low": ".", "ci95_high": ".", "class": ".",
           "kiv2_depth": ".", "flank_depth": ".", "mapq0_share": ".", "message": "", "advice": ADVICE}
    reason = b["reason"]
    if reason == "no_data":
        out["message"] = "not assessed: the BAM has no chr6, so LPA has no reads"
    else:
        kdepth, fdepth, share = b["target_depth"], b["control_depth"], b["mapq0_share"]
        out.update(kiv2_depth=f"{kdepth:.2f}", flank_depth=f"{fdepth:.2f}", mapq0_share=f"{share:.3f}")
        if reason == "low_control_depth":
            out["status"] = "too little depth"
            out["message"] = (f"too little depth: {fdepth:.1f}x over the single-copy LPA exons, under the "
                              f"{r['_min_control']:g}x the count needs")
        elif reason == "low_mapq0_share":
            out["message"] = (f"not assessed: reads of mapping quality 0 make up {100 * share:.1f}% of the depth "
                              "over the KIV-2 exons, so the aligner left out the reads that fit several KIV-2 "
                              "copies (as novoalign does by default) and the count would come out low; realign "
                              "with minimap2 (step 02) or BWA-MEM")
        elif reason == "below_floor":
            out["message"] = (f"not assessed: the KIV-2 exons have {kdepth:.1f}x against {fdepth:.1f}x over the "
                              "single-copy exons, too little for any repeat count (a deletion or a mapping problem)")
        else:
            est, (lo, hi) = b["estimate"], b["ci95"]
            out.update(status="estimated", kiv2_copies=f"{est:.1f}", ci95_low=f"{lo:.1f}", ci95_high=f"{hi:.1f}",
                       **{"class": b["class"]},
                       message=(f"about {est:.0f} KIV-2 repeats in total (95% interval {lo:.1f} to {hi:.1f}): "
                                f"{b['reading']}"))
    line = (f"LPA KIV-2: {out['message']} (KIV-2 exons {out['kiv2_depth']}x, single-copy exons "
            f"{out['flank_depth']}x, MAPQ 0 share {out['mapq0_share']})")
    return [(k, out[k]) for k in ORDER], line


def assess(m, q1):
    """The result dict for the mean depths of read_bedcov: m over every read,
    q1 over reads of MAPQ 1 or more."""
    out = dict(table(cn.assess(row(), m, q1))[0])
    del out["advice"]
    return out


ORDER = ("status", "kiv2_copies", "ci95_low", "ci95_high", "class", "kiv2_depth", "flank_depth",
         "mapq0_share", "message", "advice")


NULL_LIMIT = ("A short-read fraction cannot say which chromosome or which repeat carries it, so it cannot say "
              "whether it sits on the small or the large isoform.")


def null_allele(counts, kiv2_copies=".", why_not=None):
    """[(metric, value)] of the KIV-2 4925G>A table from
    traits_report.read_marker_pileup's counts (None: unreadable), with the
    count N of step 46's table ('.' when there is none). why_not: the reason
    the reads cannot be trusted (the aligner left out the reads of MAPQ 0),
    which makes it 'not assessed' whatever the counts."""
    import traits_report as tr  # bin/: the copy-aware marker catalog
    m = tr.MARKERS["LPA_KIV2_4925"]
    a = tr.assess_marker("LPA_KIV2_4925", None if why_not else counts)
    st = a["status"]
    try:
        total = float(kiv2_copies) + 2
    except ValueError:
        total = None
    copies = "."
    if st in ("present", "uncertain") and total is not None:
        copies = f"{float(a['fraction']) * total:.1f}"
    name = f"{m['name']} (a splice variant that lowers the Lp(a) made from its chromosome, Coassin 2017)"
    if st == "not assessed":
        msg = f"{name}: not assessed, " + (why_not or "the reads at its seven GRCh38 copies could not be read.")
    else:
        msg = (f"{name}: {st}, {a['alt_reads']} of {a['depth']} reads at its seven GRCh38 copies carry it"
               + (f" (about {copies} of the {total:.0f} copies counted)" if copies != "." else "")
               + (f", under the {m['min_depth']} reads a reading needs" if st == "too little depth" else "")
               + ".")
        if st == "present":
            msg += (" Carriers had a lower Lp(a) than non-carriers with the same number of repeats, most of "
                    "all with few repeats. " + NULL_LIMIT)
        elif st == "not seen":
            msg += (" Not seen in the reads is not proof of absence: one carrying repeat among many gives few "
                    "reads, and a short-read fraction never says which chromosome or repeat holds a variant.")
        else:
            msg += " " + NULL_LIMIT
    return [("marker", m["name"]), ("rsid", m["dbsnp"][0] + " (first KIV-2 copy)"), ("status", st),
            ("alt_reads", a["alt_reads"]), ("depth", a["depth"]), ("fraction", a["fraction"]), ("copies", copies),
            ("kiv2_copies", kiv2_copies), ("message", msg)]


def write_null_allele(table_path, pileup_path, out_path):
    """Write the KIV-2 4925G>A table beside step 46's table (read for the
    count N; a table that is missing or unreadable gives no N). Returns the message."""
    import traits_report as tr
    t = {}
    try:
        with open(table_path, errors="replace") as f:
            t = dict(ln.rstrip("\n").split("\t", 1) for ln in f if "\t" in ln)
    except OSError:
        pass
    n = t.get("kiv2_copies", ".")
    why_not = None
    try:
        share = float(t.get("mapq0_share", "."))
    except ValueError:
        share = None
    if share is not None and share < row()["_mapq0_min"]:
        # The reads that fit several copies are the ones a carrying repeat
        # gives; without them the share at the seven copies is biased.
        why_not = (f"reads of mapping quality 0 make up {100 * share:.1f}% of the KIV-2 depth, so the aligner "
                   "left out the reads that fit several copies, the ones a carrying repeat gives.")
    try:
        counts = tr.read_marker_pileup(pileup_path, "LPA_KIV2_4925")
    except (OSError, ValueError) as e:
        print(f"WARNING: {e}; KIV-2 4925G>A is not assessed", file=sys.stderr)
        counts = None
    rows = null_allele(counts, n, why_not)
    tr.write_metric_table(out_path, rows)
    return dict(rows)["message"]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lpa_kiv2.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the regions")
    c = sub.add_parser("estimate", help="count the KIV-2 repeats from samtools bedcov over the regions")
    c.add_argument("--bedcov", required=True)
    c.add_argument("--bedcov-mapq1", required=True, help="samtools bedcov -j -Q 1 over the same regions")
    c.add_argument("--out", required=True)
    n = sub.add_parser("null-allele", help="the KIV-2 4925G>A table beside the count")
    n.add_argument("--pileup", required=True, help="samtools mpileup at traits_report.py --marker-regions LPA_KIV2_4925")
    n.add_argument("--table", required=True, help="the table estimate wrote (for the count)")
    n.add_argument("--out", required=True, help="the 4925G>A table")
    args = ap.parse_args(argv)
    if args.cmd == "null-allele":
        try:
            print("LPA " + write_null_allele(args.table, args.pileup, args.out))
        except OSError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        return 0
    if args.cmd == "bed":
        print("\n".join(bed_lines()))
        return 0
    return cn.main(["call", "--locus", LOCUS, "--bedcov", args.bedcov, "--bedcov-mapq1", args.bedcov_mapq1,
                    "--table", args.out])


if __name__ == "__main__":
    sys.exit(main())
