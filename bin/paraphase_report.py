#!/usr/bin/env python3
# remnant: kept for collect_summary and couple_report; deleted in lane R4
"""paraphase_report.py: what is left of step 35b's two tables once they
moved to `og paraphase check-bam` and `og paraphase table`
(rust/crates/og-steps/src/paraphase.rs). Only the rule the summary and the
couple report still import is kept: SMN1 from step 35's Parascopy table, set
beside SMN1 from Paraphase, so both reports judge agreement the same way.
`og dev remnants` (tests/parity/remnants.tsv) checks that nothing else comes
back. Standard library only.
"""
import re

# SMN1 exon 7, c.840 (the base that tells SMN1 from SMN2), GRCh38 1-based.
SMN1_EXON7 = ("chr5", 70951946)


def parascopy_smn1(path):
    """SMN1 copy number from step 35's table (<s>_smn_copy_number.tsv): the
    paralog-specific value of the row, or of the homologous region, that holds
    SMN1 exon 7. {'cn': int or None, 'text': str}; cn is None unless the
    filter is PASS and the quality at least 20 (Parascopy's threshold)."""
    chrom, pos = SMN1_EXON7
    with open(path, errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            regions = []
            try:
                regions.append((r["chrom"], int(r["start"]) + 1, int(r["end"])))
            except (KeyError, ValueError):
                continue
            for h in (r.get("homologous_regions") or "*").split(","):
                m = re.match(r"^([^:]+):(\d+)-(\d+)", h.strip())
                if m:
                    regions.append((m.group(1), int(m.group(2)), int(m.group(3))))
            idx = next((i for i, (c, s, e) in enumerate(regions) if c == chrom and s <= pos <= e), None)
            if idx is None:
                continue
            cns = (r.get("psCN") or "").split(",")
            quals = (r.get("psCN_qual") or "").split(",")
            cn_txt = cns[idx] if idx < len(cns) else "?"
            q_txt = quals[idx] if idx < len(quals) else "*"
            filt = r.get("psCN_filter") or "."
            try:
                q = float(q_txt)
            except ValueError:
                q = None
            ok = filt == "PASS" and q is not None and q >= 20 and cn_txt.isdigit()
            return {"cn": int(cn_txt) if ok else None,
                    "text": f"{cn_txt} (quality {q_txt}, filter {filt})" + ("" if ok else ", not resolved")}
    return {"cn": None, "text": "no row covers SMN1 exon 7"}


def compare_smn1(paraphase_cn, parascopy):
    """'agree', 'disagree', or None when either tool has no SMN1 call."""
    if paraphase_cn is None or parascopy is None or parascopy.get("cn") is None:
        return None
    return "agree" if int(paraphase_cn) == parascopy["cn"] else "disagree"
