#!/usr/bin/env python3
"""Step 42 CLI goldens, rule-allele controls and og cn marker outputs.

Internal validation, phase and bounded-reader assertions live in
og-steps::traits::tests. The parity harness checks the old Python at oracle.ref.
"""
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import traits_report as tr  # noqa: E402
TABLE = os.path.join(REPO, "assets", "traits.tsv")
RULES = os.path.join(REPO, "assets", "genotype_rules.tsv")
CHECK = os.path.join(REPO, "scripts", "ci", "check-rule-alleles.sh")
OG = os.environ.get("OG", "")
SMN_HEADER = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n"
SMN1_NOT_PROOF = "This marker does not prove 2+0, and its absence does not rule 2+0 out."
FAILS = 0
APOE_GOLDEN = {
    ("T/T", "C/C"): "e3d7512b649f1d23f5acbcaa08da19bd6a009d5b39ac57cb6856591c50e7cfda",
    ("T/T", "C/T"): "62e910d8ec7ac45ab15e9c9ed85c93b128ac41703d1149a332082dc0d7b5f389",
    ("T/T", "T/T"): "41b4ad51c67028bb19af1b7ab7e4f414f5c0354e4c4a7c3c19309186b429f960",
    ("T/C", "C/C"): "08c5ed9e2898f2ae29cbd49d5ac97394758fa852f9d35e3c7173b2aa8d3111ad",
    ("T/C", "C/T"): "555c9cc0e8b59b10bdc5bdc2e33d93fb4ccbee1e6ea12916b17e73ffd6fdb8d8",
    ("T/C", "T/T"): "043f51b41486f217cac3a71e705e0e29c1a1ced9745e60abde2230df2fe901db",
    ("C/C", "C/C"): "ca7f36540cf1269b84649b4ad6ace19e0abc03eab7984d495de584540bb11dcc",
    ("C/C", "C/T"): "db499d1bfebedf77ed34a9d0525030d34be6336e685758a8efd53125197e3d00",
    ("C/C", "T/T"): "16e5214f90c5f2c401e5849dffd76a08dc3fc5ba2deb324d7fb02c76a92113be",
    ("T/C", "."): "34624db16a92981c0882edb3d8ce8d291162aa78b7c5208eba2360141e6ccbd8",
    ("T/G", "C/C"): "1d0aee0b8d241561ebd7477d6f4a8be274a48157a297abecd444c6c6a0a6abac",
}
def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def dbsnp_checks(rules):
    """10. scripts/ci/check-rule-alleles.sh on the rules and on planted copies. Runs first,
    so a rules table that breaks the other checks still shows what dbSNP says."""
    work = tempfile.mkdtemp()
    try:
        def dbsnp(path, *extra):
            r = subprocess.run(["bash", CHECK, *extra, path], capture_output=True, text=True)
            return r.returncode, r.stdout + r.stderr

        rc, log = dbsnp(RULES)
        print("  " + log.strip().replace("\n", "\n  "))
        n = sum(len(v["sites"]) for v in rules.values())
        check(f"check-rule-alleles.sh: every rule site matches dbSNP ({n} sites)",
              rc == 0 and f"{n} rule sites checked" in log, log)
        text = open(RULES).read()
        site = "site\tTAS2R38\trs1726866\tchr7\t141972905\t"

        def planted(desc, new, want):
            f = os.path.join(work, "planted.tsv")
            write(f, text.replace(site + "G\tA", site + new))
            rc, log = dbsnp(f)
            check(f"check-rule-alleles.sh control: {desc} fails, naming it", rc == 1 and want in log, (rc, log))

        planted("ref and alt swapped", "A\tG", "rs1726866: ref and alt are swapped")
        planted("the minus-strand bases", "C\tT", "rs1726866: the alleles are the minus-strand bases")
        f = os.path.join(work, "planted.tsv")
        write(f, text.replace("\trs1726866\tchr7\t141972905\t", "\trs1726867\tchr7\t141972905\t"))
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: an rsID the dbSNP table lacks fails",
              rc == 1 and "rs1726867: not in the pinned dbSNP table" in log, (rc, log))
        write(f, text.replace("\trs1726866\tchr7\t141972905\t", "\trs1726866\tchr7\t141972906\t"))
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: a wrong position fails",
              rc == 1 and "dbSNP at chr7:141972905" in log, (rc, log))
        write(f, "# no sites\n")
        rc, log = dbsnp(f)
        check("check-rule-alleles.sh control: a rules file with no site fails", rc == 1 and "no site line" in log,
              (rc, log))
        root = os.path.join(work, "root")
        os.makedirs(os.path.join(root, "scripts", "ci"))
        shutil.copy(os.path.join(REPO, "versions.env"), root)
        copy = os.path.join(root, "scripts", "ci", "check-rule-alleles.sh")
        write(copy, open(CHECK).read().replace("\n#: # build ", "\n#: # build x"))
        r = subprocess.run(["bash", copy, RULES], capture_output=True, text=True)
        check("check-rule-alleles.sh control: a dbSNP table of another build than versions.env pins fails",
              r.returncode == 1 and "versions.env pins" in r.stderr, (r.returncode, r.stdout, r.stderr))
    finally:
        shutil.rmtree(work)


MULTI_SNP = {
    "APOL1": [(("G/G", "G/G"), "G1/G1"), (("A/G", "T/G"), "G0/G1")],
    "SERPINA1": [(("T/T", "T/T"), "PiZZ"), (("T/A", "C/T"), "PiSZ"), (("A/A", "C/C"), "PiSS")],
    "HFE": [(("A/A", "C/C"), "C282Y/C282Y"), (("G/A", "C/G"), "C282Y/H63D"), (("G/G", "G/G"), "H63D/H63D")],
    "GC": [(("C/C", "G/G"), "Gc1S/Gc1S"), (("A/C", "G/T"), "Gc1S/Gc2"), (("A/A", "T/T"), "Gc2/Gc2")],
    "CAD9p21": [(("G/G", "C/C"), "two risk copies"), (("A/G", "G/C"), "one risk copy")],
    "MAPT": [(("G/G", "G/G"), "H2/H2"), (("A/G", "A/G"), "H1/H2")],
}
# Independent labels for homozygous GRCh38 reference blocks.
MULTI_SNP_REF = {"APOL1": "G0/G0", "SERPINA1": "PiMM", "HFE": "wt/wt", "GC": "Gc1F/Gc1F",
                 "CAD9p21": "no risk copy", "MAPT": "H1/H1"}



def says_not_proof(message):
    """The SMN1 marker's message must say the marker does not prove 2+0."""
    return "does not prove 2+0" in message


def marker_errors(name, m=None, dbsnp=None):
    """Problems of one marker (empty when it is consistent): the dbSNP record
    is one of its copies with its ref and alt, on the plus strand and not
    swapped; every copy sits in the middle of the same context, whose middle
    base is the ref."""
    m = tr.MARKERS[name] if m is None else m
    rs, chrom, pos, ref, alt = m["dbsnp"] if dbsnp is None else dbsnp
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    errs = []
    if chrom != m["chrom"] or pos not in m["copies"]:
        errs.append(f"{name}: dbSNP places {rs} at {chrom}:{pos}, not at one of the marker's copies")
    if (m["ref"], m["alt"]) == (alt, ref):
        errs.append(f"{name}: ref and alt are swapped (marker {m['ref']}>{m['alt']}, dbSNP {ref}>{alt})")
    elif (m["ref"], m["alt"]) == (comp.get(ref), comp.get(alt)):
        errs.append(f"{name}: the alleles are the minus-strand bases (marker {m['ref']}>{m['alt']}, "
                    f"dbSNP {ref}>{alt}): write the plus-strand bases")
    elif (m["ref"], m["alt"]) != (ref, alt):
        errs.append(f"{name}: marker {m['ref']}>{m['alt']}, dbSNP {ref}>{alt}")
    ctx = m["context"]
    if len(ctx) % 2 != 1 or ctx[len(ctx) // 2] != m["ref"]:
        errs.append(f"{name}: the middle base of the context is not the ref {m['ref']}")
    if len(set(m["copies"])) != len(m["copies"]):
        errs.append(f"{name}: a copy is listed twice")
    return errs


def marker_checks():
    """11. The copy-aware markers."""
    for name, m in tr.MARKERS.items():
        check(f"marker {name}: the pinned dbSNP record agrees with the marker", marker_errors(name) == [],
              marker_errors(name))
        rs, chrom, pos, ref, alt = m["dbsnp"]
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        e = marker_errors(name, dict(m, ref=alt, alt=ref))
        check(f"marker {name} control: ref and alt swapped is caught", any("swapped" in x for x in e), e)
        e = marker_errors(name, dict(m, ref=comp[ref], alt=comp[alt],
                                        context=m["context"][:10] + comp[ref] + m["context"][11:]))
        check(f"marker {name} control: the minus-strand bases are caught", any("minus-strand" in x for x in e), e)
        e = marker_errors(name, dbsnp=(rs, chrom, pos + 1, ref, alt))
        check(f"marker {name} control: a dbSNP position off the copies is caught",
              any("not at one of the marker's copies" in x for x in e), e)
        e = marker_errors(name, dict(m, context=m["context"][:10] + alt + m["context"][11:]))
        check(f"marker {name} control: a context whose middle base is not the ref is caught",
              any("middle base" in x for x in e), e)
    check("LPA_KIV2_4925: seven copies (six KIV-2 units and KIV-1)", len(tr.MARKERS["LPA_KIV2_4925"]["copies"]) == 7)
    if not OG:
        print("skip: OG is not set (OG=rust/target/release/og): og cn marker-regions and marker-read did not run")
        return
    for name, m in tr.MARKERS.items():
        res = og_cn("marker-regions", name)
        check(f"marker {name}: og cn marker-regions prints one region per copy of MARKERS",
              res.returncode == 0 and res.stdout.split() == [f"{m['chrom']}:{p}-{p}" for p in m["copies"]], res)

    sm = tr.MARKERS["SMN1_27134"]
    work = tempfile.mkdtemp()
    try:
        def pileup(t_reads, g_reads, extra=""):
            p = os.path.join(work, "p.tsv")
            write(p, "".join(f"chr5\t{pos}\tN\t{t_reads + g_reads}\t^I{'T' * t_reads}{'g' * g_reads}$\t.\n"
                             for pos in sm["copies"]) + extra)
            return p

        def smn(pscn="2,2", pqual="30,30", pfilt="PASS", agcn="4", homolog="chr5:70049000-70078000"):
            p = os.path.join(work, "smn.tsv")
            write(p, SMN_HEADER + f"chr5\t70925000\t70953000\tSMN1\tPASS\t{agcn}\t50\t{pfilt}\t{pscn}\t{pqual}\t{homolog}\n")
            return p

        def table(p, smn_path):
            out = os.path.join(work, "o.tsv")
            if os.path.exists(out):
                os.remove(out)
            res = og_cn("marker-read", "SMN1_27134", "--pileup", p, "--copy-table", smn_path, "--out", out)
            return res, dict(x.split("\t", 1) for x in open(out).read().splitlines()[1:]) if os.path.exists(out) else {}

        res, t = table(pileup(15, 5), smn())
        check("SMN1 marker: a quarter of the reads with G is present, exit 0",
              res.returncode == 0 and t.get("status") == "present" and t.get("alt_reads") == "10"
              and t.get("depth") == "40" and t.get("fraction") == "0.2500", (res, t))
        check("SMN1 marker: copies = fraction x agCN (0.25 x 4)", t.get("copies") == "1.0" and t.get("agcn") == "4", t)
        check("SMN1 marker: the SMN1 psCN at exon 7 is read", t.get("smn1_pscn") == "2 (quality 30, filter PASS)", t)
        check("SMN1 marker: present with two copies says it raises the chance of 2+0 and does not prove it",
              "higher chance" in t.get("message", "") and says_not_proof(t.get("message", "")), t)
        res, t = table(pileup(20, 0), smn())
        check("SMN1 marker: reads with T only (the reference) is not seen, copies '.'",
              t.get("status") == "not seen" and t.get("copies") == "." and "lower" in t.get("message", ""), t)
        res, t = table(pileup(20, 0, extra="chr5\t70952074\tN\t2\tGG\t.\n"), smn())
        check("SMN1 marker: 2 reads with G is uncertain", t.get("status") == "uncertain", t)
        res, t = table(pileup(3, 1), smn())
        check("SMN1 marker: under 10 reads is too little depth", t.get("status") == "too little depth", t)
        res, t = table(pileup(15, 5, extra="chr5\t70952075\tN\t30\t" + "G" * 30 + "\t.\n"), smn())
        check("SMN1 marker: reads at another position are not counted", t.get("alt_reads") == "10", t)
        bad = os.path.join(work, "bad.tsv")
        write(bad, "nonsense\n")
        res, t = table(bad, smn())
        check("SMN1 marker: an unreadable pileup is not assessed, exit 0",
              res.returncode == 0 and t.get("status") == "not assessed", (res, t))
        res, t = table(pileup(15, 5), smn(pscn="1,3"))
        check("SMN1 marker: another SMN1 copy number says the marker refines a two-copy result only",
              "SMN1 has 1 copy here" in t.get("message", ""), t)
        res, t = table(pileup(15, 5), smn(pqual="12,30"))
        check("SMN1 marker: a psCN of quality under 20 is not resolved",
              t.get("smn1_pscn") == "2 (quality 12, filter PASS), not resolved"
              and "SMN1 copy number: 2 (quality 12, filter PASS), not resolved;" in t.get("message", ""), t)
        # the SMN table with SMN2's region first: exon 7 of SMN1 is the homologous region
        p = os.path.join(work, "smn2.tsv")
        write(p, SMN_HEADER + "chr5\t70049000\t70078000\tSMN1\tPASS\t4\t50\tPASS\t2,2\t30,40\tchr5:70925000-70953000\n")
        res, t = table(pileup(15, 5), p)
        check("SMN1 marker: exon 7 found through the homologous region (its own psCN and quality)",
              t.get("smn1_pscn") == "2 (quality 40, filter PASS)" and t.get("agcn") == "4", t)
        res, t = table(pileup(15, 5), os.path.join(work, "missing.tsv"))
        check("SMN1 marker: no SMN table says so", "no SMN table" in t.get("smn1_pscn", ""), t)
        res, t = table(pileup(15, 5), work)
        check("SMN1 marker: an SMN table that cannot be read keeps the marker's reading, exit 0",
              res.returncode == 0 and t.get("status") == "present" and t.get("agcn") == "."
              and "could not be read" in t.get("smn1_pscn", ""), (res, t))

        # Every status with every psCN: the message says the marker does not prove 2+0.
        msgs = []
        for counts in [(30, 10), (40, 0), (38, 2), (3, 1), None]:
            if counts is None:
                p = bad
            else:   # every read at the SMN1 copy
                p = os.path.join(work, "all.tsv")
                write(p, f"chr5\t{sm['copies'][0]}\tN\t{sum(counts)}\t{'T' * counts[0]}{'G' * counts[1]}\t.\n")
            for path in ("2,2", "1,3", "5,5", None):
                smn_path = (os.path.join(work, "missing.tsv") if path is None
                            else smn(pqual=path) if path == "5,5" else smn(pscn=path))
                msgs.append(table(p, smn_path)[1].get("message", ""))
        check(f"SMN1 marker: every message ({len(msgs)}) says it does not prove 2+0",
              len(msgs) == 20 and all(says_not_proof(x) for x in msgs), [x for x in msgs if not says_not_proof(x)])
        check("SMN1 marker control: a message without the sentence fails that test",
              not says_not_proof(msgs[0].replace(SMN1_NOT_PROOF, "")))
    finally:
        shutil.rmtree(work)



def traits_tables(genotypes, out, rules_out, sample="S"):
    """Use og in the Rust job. Unit shards without og use the pinned oracle."""
    args = ["--table", TABLE, "--rules", RULES, "--genotypes", genotypes,
            "--sample", sample, "--out", out, "--rules-out", rules_out]
    if OG:
        subprocess.run([OG, "traits", "table", *args], check=True, capture_output=True, text=True)
    else:
        with tempfile.TemporaryDirectory(prefix="traits-oracle-") as work:
            ref = "669a8196cf4c377200f9b4af9c235ef85fb774a1"
            if subprocess.run(["git", "-C", REPO, "cat-file", "-e", ref + "^{commit}"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode:
                subprocess.run(["git", "-C", REPO, "fetch", "--depth=1", "origin", ref], check=True)
            source = subprocess.check_output(["git", "-C", REPO, "show", ref + ":bin/traits_report.py"])
            path = os.path.join(work, "traits_report.py")
            with open(path, "wb") as f:
                f.write(source)
            subprocess.run([sys.executable, path, *args], check=True, capture_output=True, text=True)


def og_cn(*args):
    return subprocess.run([OG, "cn", *args], capture_output=True, text=True)



def plain_tsv_checks():
    """Reject CSV quoting before it can produce a table the collector cannot read."""
    import collect_summary as cs
    with tempfile.TemporaryDirectory() as work:
        folder = os.path.join(work, "traits")
        os.makedirs(folder)
        g = os.path.join(work, "g.tsv")
        write(g, "")
        out, rout = os.path.join(folder, "S_traits.tsv"), os.path.join(folder, "S_genotype_rules.tsv")
        source = open(TABLE).read()
        rules_source = open(RULES).read()
        for name, text, rules_text in [
            ("quoted header", source.replace("chrom\t", '\"chrom\"\t', 1), rules_source),
            ("quoted catalog field", source.replace("\tMCM6\t", '\t\"MCM6\"\t', 1), rules_source),
            ("quoted rule text", source, rules_source.replace("\tAPOE\tAPOE\t", '\tAPOE\t\"APOE\"\t', 1)),
            ("embedded tab", source.replace("\tMCM6\t", '\t\"MCM6\textra\"\t', 1), rules_source),
        ]:
            # Keep complete, collector-readable outputs from the previous run.
            traits_tables(g, out, rout)
            previous = open(out, "rb").read(), open(rout, "rb").read()
            bad, badrules = os.path.join(work, "bad.tsv"), os.path.join(work, "badrules.tsv")
            write(bad, text)
            write(badrules, rules_text)
            res = subprocess.run([OG, "traits", "table", "--table", bad, "--rules", badrules,
                                  "--genotypes", g, "--sample", "S", "--out", out,
                                  "--rules-out", rout], capture_output=True, text=True)
            check(f"plain TSV control: {name} is refused", res.returncode == 1, res.stderr)
            check(f"plain TSV control: {name} preserves both outputs",
                  previous == (open(out, "rb").read(), open(rout, "rb").read()))
            try:
                _, collected = cs.sec_traits(work, "S")
                readable = bool(collected["rows"])
            except cs.Unreadable:
                readable = False
            check(f"plain TSV control: {name} leaves collector-readable output", readable)


def main():
    rules = {}
    with open(RULES) as f:
        for cells in csv.reader((line for line in f if not line.startswith("#")), delimiter="\t"):
            if cells and cells[0] == "site":
                rules.setdefault(cells[1], {"sites": []})["sites"].append(cells[2:])
    dbsnp_checks(rules)
    marker_checks()
    if not OG:
        check("OG_REQUIRED needs the og binary", not os.environ.get("OG_REQUIRED"))
        print("skip: OG is not set: traits CLI and cn marker checks did not run")
        return int(bool(FAILS))
    plain_tsv_checks()
    with tempfile.TemporaryDirectory() as work:
        def table(records):
            g = os.path.join(work, "g.tsv")
            write(g, "".join("\t".join(rec) + "\n" for rec in records))
            os.makedirs(os.path.join(work, "traits"), exist_ok=True)
            out, rout = os.path.join(work, "traits", "S_traits.tsv"), os.path.join(work, "traits", "S_genotype_rules.tsv")
            res = subprocess.run([OG, "traits", "table", "--table", TABLE, "--rules", RULES,
                                  "--genotypes", g, "--sample", "S", "--out", out,
                                  "--rules-out", rout], capture_output=True, text=True)
            check("traits table exits 0", res.returncode == 0, res.stderr)
            if res.returncode:
                return [], []
            with open(out) as f:
                traits = list(csv.DictReader(f, delimiter="\t"))
            with open(rout) as f:
                lines = list(csv.DictReader(f, delimiter="\t"))
            return traits, lines
        for (a, b), want in APOE_GOLDEN.items():
            records = [("chr19", pos, ref, alt, gt) for pos, ref, alt, gt in
                       [("44908684", "T", "C", a), ("44908822", "C", "T", b)] if gt != "."]
            traits, lines = table(records)
            apoe = next(r for r in traits if "+" in r["rsid"])
            got = hashlib.sha256("\t".join(apoe[k] for k in tr.OUT_COLUMNS).encode()).hexdigest()
            check(f"APOE golden {a}, {b}", got == want, got)
        for rid, cases in MULTI_SNP.items():
            for gts, label in cases:
                records = [(s[1], s[2], s[3], s[4], gt) for s, gt in zip(rules[rid]["sites"], gts)]
                _, lines = table(records)
                got = next(r for r in lines if r["rule"] == rid)
                check(f"{rid} {gts}: {label}", got["status"] == "called" and got["result"] == label, got)
        import collect_summary as cs
        import render_report as rr
        check("every multi-SNP rule is in the catalog", set(MULTI_SNP) <= set(rules))
        for rid, cases in MULTI_SNP.items():
            sites = rules[rid]["sites"]
            records = [(s[1], s[2], s[3], ".", f"{s[3]}/{s[3]}", "40") for s in sites]
            _, lines = table(records)
            got = next(r for r in lines if r["rule"] == rid)
            check(f"{rid}, reference blocks only: {MULTI_SNP_REF[rid]}",
                  (got["status"], got["result"]) == ("called", MULTI_SNP_REF[rid]), got)
            gts, want = cases[0]
            a, b = gts[0].split("/")
            s0 = sites[0]
            flip = s0[3] if b == s0[4] else s0[4]
            flipped = (f"{a}/{flip}",) + tuple(gts[1:])
            records = [(s[1], s[2], s[3], s[4], gt) for s, gt in zip(sites, flipped)]
            _, lines = table(records)
            got = next(r for r in lines if r["rule"] == rid)
            check(f"control: {rid} with one allele flipped is not {want}", got["result"] != want, got)
            # Collect and render the actual CLI output, independently for each rule.
            records = [(s[1], s[2], s[3], s[4], gt) for s, gt in zip(sites, gts)]
            _, lines = table(records)
            read = cs.genotype_rules(work, "S") or []
            page = "\n".join(rr.rules_html({"rules": read}, "S"))
            got = next(r for r in lines if r["rule"] == rid)
            check(f"HTML report: {rid} with result {want}", rr.E(got["trait"]) in page and rr.E(want) in page, page[:300])
        records = [(s[1], s[2], s[3], s[4], gt) for s, gt in zip(rules["APOL1"]["sites"], ("A/G", "T/G"))]
        _, lines = table(records)
        reading = next(r for r in lines if r["rule"] == "APOL1")["reading"]
        check("APOL1 G0/G1 says G2 is not read and the other copy could carry it",
              "G2 is not read" in reading and "G1/G2" in reading, reading)
    print("RESULT:", "ALL PASS" if not FAILS else f"{FAILS} FAILED")
    return int(bool(FAILS))


if __name__ == "__main__":
    sys.exit(main())
