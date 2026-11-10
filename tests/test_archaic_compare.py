#!/usr/bin/env python3
"""bin/archaic_compare.py: merging, Chen's file format, and the comparison table."""
import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "archaic_compare.py")
sys.path.insert(0, os.path.join(HERE, "..", "bin"))
import archaic_compare as ac  # noqa: E402

CHEN = (
    "chr start   end LOD maxLOD  size    pop anc ID   \n"
    "22\t100\t200\t9\t9\t100\tCEU\tEUR\tA\n"
    "22\t150\t300\t9\t9\t150\tCEU\tEUR\tA\n"   # overlaps the first: union 200 bp
    "22\t1000\t1100\t9\t9\t100\tCEU\tEUR\tA\n"
    "21\t10\t60\t9\t9\t50\tCEU\tEUR\tA\n"
    "22\t100\t150\t9\t9\t50\tGBR\tEUR\tB\n"
    "22\t100\t900\t9\t9\t800\tYRI\tAFR\tC\n"
)
MINE = (
    "chrom\tstart\tend\tlength_bp\tancestry\tmatched\tlod\tgenes\n"
    "chr22\t100\t400\t300\tNeanderthal\tAltai\t5\t\n"
    "chr22\t5000\t6000\t1000\tDenisovan\tDenisova\t5\t\n"
)


class Compare(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.chen = os.path.join(self.d, "chen.txt")
        self.mine = os.path.join(self.d, "mine.tsv")
        open(self.chen, "w").write(CHEN)
        open(self.mine, "w").write(MINE)

    def run_cmp(self, *extra):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ac.main(["compare", "--chen-file", self.chen, "--sample", "A", "--group", "EUR",
                     "--chroms", "21 22", "--mine", f"x={self.mine}", *extra])
        return buf.getvalue()

    def test_merged_bp(self):
        self.assertEqual(ac.merged_bp([(100, 200), (150, 300), (1000, 1100)]), 300)
        self.assertEqual(ac.merged_bp([]), 0)

    def test_person_total_is_the_union(self):
        out = self.run_cmp()
        # chr22: 100-300 and 1000-1100 = 300; chr21: 50
        self.assertIn("| 22 | 300 |", out)
        self.assertIn("| 21 | 50 |", out)
        self.assertIn("| all | 350 |", out)

    def test_only_neanderthal_rows_of_mine_count(self):
        out = self.run_cmp()
        self.assertIn("x: 300 bp, 0.86 times the released A figure;", out)

    def test_overlap_with_the_released_calls(self):
        out = self.run_cmp()
        # mine chr22 100-400; released chr22 100-300 and 1000-1100: 200 bp shared
        self.assertIn("200 bp under the released calls, 100 bp only here, 150 bp only released", out)

    def test_intersect_bp(self):
        self.assertEqual(ac.intersect_bp([(0, 10), (20, 30)], [(5, 25)]), 10)
        self.assertEqual(ac.intersect_bp([(0, 10)], []), 0)

    def test_all_row_maximum_is_a_person_total(self):
        # A is largest on chr21 (100), B on chr22 (100): summed maxima would be 200, no person has it
        open(self.chen, "w").write(
            "chr start end LOD maxLOD size pop anc ID\n"
            "21\t0\t100\t9\t9\t100\tCEU\tEUR\tA\n22\t0\t10\t9\t9\t10\tCEU\tEUR\tA\n"
            "21\t0\t10\t9\t9\t10\tGBR\tEUR\tB\n22\t0\t100\t9\t9\t100\tGBR\tEUR\tB\n")
        out = self.run_cmp()
        self.assertIn("| all | 110 | 110 | 110 |", out)

    def test_group_spread_uses_the_group_only(self):
        out = self.run_cmp()
        self.assertIn("2 EUR people", out)  # A and B, not the AFR person C
        self.assertIn("highest 350", out)

    def test_unknown_sample_stops(self):
        with self.assertRaises(SystemExit):
            ac.main(["compare", "--chen-file", self.chen, "--sample", "Z", "--chroms", "22",
                     "--mine", f"x={self.mine}"])

    def test_empty_group_stops(self):
        with self.assertRaises(SystemExit) as r:
            ac.main(
                ["compare", "--chen-file", self.chen, "--sample", "A", "--group", "NOPE", "--chroms", "22",
                 "--mine", f"x={self.mine}"])
        self.assertIn("NOPE", str(r.exception))

    def test_wrong_segments_file_stops(self):
        open(self.mine, "w").write("a\tb\n")
        with self.assertRaises(SystemExit):
            self.run_cmp()

    def test_genomes_union_row(self):
        hdr = "ID\tchrom\tstart\tend\tslod\n"
        rows = {
            "Altai": "SAMPLE\t22\t100000\t300000\t5\n",
            "Vindija": "SAMPLE\t22\t1000000\t1100000\t5\n",
            "Denisova": "SAMPLE\t22\t5000000\t5100000\t5\n",
        }
        for n, body in rows.items():
            open(os.path.join(self.d, f"chr22.{n}.txt"), "w").write(hdr + body)
        spans = os.path.join(self.d, "spans.tsv")
        open(spans, "w").write("chrom\tstart\tend\n22\t1\t9999\n")
        chen = os.path.join(self.d, "chen_big.txt")
        open(chen, "w").write(
            "chr start end LOD maxLOD size pop anc ID\n"
            "22\t100000\t300000\t9\t9\t200000\tCEU\tEUR\tA\n"
            "22\t1000000\t1100000\t9\t9\t100000\tCEU\tEUR\tA\n")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ac.main(["genomes", "--chen-file", chen, "--sample", "A",
                     "--ibdmix-dir", self.d, "--spans", spans])
        out = buf.getvalue()
        # released A on chr22: 300,000 bp; Altai and Vindija together are exactly that
        self.assertIn("| Altai+Vindija | 300000 | 300000 | 0 | 0 |", out)
        self.assertIn("| Denisova | 100000 | 0 | 100000 | 300000 |", out)
        self.assertIn("Released A: 300000 bp on these 1 chromosomes.", out)

    def test_cli(self):
        r = subprocess.run([sys.executable, "-I", SCRIPT, "compare", "--chen-file", self.chen, "--sample", "A",
                            "--chroms", "22", "--mine", f"x={self.mine}"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == "__main__":
    unittest.main()
