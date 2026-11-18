import importlib.util
from pathlib import Path
import tempfile
import unittest


spec = importlib.util.spec_from_file_location(
    "integrate_mito", Path(__file__).resolve().parents[1] / "workflow/scripts/integrate_mito.py")
mito = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mito)


def paf(query, qlen, qs, qe, target, tlen, ts, te, matches=None):
    span = qe - qs
    return "\t".join(map(str, [query, qlen, qs, qe, "+", target, tlen, ts, te,
                                span if matches is None else matches, span, 60])) + "\n"


class MitoIntegrationTests(unittest.TestCase):
    def run_integration(self, assembly, mitochondria, alignments):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        (root / "assembly.fa").write_text(assembly)
        (root / "mito.fa").write_text(mitochondria)
        (root / "matches.paf").write_text(alignments)
        removed = mito.integrate(root / "assembly.fa", root / "mito.fa", root / "matches.paf",
                                 "SAMPLE", root / "final.fa", root / "removed.txt", root / "report.tsv")
        return removed, mito.read_fasta(root / "final.fa"), (root / "report.tsv").read_text()

    def test_circular_split_removes_target_and_preserves_nuclear_flanks(self):
        mitochondrial = "A" * 100
        assembly = ">circle\n" + mitochondrial + "\n>nuclear\n" + "C" * 100 + mitochondrial + "G" * 100 + "\n"
        alignments = (paf("oatk", 100, 0, 56, "circle", 100, 44, 100)
                      + paf("oatk", 100, 56, 100, "circle", 100, 0, 44)
                      + paf("oatk", 100, 0, 100, "nuclear", 300, 100, 200))
        removed, records, report = self.run_integration(assembly, ">oatk\n" + mitochondrial + "\n", alignments)
        self.assertEqual(removed, {"circle"})
        self.assertEqual([h for h, _ in records], ["nuclear", "SAMPLE_mito"])
        self.assertEqual(records[0][1], "C" * 100 + mitochondrial + "G" * 100)
        self.assertIn("circle\toatk\t1.00000000\t1.00000000\tremove", report)

    def test_overlapping_hits_do_not_inflate_query_coverage(self):
        alignments = (paf("oatk", 100, 0, 60, "target", 100, 0, 60)
                      + paf("oatk", 100, 0, 60, "target", 100, 40, 100))
        removed, records, _ = self.run_integration(
            ">target\n" + "A" * 100 + "\n", ">oatk\n" + "A" * 100 + "\n", alignments)
        self.assertEqual(removed, set())
        self.assertEqual(len(records), 2)

    def test_low_identity_does_not_remove_target(self):
        removed, _, _ = self.run_integration(
            ">target\n" + "A" * 100 + "\n", ">oatk\n" + "A" * 100 + "\n",
            paf("oatk", 100, 0, 100, "target", 100, 0, 100, matches=90))
        self.assertEqual(removed, set())

    def test_multiple_mitochondrial_records_receive_unique_names(self):
        _, records, _ = self.run_integration(
            ">nuclear\nAAAA\n", ">oatk1\nCCCC\n>oatk2\nGGGG\n", "")
        self.assertEqual([h for h, _ in records], ["nuclear", "SAMPLE_mito_1", "SAMPLE_mito_2"])

    def test_unknown_alignment_ids_fail(self):
        with self.assertRaisesRegex(ValueError, "PAF ID absent"):
            self.run_integration(
                ">target\nAAAA\n", ">oatk\nAAAA\n", paf("oatk", 4, 0, 4, "wrong", 4, 0, 4))


if __name__ == "__main__":
    unittest.main()
