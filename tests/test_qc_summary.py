import csv
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest


spec = importlib.util.spec_from_file_location("qc_summary", Path(__file__).resolve().parents[1] / "workflow/scripts/qc_summary.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class QcSummaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.assembly = self.root / "assembly.fasta"
        self.assembly.write_text(">nuclear description\nAAAA\n>mito\nGCNN\n")

    def record(self, stage="raw"):
        busco = self.root / stage / "busco"
        merqury = self.root / stage / "merqury"
        busco.mkdir(parents=True)
        merqury.mkdir()
        (busco / "short_summary.specific.test.txt").write_text(
            "# The lineage dataset is: test_odb12 (number of BUSCOs: 100)\n"
            f"# Summarized benchmarking in BUSCO notation for file {self.assembly}\n"
            "C:99.0%[S:98.0%,D:1.0%],F:0.5%,M:0.5%,n:100,E:0.0%\n")
        label = f"S.{stage}.hifiasm"
        (merqury / "merqury.qv").write_text(f"{label}\t0\t3\t+inf\t0\n")
        (merqury / "completeness.stats").write_text(f"{label}\tall\t3\t4\t75\n")
        (merqury / "inputs.json").write_text(json.dumps({"assembly": str(self.assembly),
            "assembly_sha256": qc.digest(self.assembly), "label": label, "k": 3}))
        mapping = self.root / "flagstat.txt"
        mapping.write_text("500 + 0 in total (QC-passed reads + QC-failed reads)\n"
                           "100 + 0 paired in sequencing\n80 + 0 properly paired (80.00% : N/A)\n")
        return {"sample": "S", "assembler": "hifiasm", "stage": stage, "assembly": str(self.assembly),
                "busco": str(busco), "merqury": str(merqury), "mapping": str(mapping)}

    def test_fasta_metrics_include_short_organellar_records_and_explicit_gc_denominator(self):
        stats = qc.fasta_statistics(self.assembly)
        self.assertEqual((stats["sequence_count"], stats["total_bp"], stats["sequence_n50_bp"]), (2, 8, 4))
        self.assertEqual(stats["gap_n_bp"], 2)
        self.assertAlmostEqual(stats["gc_percent_acgt"], 100 / 3)
        self.assembly.write_text(">same\nAAAA\n>same\nCCCC\n")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            qc.fasta_statistics(self.assembly)

    def test_native_label_order_and_infinite_qv_are_preserved(self):
        record = self.record()
        rows = qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")
        self.assertTrue(math.isinf(rows[0]["merqury_qv"]))
        with (self.root / "table.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["merqury_qv"], "+inf")
        self.assertEqual(row["mapping_status"], "NOT_APPLICABLE_RAW")
        self.assertEqual(row["properly_paired_reads_percent"], "NA")

    def test_proper_pairing_uses_paired_reads_not_all_alignments(self):
        record = self.record("final")
        rows = qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")
        self.assertEqual(rows[0]["properly_paired_reads_percent"], 80)
        self.assertEqual(rows[0]["tidk_status"], "NOT_REQUESTED")
        Path(record["mapping"]).write_text("0 + 0 paired in sequencing\n0 + 0 properly paired\n")
        self.assertIsNone(qc.mapping_summary(record["mapping"])["properly_paired_reads_percent"])

    def test_selected_busco_specific_report_has_matching_reference_identity(self):
        record = self.record()
        (Path(record["busco"]) / "short_summary.generic.other.txt").write_text("not the selected report")
        report = qc.busco_summary(record["busco"], self.assembly)
        self.assertEqual(report["busco_lineage"], "test_odb12")
        with self.assertRaisesRegex(ValueError, "different FASTA"):
            qc.busco_summary(record["busco"], self.root / "different.fasta")

    def test_wrong_merqury_assembly_and_duplicate_row_identity_fail(self):
        record = self.record()
        with self.assertRaisesRegex(ValueError, "Duplicate QC identity"):
            qc.collect([record, record], self.root / "table.tsv", self.root / "provenance.json")
        source = Path(record["merqury"]) / "inputs.json"
        metadata = json.loads(source.read_text())
        metadata["assembly"] = str(self.root / "different.fasta")
        source.write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, "provenance"):
            qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")

    def test_merqury_rejects_changed_assembly_at_same_path_or_missing_hash(self):
        record = self.record()
        metadata = Path(record["merqury"]) / "inputs.json"
        old_content = self.assembly.read_text()
        self.assembly.write_text(old_content.replace("AAAA", "AAAT"))
        with self.assertRaisesRegex(ValueError, "provenance"):
            qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")
        self.assembly.write_text(old_content)
        receipt = json.loads(metadata.read_text())
        del receipt["assembly_sha256"]
        metadata.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "provenance"):
            qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")

    def test_optional_untestable_heterozygosity_is_not_imputed_as_zero(self):
        record = self.record("final")
        summary = self.root / "heterozygosity.json"
        summary.write_text(json.dumps({"sample": "S", "assembler": "hifiasm", "status": "UNTESTABLE_ZERO_CALLABLE", "genotype_model": "diploid",
                                       "callable_bases": 0, "nuclear_acgt_bases": 6,
                                       "heterozygous_snps": 0, "heterozygosity_percent": None}))
        record["heterozygosity"] = str(summary)
        receipt = {"state": "completed", "inputs": {"reference": {"path": str(self.assembly), "sha256": qc.digest(self.assembly)}},
                   "outputs": {summary.name: {"sha256": qc.digest(summary)}}}
        (summary.parent / "receipt.json").write_text(json.dumps(receipt))
        rows = qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")
        self.assertEqual(rows[0]["callable_acgt_percent"], 0)
        self.assertIsNone(rows[0]["heterozygosity_percent"])
        self.assertEqual(rows[0]["heterozygosity_status"], "UNTESTABLE_ZERO_CALLABLE")
        wrong = json.loads(summary.read_text())
        wrong["sample"] = "OTHER"
        summary.write_text(json.dumps(wrong))
        receipt["outputs"][summary.name]["sha256"] = qc.digest(summary)
        (summary.parent / "receipt.json").write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "provenance"):
            qc.collect([record], self.root / "wrong.tsv", self.root / "wrong.json")

    def test_tidk_summary_requires_matching_assembly_receipts(self):
        record = self.record("final")
        folder = self.root / "tidk"
        folder.mkdir()
        summary = folder / "assembly_summary.tsv"
        summary.write_text("screen_status\tbest_motif\tinput_sequences\tdetected_ends\tpossible_sequence_ends\n"
                           "CANDIDATE_IDENTIFIED\tAACCCT\t2\t2\t4\n")
        config = folder / "run_config.json"
        config.write_text(json.dumps({"input_path": str(self.assembly), "input_sha256": qc.digest(self.assembly)}))
        receipt_path = folder / "complete.json"
        receipt = {"output_sha256": {summary.name: qc.digest(summary), config.name: qc.digest(config)}}
        receipt_path.write_text(json.dumps(receipt))
        record["tidk"] = str(summary)
        rows = qc.collect([record], self.root / "table.tsv", self.root / "provenance.json")
        self.assertEqual(rows[0]["tidk_detected_ends"], 2)
        metadata = json.loads(config.read_text())
        metadata["input_path"] = str(self.root / "other.fasta")
        config.write_text(json.dumps(metadata))
        receipt["output_sha256"][config.name] = qc.digest(config)
        receipt_path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "provenance"):
            qc.collect([record], self.root / "wrong.tsv", self.root / "wrong.json")


if __name__ == "__main__":
    unittest.main()
