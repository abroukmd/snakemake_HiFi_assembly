"""Known-answer checks for nQuire completion, missing data, and model ranking."""

import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/run_nquire.py"
SPEC = importlib.util.spec_from_file_location("run_nquire", SCRIPT)
NQUIRE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NQUIRE)
NATIVE_NQUIRE = os.environ.get("NQUIRE_TEST_EXECUTABLE") or shutil.which("nQuire")
SAMTOOLS = os.environ.get("SAMTOOLS_TEST_EXECUTABLE") or shutil.which("samtools")


def bin_data(positions, extended=False):
    header = struct.pack(">IBI", positions, int(extended), 200)
    if extended:
        header += struct.pack(">I", 123456)
    records = []
    for index in range(positions):
        record = struct.pack(">II", 0, index) if extended else b""
        records.append(record + struct.pack(">HBHH", 20, 2, 8, 12))
    return header + b"".join(records)


class NquireTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="hifi-nquire-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def result(self, content=None):
        path = self.root / "model.results"
        path.write_text(content if content is not None else
                        "\t".join(NQUIRE.HEADER) + "\nS-denoised.bin\t100\t98\t90\t80\t2\t10\t20\n")
        return path

    def fake_executable(self, raw_positions=3, denoised_positions=2, fail_step=None,
                        tie=False, invalid_model=False):
        path = self.root / "nQuire"
        settings = dict(raw=raw_positions, denoised=denoised_positions, fail=fail_step,
                        tie=tie, invalid_model=invalid_model)
        path.write_text("#!" + sys.executable + "\n" +
                        "settings = " + repr(settings) + "\n" + '''
from pathlib import Path
import struct
import sys
step = sys.argv[1]
with Path("calls.txt").open("a") as out:
    out.write(step + "\\n")
if step == settings["fail"]:
    print("forced native failure", file=sys.stderr)
    raise SystemExit(7)
def write_bin(name, positions):
    Path(name).write_bytes(struct.pack(">IBI", positions, 0, 200) +
                           struct.pack(">HBHH", 20, 2, 8, 12) * positions)
if step == "create":
    write_bin(sys.argv[sys.argv.index("-o") + 1] + ".bin", settings["raw"])
elif step == "denoise":
    write_bin(sys.argv[sys.argv.index("-o") + 1] + ".bin", settings["denoised"])
elif step == "lrdmodel":
    print("file\\tfree\\tdip\\ttri\\ttet\\td_dip\\td_tri\\td_tet")
    values = "100\\t98\\t90\\t80\\t2\\t10\\t20"
    if settings["tie"]:
        values = "100\\t98\\t98\\t80\\t2\\t2\\t20"
    if settings["invalid_model"]:
        values = "nan\\tnan\\tnan\\tnan\\tnan\\tnan\\tnan"
    print(sys.argv[-1] + "\\t" + values)
else:
    raise SystemExit(3)
''')
        path.chmod(0o755)
        return path

    def run_fake(self, **settings):
        bam = self.root / "BAM path with spaces ' and $literal.bam"
        bam.write_bytes(b"fixture BAM delegated to fake executable")
        outdir = self.root / "evaluation"
        result = NQUIRE.run(self.fake_executable(**settings), bam, outdir, "S", 2)
        return result, outdir

    def test_bin_header_position_counts_and_extended_records(self):
        for extended in (False, True):
            path = self.root / "native.bin"
            path.write_bytes(bin_data(2, extended))
            result = NQUIRE.read_bin(path)
            self.assertEqual(result["positions"], 2)
            self.assertEqual(result["extended"], extended)
            self.assertEqual(result["minimum_allele_fraction"], .2)
            self.assertEqual(result["bam_header_hash"], 123456 if extended else None)

    def test_truncated_header_records_and_extra_records_are_rejected(self):
        for content in (b"", bin_data(1)[:8], bin_data(1)[:-1],
                        bin_data(0) + b"extra", struct.pack(">IBI", 1, 2, 200),
                        struct.pack(">IBI", 1, 0, 200) + struct.pack(">HBHH", 0, 2, 1, 1)):
            path = self.root / "bad.bin"
            path.write_bytes(content)
            with self.assertRaises(ValueError):
                NQUIRE.read_bin(path)

    def test_expected_model_has_minimum_delta(self):
        result = NQUIRE.parse_results(self.result(), "S-denoised.bin")
        self.assertEqual(result["best_supported_fixed_model"], "2")
        self.assertEqual(result["status"], "ranked_fixed_models")

    def test_header_only_extra_rows_wrong_label_and_nonfinite_values_fail(self):
        valid = self.result().read_text()
        cases = ["\t".join(NQUIRE.HEADER) + "\n", valid + valid.splitlines()[1] + "\n",
                 valid.replace("S-denoised.bin", "other-denoised.bin"),
                 valid.replace("100\t98", "nan\t98"), valid.replace("100\t98", "inf\t98"),
                 valid.replace("d_dip", "wrong_header"), valid.replace("\t2\t10\t20", "\t3\t10\t20")]
        for content in cases:
            with self.subTest(content=content), self.assertRaises(ValueError):
                NQUIRE.parse_results(self.result(content), "S-denoised.bin")

    def test_zero_raw_positions_skip_both_em_steps_and_publish_na_status(self):
        summary, output = self.run_fake(raw_positions=0)
        self.assertEqual(summary["status"], "untestable_zero_raw_positions")
        self.assertIsNone(summary["best_supported_fixed_model"])
        self.assertFalse(summary["lrdmodel_executed"])
        self.assertEqual((output / "calls.txt").read_text(), "create\n")
        self.assertFalse((output / "S-denoised.bin").exists())
        self.assertEqual((output / "S_ploidy.tsv").read_text(), "S\tNA\n")
        self.assertIn("NOT_RUN", (output / "S_nquire.results").read_text())
        self.assertTrue((output / "nquire.done").exists())

    def test_zero_denoised_positions_skip_model_and_publish_na_status(self):
        summary, output = self.run_fake(denoised_positions=0)
        self.assertEqual(summary["status"], "untestable_zero_denoised_positions")
        self.assertEqual((output / "calls.txt").read_text(), "create\ndenoise\n")
        self.assertIsNone(summary["fixed_model_deltas"])
        self.assertEqual((output / "S_ploidy.tsv").read_text(), "S\tNA\n")

    def test_ties_are_explicitly_ambiguous_with_na_model(self):
        summary, output = self.run_fake(tie=True)
        self.assertEqual(summary["status"], "ambiguous_tie")
        self.assertEqual(summary["minimum_delta_models"], ["2", "3"])
        self.assertIsNone(summary["best_supported_fixed_model"])
        self.assertEqual((output / "S_ploidy.tsv").read_text(), "S\tNA\n")

    def test_positive_native_counts_and_reports_are_validated_and_hashed(self):
        summary, output = self.run_fake()
        self.assertEqual((summary["raw_bin"]["positions"], summary["denoised_bin"]["positions"]), (3, 2))
        self.assertEqual((output / "S_ploidy.tsv").read_text(), "S\t2\n")
        self.assertTrue(summary["lrdmodel_executed"])
        self.assertEqual(summary["state"], "completed")
        self.assertEqual(summary["implementation"]["sha256"], NQUIRE.sha256(SCRIPT))
        receipt = json.loads((output / "S_nquire.summary.json").read_text())
        for name, artifact in receipt["outputs"].items():
            self.assertEqual(NQUIRE.sha256(output / name), artifact["sha256"])
        reread = NQUIRE.validate_analysis(output / "S_nquire.results", "S", output / "S.bin",
                                         output / "S-denoised.bin")
        self.assertEqual(reread["best_supported_fixed_model"], "2")

    def test_failed_rerun_archives_previous_outputs_and_keeps_diagnostics(self):
        output = self.root / "evaluation"
        output.mkdir()
        (output / "nquire.done").touch()
        (output / "old.results").write_text("retained previous report\n")
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_fake(fail_step="denoise")
        self.assertFalse((output / "nquire.done").exists())
        previous = list(self.root.glob(".evaluation.previous.*"))
        self.assertEqual(len(previous), 1)
        self.assertTrue((previous[0] / "nquire.done").exists())
        failed = list(self.root.glob(".S.nquire.*"))
        self.assertEqual(len(failed), 1)
        self.assertIn("forced native failure", (failed[0] / "S.nquire.log").read_text())
        self.assertEqual(json.loads((failed[0] / "failure.json").read_text())["state"], "failed_incomplete")

    def test_native_nan_and_increased_denoised_count_cannot_complete(self):
        for settings in ({"invalid_model": True}, {"denoised_positions": 4}):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                self.run_fake(**settings)
            self.assertFalse((self.root / "evaluation/nquire.done").exists())


@unittest.skipUnless(NATIVE_NQUIRE and SAMTOOLS,
                     "Set NQUIRE_TEST_EXECUTABLE and provide samtools for the native empty-BAM regression")
class NativeNquireTests(unittest.TestCase):
    def test_empty_native_bam_never_calls_em_or_claims_diploidy(self):
        with tempfile.TemporaryDirectory(prefix="hifi-nquire-native-test-") as temporary:
            root = Path(temporary)
            sam, bam = root / "empty.sam", root / "empty.bam"
            sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:fixture\tLN:1000\n")
            subprocess.run([SAMTOOLS, "view", "-b", "-o", str(bam), str(sam)], check=True,
                           capture_output=True, timeout=10)
            output = root / "nQuire"
            result = subprocess.run([sys.executable, str(SCRIPT), "--executable", NATIVE_NQUIRE,
                                     "--bam", str(bam), "--outdir", str(output), "--sample", "EMPTY",
                                     "--threads", "1"], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            receipt = json.loads((output / "EMPTY_nquire.summary.json").read_text())
            self.assertEqual(receipt["raw_bin"]["positions"], 0)
            self.assertEqual([command["step"] for command in receipt["commands"]], ["create"])
            self.assertEqual(receipt["status"], "untestable_zero_raw_positions")
            self.assertEqual((output / "EMPTY_ploidy.tsv").read_text(), "EMPTY\tNA\n")
            self.assertNotIn("problem allocating Data in EM", (output / "EMPTY.nquire.log").read_text())


if __name__ == "__main__":
    unittest.main()
