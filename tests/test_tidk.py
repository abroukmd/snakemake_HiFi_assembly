"""Known-answer checks for the optional exact-repeat TIDK module."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/tidk_telomeres.py"
SPEC = importlib.util.spec_from_file_location("tidk_telomeres", SCRIPT)
TIDK_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TIDK_MODULE)
TIDK_EXECUTABLE = os.environ.get("TIDK_TEST_EXECUTABLE") or shutil.which("tidk")


class RepeatArrayTests(unittest.TestCase):
    def test_primitive_rotation_and_reverse_complement_share_one_candidate(self):
        self.assertEqual(TIDK_MODULE.primitive("TTAGGGTTAGGG"), "TTAGGG")
        self.assertEqual(TIDK_MODULE.variants("TTAGGG"), TIDK_MODULE.variants("CCCTAA"))

    def test_end_window_crossing_counts_clipped_density(self):
        sequence = "AC" * 4 + "TTAGGG" * 6 + "AC" * 20
        spans, left, right, terminal, internal = TIDK_MODULE.analyse(sequence, "TTAGGG", 5, 12)
        self.assertEqual(spans, [[8, 44]])
        self.assertEqual(left, [8, 44])
        self.assertIsNone(right)
        self.assertEqual(terminal, 4)
        self.assertEqual(internal, 32)

    def test_short_sequence_windows_do_not_overlap(self):
        sequence = "TTAGGG" * 6
        end_size = min(10000, len(sequence) // 2)
        spans, left, right, terminal, internal = TIDK_MODULE.analyse(sequence, "TTAGGG", 5, end_size)
        self.assertEqual(spans, [[0, 36]])
        self.assertEqual(left, right)
        self.assertEqual(terminal, len(sequence))
        self.assertEqual(internal, 0)

    def test_discovery_collapses_equivalent_motifs(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "explore.tsv"
            source.write_text("canonical_repeat_unit\tcount_repeat_runs_gt_20\nTTAGGG\t8\nCCCTAA\t7\nTAGGGT\t9\n")
            candidates = TIDK_MODULE.parse_candidates(source, 10)
            self.assertEqual(len(candidates), 1)
            self.assertEqual(candidates[0]["explore_score"], 9)

    def test_tool_failure_preserves_diagnostics_without_publishing_completion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "fixture.fa"
            fasta.write_text(">sequence\n" + "TTAGGG" * 30 + "\n")
            executable = root / "tidk"
            executable.write_text('#!/bin/sh\nif [ "$1" = --version ]; then echo "tidk 0.2.7"; exit 0; fi\necho forced-explore-failure >&2\nexit 7\n')
            executable.chmod(0o755)
            output = root / "screen"
            result = subprocess.run([sys.executable, str(SCRIPT), str(fasta), "--tidk", str(executable),
                                     "-o", str(output)], text=True, capture_output=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((output / "complete.json").exists())
            staging = list(root.glob(".screen.tidk.*"))
            self.assertEqual(len(staging), 1)
            self.assertIn(str(staging[0]), result.stderr)
            self.assertIn("forced-explore-failure", (staging[0] / "tidk_explore.log").read_text())


@unittest.skipUnless(TIDK_EXECUTABLE, "Set TIDK_TEST_EXECUTABLE to test the installed 0.2.7 binary")
class NativeTidkTests(unittest.TestCase):
    def test_terminal_fixture_is_published_with_hashes_and_explicit_status(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "fixture.fa"
            # A deterministic non-telomeric interior and exact repeats at both ends.
            fasta.write_text(">chromosome\n" + "TTAGGG" * 30 + "ACGTCAGT" * 200 + "CCCTAA" * 30 + "\n"
                             ">mito_hint_only\nACGTCAGTACGTCAGTACGTCAGT\n")
            output = root / "screen"
            environment = dict(os.environ, RAYON_NUM_THREADS="1")
            result = subprocess.run([sys.executable, str(SCRIPT), str(fasta), "--tidk", TIDK_EXECUTABLE,
                                     "-o", str(output), "-e", "200", "-d", "0.2", "-n", "10"],
                                    env=environment, text=True, capture_output=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = list(csv.DictReader((output / "assembly_summary.tsv").read_text().splitlines(), delimiter="\t"))[0]
            self.assertEqual(summary["screen_status"], "CANDIDATE_IDENTIFIED")
            self.assertEqual(summary["retained_sequences"], "2")
            self.assertEqual(summary["both_ends"], "1")
            inclusion = list(csv.DictReader((output / "sequence_inclusion.tsv").read_text().splitlines(), delimiter="\t"))
            self.assertEqual(inclusion[1]["organellar_name_hint"], "1")
            self.assertEqual(inclusion[1]["included"], "1")
            config = json.loads((output / "run_config.json").read_text())
            self.assertEqual(config["input_sha256"], TIDK_MODULE.sha256(fasta))
            self.assertEqual(config["rayon_threads"], "1")
            receipt = json.loads((output / "complete.json").read_text())
            for filename, expected in receipt["output_sha256"].items():
                self.assertEqual(TIDK_MODULE.sha256(output / filename), expected)
            self.assertEqual(list(root.glob(".screen.tidk.*")), [])
            # Existing results are protected rather than silently overwritten.
            duplicate = subprocess.run([sys.executable, str(SCRIPT), str(fasta), "--tidk", TIDK_EXECUTABLE,
                                        "-o", str(output)], text=True, capture_output=True, timeout=60)
            self.assertNotEqual(duplicate.returncode, 0)
            self.assertIn("not empty", duplicate.stderr)


if __name__ == "__main__":
    unittest.main()
