"""Check MGA input preservation and its native-output contract without assembly."""

import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


spec = importlib.util.spec_from_file_location(
    "run_mga", Path(__file__).resolve().parents[1] / "workflow/scripts/run_mga.py")
mga = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mga)


class MgaHelperTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="hifi-mga-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_gzipped_fastq_strips_metadata_preserving_bases_and_qualities(self):
        reads = self.root / "reads.fastq.gz"
        with gzip.open(reads, "wt") as handle:
            handle.write("@read1 ML:B:C,1,2\nAcgTN\n+read1 ML:B:C,1,2\n!+#@I\n")
        output = self.root / "clean.fastq"
        stats = mga.clean_reads(reads, output)
        self.assertEqual(output.read_text(), "@read1\nAcgTN\n+\n!+#@I\n")
        self.assertEqual((stats["reads"], stats["bases"], stats["headers_trimmed"]), (1, 5, 1))

    def test_wrapped_fasta_and_fastq_preserve_complete_sequences(self):
        reads = self.root / "reads"
        output = self.root / "clean"
        reads.write_text(">read1 comment\nAcg\nTn\n>read2\nGG\n")
        self.assertEqual(mga.clean_reads(reads, output)["format"], "fasta")
        self.assertEqual(output.read_text(), ">read1\nAcgTn\n>read2\nGG\n")
        reads.write_text("@read1 comment\nAcg\nTn\n+\n!+#\n@I\n")
        mga.clean_reads(reads, output)
        self.assertEqual(output.read_text(), "@read1\nAcgTn\n+\n!+#@I\n")

    def test_duplicate_trimmed_ids_and_truncated_quality_fail(self):
        reads, output = self.root / "reads", self.root / "clean"
        for content, error in [(
                ">read1 first\nAAA\n>read1 second\nCCC\n", "Duplicate read identifier"),
                ("@read1\nACGT\n+\n!!!\n", "Truncated FASTQ quality"),
                ("@read1\nA\n+\n!\n\n@read2\nA\n+\n!\n", "Unexpected blank line")]:
            reads.write_text(content)
            with self.assertRaisesRegex(ValueError, error):
                mga.clean_reads(reads, output)

    def fake_executable(self, produce_polished=True, fail_exit=False):
        script = self.root / "MGA"
        script.write_text("#!/usr/bin/env python3\n"
                          "import argparse, json, os\nfrom pathlib import Path\n"
                          "p=argparse.ArgumentParser()\n"
                          "p.add_argument('--reads')\np.add_argument('--output')\n"
                          "p.add_argument('--threads')\na=p.parse_args()\n"
                          "out=Path(a.output)\nout.mkdir()\n"
                          "(out/'native_affinity.json').write_text(json.dumps(sorted(os.sched_getaffinity(0))))\n"
                          "(out/'graph.gfa').write_text('H\\tVN:Z:1.0\\n')\n"
                          "(out/'graph.link').symlink_to(out/'graph.gfa')\n"
                          "(out/'reads.link').symlink_to(Path(a.reads))\n"
                          + ("scratch=Path(a.reads).parent/'mapping.temporary'\n"
                             "scratch.write_text('partial mapping evidence\\n')\n"
                             "(out/'scratch.link').symlink_to(scratch)\n"
                             "raise SystemExit(7)\n" if fail_exit else "")
                          + ("(out/'5_polishing').mkdir()\n"
                             "(out/'5_polishing/assembly.fasta').write_text('>native\\nAcgTN\\n')\n"
                             if produce_polished else ""))
        script.chmod(0o755)
        return script

    def test_run_uses_native_polished_output_and_keeps_graph_stages(self):
        reads = self.root / "source.fastq"
        reads.write_text("@read1 metadata\nACGT\n+\n!!!!\n")
        output = self.root / "output"
        (output / "5_polishing").mkdir(parents=True)
        receipt = mga.run(self.fake_executable(), reads, output, 2, "a" * 40)
        self.assertEqual((output / "assembly.fasta").read_bytes(),
                         (output / "5_polishing/assembly.fasta").read_bytes())
        self.assertTrue((output / "graph.gfa").exists())
        self.assertEqual((output / "graph.link").read_bytes(), (output / "graph.gfa").read_bytes())
        self.assertFalse(Path(mga.os.readlink(output / "graph.link")).is_absolute())
        self.assertEqual((output / "reads.link").read_text(), "@read1\nACGT\n+\n!!!!\n")
        self.assertTrue(receipt["graph_links"]["prepared_reads_retained"])
        self.assertEqual((output / "mga.done").read_text().strip(), "a" * 40)
        self.assertEqual(receipt["input_reads"]["headers_trimmed"], 1)
        command = json.loads((output / "mga_run.json").read_text())["command"]
        self.assertTrue(command[command.index("--reads") + 1].endswith(".fastq"))

    def test_missing_polished_output_fails_and_preserves_diagnostics(self):
        reads = self.root / "source.fasta"
        reads.write_text(">read1\nACGT\n")
        output = self.root / "output"
        with self.assertRaises(FileNotFoundError):
            mga.run(self.fake_executable(produce_polished=False), reads, output, 2, "a" * 40)
        self.assertFalse((output / "mga.done").exists())
        failed = list(self.root.glob(".mga-failed-*"))
        self.assertEqual(len(failed), 1)
        self.assertTrue((failed[0] / "native/graph.gfa").exists())

    def test_failed_native_absolute_links_keep_graph_reads_and_scratch_targets(self):
        reads = self.root / "source.fasta"
        reads.write_text(">read1 metadata\nACGT\n")
        output = self.root / "output"
        with self.assertRaises(mga.subprocess.CalledProcessError):
            mga.run(self.fake_executable(fail_exit=True), reads, output, 2, "a" * 40)
        failed, = self.root.glob(".mga-failed-*")
        self.assertEqual((failed / "native/graph.link").read_bytes(),
                         (failed / "native/graph.gfa").read_bytes())
        self.assertEqual((failed / "native/reads.link").read_text(), ">read1\nACGT\n")
        self.assertEqual((failed / "native/scratch.link").read_text(), "partial mapping evidence\n")
        for name in ["graph.link", "reads.link", "scratch.link"]:
            self.assertFalse(Path(mga.os.readlink(failed / "native" / name)).is_absolute())
        receipt = json.loads((failed / "mga_run.json").read_text())
        self.assertEqual(len(receipt["failed_workspace_links"]), 3)
        self.assertTrue(all(item["target_existed_before_preservation"]
                            for item in receipt["failed_workspace_links"]))
        self.assertFalse((output / "mga.done").exists())

    def test_existing_output_is_preserved_and_refused_before_execution(self):
        reads = self.root / "source.fasta"
        reads.write_text(">read1\nACGT\n")
        output = self.root / "output"
        output.mkdir()
        (output / "graph.gfa").write_text("existing graph evidence\n")
        with mock.patch.object(mga.subprocess, "run") as execute:
            with self.assertRaises(FileExistsError):
                mga.run(self.fake_executable(), reads, output, 2, "a" * 40)
            execute.assert_not_called()
        self.assertEqual((output / "graph.gfa").read_text(), "existing graph evidence\n")

    def test_native_child_is_cpu_bounded_without_changing_parent_affinity(self):
        inherited = mga.os.sched_getaffinity(0)
        reads = self.root / "source.fasta"
        reads.write_text(">read1\nACGT\n")
        output = self.root / "output"
        receipt = mga.run(self.fake_executable(), reads, output, 1, "a" * 40)
        child = set(json.loads((output / "native_affinity.json").read_text()))
        self.assertLessEqual(len(child), 1)
        self.assertTrue(child <= inherited)
        self.assertEqual(child, set(receipt["cpu_allocation"]["selected_cpus"]))
        self.assertEqual(mga.os.sched_getaffinity(0), inherited)

    def test_missing_affinity_support_fails_clearly(self):
        with mock.patch.object(mga.os, "sched_getaffinity", None):
            with self.assertRaisesRegex(RuntimeError, "CPU affinity support"):
                mga.native_cpu_allocation(1)

    def test_dangling_and_unknown_scratch_links_are_rejected(self):
        scratch = self.root / "scratch"
        native = scratch / "native"
        native.mkdir(parents=True)
        prepared = scratch / "reads.fastq"
        prepared.write_text("@r\nA\n+\n!\n")
        link = native / "invalid"
        link.symlink_to(scratch / "missing")
        with self.assertRaisesRegex(ValueError, "dangling or cyclic"):
            mga.prepare_output_links(native, scratch, prepared)
        link.unlink()
        other = scratch / "temporary.txt"
        other.write_text("temporary evidence\n")
        link.symlink_to(other)
        with self.assertRaisesRegex(ValueError, "unpreserved scratch artifact"):
            mga.prepare_output_links(native, scratch, prepared)


if __name__ == "__main__":
    unittest.main()
