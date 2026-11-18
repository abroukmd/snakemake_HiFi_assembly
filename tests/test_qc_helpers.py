"""Regression tests for QC provenance and score interpretation."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "workflow/scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


genomescope = load("parse_genomescope2")
merqury = load("run_merqury")
PLOT_FIXTURE = '''spectra_cn_plot <- function(zero, y_max) {
  dat_0 = NULL
  if (zero != "") {
    dat_0=read.table(zero, header=FALSE)
    dat_0[,1]=as.factor(dat_0[,1])
    dat_0[,1]=factor(dat_0[,1], levels=rev(unique(dat_0[,1])), ordered=TRUE)
  }
  if (zero != "") {
    y_max=max(y_max, sum(dat_0[,3]*1.1))
  }
}
'''


class GenomeScopeTests(unittest.TestCase):
    def test_rank_uses_summary_fit_not_console_adups(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "models.log"
            log.write_text("GenomeScope analyzing test p=1 k=21 outdir=p1\n"
                           "Model converged het:0 kcov:20 err:0.1 model fit:0.9 len:100\n"
                           "[INFO] GenomeScope: p=2 k=21 histo=test outdir=p2\n"
                           "Model converged het:0 kcov:20 err:0.1 model fit:0.1 len:100\n")
            for p, fit in [(1, "98% 99%"), (2, "80% 90%")]:
                (root / f"p{p}").mkdir()
                (root / f"p{p}/summary.txt").write_text(
                    f"Genome Haploid Length 100 bp 110 bp\nModel Fit {fit}\n")
            rows = genomescope.summarize([log])
            self.assertEqual([r["best_by_model_fit"] for r in rows], ["YES", ""])
            self.assertEqual(rows[0]["model_fit_all_percent"], 98)
            self.assertNotIn("best", rows[0])
            self.assertIn("not_validated_ploidy", rows[0]["selection_status"])

    def test_console_only_and_failed_summaries_do_not_succeed(self):
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "models.log"
            log.write_text("GenomeScope analyzing test p=1 k=21 outdir=failed\n"
                           "Model converged het:0 kcov:20 err:0.1 model fit:0.1 len:100\n")
            with self.assertRaisesRegex(ValueError, "No converged"):
                genomescope.summarize([log])


class MerquryTests(unittest.TestCase):
    def write_reports(self, root, qv="assembly.final.hifiasm\t1\t100\t33\t0.001\n"):
        (root / "merqury.qv").write_text(qv)
        (root / "completeness.stats").write_text("assembly.final.hifiasm\tall\t90\t100\t90\n")
        (root / "merqury.spectra-asm.ln.png").write_bytes(b"\x89PNG\r\n\x1a\nfixture")
        (root / "logs").mkdir(exist_ok=True)

    def test_duplicate_basename_rows_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = "assembly.final.hifiasm\t1\t100\t33\t0.001\n"
            self.write_reports(root, row * 2)
            with self.assertRaisesRegex(ValueError, "exactly"):
                merqury.validate_reports(root, "assembly.final.hifiasm")

    def test_wrong_assembly_or_hidden_plot_failure_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_reports(root)
            with self.assertRaises(ValueError):
                merqury.validate_reports(root, "assembly.final.lja")
            (root / "logs/run.log").write_text("Execution halted\n")
            with self.assertRaisesRegex(ValueError, "nested tool failure"):
                merqury.validate_reports(root, "assembly.final.hifiasm")

    def test_hidden_meryl_option_failure_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_reports(root)
            (root / "logs/run.log").write_text("Unknown option '-dump'\n")
            with self.assertRaisesRegex(ValueError, "nested tool failure"):
                merqury.validate_reports(root, "assembly.final.hifiasm")

    def test_native_prefixed_completeness_is_retained_and_conflicts_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            native = root / "merqury.completeness.stats"
            native.write_text("S.raw.hifiasm\tall\t90\t100\t90\n")
            merqury.normalize_completeness(root)
            self.assertEqual((root / "completeness.stats").read_bytes(), native.read_bytes())
            (root / "completeness.stats").write_text("conflicting older report\n")
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                merqury.normalize_completeness(root)

    def test_unknown_plot_source_fails_without_modifying_installation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            installed, work = root / "installed", root / "work"
            source = installed / "plot/plot_spectra_cn.R"
            source.parent.mkdir(parents=True)
            source.write_text("# unexpected future or modified source\n")
            work.mkdir()
            before = source.read_bytes()
            with self.assertRaisesRegex(ValueError, "Unsupported Merqury plot source SHA256"):
                merqury.prepare_runtime(installed, work)
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(list(work.iterdir()), [])

    def test_fresh_workspace_namespaces_assembly_and_read_database(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tools = root / "bin"
            tools.mkdir()
            reads = root / "S.meryl"
            reads.mkdir()
            assembly = root / "S.fasta"
            assembly.write_text(">contig\nACGT\n")
            meryl = tools / "meryl"
            meryl.write_text("#!/bin/sh\nexit 0\n")
            meryl.chmod(0o755)
            installed = root / "share/merqury"
            for name in ("util/util.sh", "eval/spectra-cn.sh", "plot/plot_spectra_cn.R"):
                source = installed / name
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_text("# fixture\n")
            plot_source = installed / "plot/plot_spectra_cn.R"
            plot_source.write_text(PLOT_FIXTURE)
            original_plot_sha256 = merqury.sha256(plot_source)
            executable = installed / "merqury.sh"
            executable.write_text("#!" + sys.executable + "\n" + '''import os, pathlib, sys
runtime = pathlib.Path(os.environ["MERQURY"])
assert runtime.name == "merqury-runtime"
assert (runtime / "merqury.sh").samefile(__file__)
assert "dat_0[,3] > 0" in (runtime / "plot/plot_spectra_cn.R").read_text()
reads, assembly, prefix = sys.argv[1:]
assert pathlib.Path(reads).name == "illumina.reads.meryl"
assert pathlib.Path(assembly).name == "S.raw.hifiasm.fasta"
label=pathlib.Path(assembly).stem
pathlib.Path(prefix+".qv").write_text(label+"\\t1\\t100\\t33\\t0.001\\n")
pathlib.Path("completeness.stats").write_text(label+"\\tall\\t90\\t100\\t90\\n")
pathlib.Path(prefix+".spectra-asm.ln.png").write_bytes(b"\\x89PNG\\r\\n\\x1a\\nfixture")
pathlib.Path("logs").mkdir()
pathlib.Path("assembly-only.hist").write_text("1\\t0\\t1\\n2\\t0\\t0\\n")
''')
            executable.chmod(0o755)
            (tools / "merqury.sh").symlink_to(executable)
            outdir = root / "evaluation"
            outdir.mkdir()
            (outdir / "stale.meryl").mkdir()
            env = os.environ.copy()
            env["PATH"] = str(tools) + os.pathsep + env["PATH"]
            env["MERQURY"] = str(root / "unrelated-host-installation")
            # Change the source pin only inside the test interpreter. Production
            # has no environment variable or CLI switch to bypass its guard.
            bridge = '''import importlib.util, sys
from unittest.mock import patch
script, expected_hash = sys.argv[1:3]
sys.argv = [script] + sys.argv[3:]
spec = importlib.util.spec_from_file_location("run_merqury_fixture", script)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
with patch.object(module, "PINNED_PLOT_SHA256", expected_hash):
    module.main()
'''
            result = subprocess.run([sys.executable, "-c", bridge,
                                     str(ROOT / "workflow/scripts/run_merqury.py"), original_plot_sha256,
                                     "--assembly", str(assembly), "--read-db", str(reads),
                                     "--outdir", str(outdir), "--label", "S.raw.hifiasm",
                                     "--threads", "2"], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((outdir / "stale.meryl").exists())
            self.assertTrue((outdir / "S.raw.hifiasm.fasta").samefile(assembly))
            self.assertTrue((outdir / "illumina.reads.meryl").samefile(reads))
            self.assertTrue((outdir / "inputs.json").is_file())
            import json
            receipt = json.loads((outdir / "inputs.json").read_text())
            self.assertEqual(receipt["assembly_sha256"], merqury.sha256(assembly))
            self.assertEqual(receipt["tools"]["merqury_root"], str(installed))
            self.assertIn(str(executable), receipt["tools"]["source_sha256"])
            adaptation = receipt["tools"]["plot_adaptation"]
            self.assertEqual(adaptation["original_sha256"], original_plot_sha256)
            self.assertEqual(merqury.sha256(plot_source), original_plot_sha256)
            patched = outdir / adaptation["patched_path"]
            self.assertFalse(patched.is_symlink())
            self.assertEqual(merqury.sha256(patched), adaptation["patched_sha256"])
            self.assertIn("dat_0=NULL", patched.read_text())
            self.assertIn("if (!is.null(dat_0))", patched.read_text())
            self.assertEqual((outdir / "assembly-only.hist").read_text(), "1\t0\t1\n2\t0\t0\n")
            self.assertEqual((outdir / "merqury.qv").read_text(), "S.raw.hifiasm\t1\t100\t33\t0.001\n")


if __name__ == "__main__":
    unittest.main()
