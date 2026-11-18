"""Exercise the runner contract without executing expensive alignment tools."""
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


RUNNER = Path(__file__).resolve().parents[1] / "workflow/scripts/run_chromeister.sh"
PNG_SIGNATURE = bytes.fromhex("89504e470d0a1a0a")


class ChromeisterRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        tools = self.root / "bin"
        tools.mkdir()
        programs = {
            "CHROMEISTER": '''#!/usr/bin/env python3
import os, pathlib, sys, time
a=sys.argv
ref=pathlib.Path(a[a.index('-query')+1])
content=ref.read_text()
time.sleep(0.1)
assert ref.read_text() == content
if os.environ.get('MOCK_NO_MATRIX'):
    sys.exit(0)
pathlib.Path(a[a.index('-out')+1]).write_text(str(ref)+'\\n'+content)
''',
            "compute_score.R": "#!/bin/sh\nexit 0\n",
            "Rscript": '''#!/usr/bin/env python3
import os, pathlib, sys
if os.environ.get('MOCK_NO_PNG'):
    print('score: 0.9')
    sys.exit(0)
png=sys.argv[3] if sys.argv[1]=='-e' else sys.argv[2]+'.filt.png'
data=b'not a PNG' if os.environ.get('MOCK_BAD_PNG') else bytes.fromhex('89504e470d0a1a0a')+b'fixture'
pathlib.Path(png).write_bytes(data)
print('score: 0.9')
''',
        }
        for name, source in programs.items():
            path = tools / name
            path.write_text(source)
            path.chmod(0o755)
        self.env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"])
        self.query = self.root / "assembly.fa"
        self.query.write_text(">assembly\nACGT\n")
        self.ref = self.root / "reference.fna.gz"
        with gzip.open(self.ref, "wt") as out:
            out.write(">reference\nACGT\n")

    def command(self, label):
        prefix = self.root / label
        return ["bash", str(RUNNER), str(self.query), str(self.ref),
                str(prefix) + ".mat", str(prefix) + ".mat.filt.png",
                str(prefix) + ".score", str(prefix) + ".done",
                str(prefix) + ".log", str(prefix) + ".err", "0"]

    def test_parallel_jobs_use_separate_references_and_real_images(self):
        processes = [subprocess.Popen(self.command(label), env=self.env)
                     for label in ("hifiasm", "lja")]
        self.assertEqual([p.wait(timeout=10) for p in processes], [0, 0])
        prepared_refs = []
        for label in ("hifiasm", "lja"):
            self.assertTrue((self.root / f"{label}.done").exists())
            png = (self.root / f"{label}.mat.filt.png").read_bytes()
            self.assertTrue(png.startswith(PNG_SIGNATURE))
            prepared_refs.append((self.root / f"{label}.mat").read_text().splitlines()[0])
        self.assertNotEqual(*prepared_refs)
        self.assertTrue(self.ref.exists())
        self.assertEqual(list(self.root.glob(".chromeister.*")), [])

    def test_invalid_image_does_not_leave_success_marker(self):
        marker = self.root / "bad.done"
        marker.touch()
        env = dict(self.env, MOCK_BAD_PNG="1")
        result = subprocess.run(self.command("bad"), env=env, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(marker.exists())
        self.assertIn("Invalid PNG", (self.root / "bad.err").read_text())

    def test_stale_matrix_cannot_replace_missing_new_alignment(self):
        command = self.command("stale_matrix")
        subprocess.run(command, env=self.env, timeout=10, check=True)
        result = subprocess.run(command, env=dict(self.env, MOCK_NO_MATRIX="1"), timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / "stale_matrix.done").exists())
        self.assertFalse((self.root / "stale_matrix.mat").exists())
        self.assertFalse((self.root / "stale_matrix.mat.filt.png").exists())
        self.assertIn("fresh matrix", (self.root / "stale_matrix.err").read_text())

    def test_stale_image_cannot_replace_missing_new_plot(self):
        command = self.command("stale_image")
        subprocess.run(command, env=self.env, timeout=10, check=True)
        result = subprocess.run(command, env=dict(self.env, MOCK_NO_PNG="1"), timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / "stale_image.done").exists())
        self.assertFalse((self.root / "stale_image.mat.filt.png").exists())
        self.assertIn("fresh PNG", (self.root / "stale_image.err").read_text())

    def test_optional_empty_haplotigs_are_explicitly_labelled(self):
        self.query.write_text(">hap_empty\nN\n")
        command = self.command("empty")
        command[-1] = "1"
        subprocess.run(command, env=self.env, timeout=10, check=True)
        self.assertTrue((self.root / "empty.done").exists())
        self.assertEqual((self.root / "empty.score").read_text(), "SKIPPED hap_empty\n")
        self.assertTrue((self.root / "empty.mat").read_text().startswith("# SKIPPED"))
        self.assertTrue((self.root / "empty.mat.filt.png").read_bytes().startswith(PNG_SIGNATURE))


if __name__ == "__main__":
    unittest.main()
