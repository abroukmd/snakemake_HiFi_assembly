"""Construct the optional MGA branch without installing or executing MGA."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake runtime required")
class MgaDagTests(unittest.TestCase):
    def test_raw_reads_target_requires_pinned_install_and_native_polishing(self):
        repo = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory(prefix="hifi-mga-dag-") as tmp:
            root = Path(tmp)
            reads = root / "reads.fasta"
            reads.write_text(">fixture\nACGT\n")
            config = root / "config.json"
            config.write_text(json.dumps({"project_dir": str(root / "project")}))
            snakefile = root / "Snakefile"
            snakefile.write_text(
                "import os\n"
                "project_dir = config['project_dir']\n"
                "def outpath(*parts):\n    return os.path.join(project_dir, *parts)\n"
                f"sample_paths = {{'TEST': {str(reads)!r}}}\n"
                "rule all:\n"
                "    input: outpath('Assemblies/TEST/MGA/assembly.fasta')\n"
                f"include: {str(repo / 'workflow/rules/mga.smk')!r}\n")
            result = subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile),
                 "--configfile", str(config), "--dry-run", "--cores", "4", "--nocolor"],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertIn("rule install_mga:", result.stdout)
            self.assertIn("rule mga:", result.stdout)
            self.assertIn("5_polishing/assembly.fasta", result.stdout)
            self.assertIn("f0f26da124e5e1362875dbf15c0de1bec2fd34f8", result.stdout)
            self.assertNotIn("rule lja:", result.stdout)


if __name__ == "__main__":
    unittest.main()
