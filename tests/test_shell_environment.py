"""Execute one real workflow rule with a controlled lightweight tool fixture."""

import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SNAKEFILE = Path(__file__).resolve().parents[1] / "workflow" / "Snakefile"


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake runtime required")
class ShellEnvironmentTests(unittest.TestCase):
    def test_rule_python_disables_user_site_before_running_tool(self):
        with tempfile.TemporaryDirectory(prefix="hifi-shell-environment-test-") as temporary:
            root = Path(temporary)
            reads = root / "reads.fastq.gz"
            with gzip.open(reads, "wt") as handle:
                handle.write("@fixture\nACGTACGT\n+\nIIIIIIII\n")
            sheet = root / "samples.tsv"
            with sheet.open("w") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "HiFi", "Ref", "fq1", "fq2"])
                writer.writerow(["TEST", *([str(reads)] * 4)])
            project = root / "project"
            config = root / "config.json"
            config.write_text(json.dumps({"project_dir": str(project), "sample_sheet": str(sheet),
                                          "oatk_db": str(reads)}))
            tools = root / "tools"
            tools.mkdir()
            nanoplot = tools / "NanoPlot"
            nanoplot.write_text(f"#!{sys.executable}\n" + '''import argparse, gzip, json, os, pathlib, sys
assert os.environ.get("PYTHONNOUSERSITE") == "1", "Rule did not set PYTHONNOUSERSITE=1"
assert sys.flags.no_user_site == 1, "Tool interpreter still permits user-site imports"
parser = argparse.ArgumentParser()
parser.add_argument("-o", type=pathlib.Path, required=True)
parser.add_argument("--fastq", required=True)
args, unused = parser.parse_known_args()
with gzip.open(args.fastq, "rt") as handle:
    assert handle.read().startswith("@fixture\\nACGTACGT\\n")
args.o.mkdir(parents=True, exist_ok=True)
(args.o / "NanoStats.txt").write_text("Number of reads: 1\\n")
(args.o / "NanoPlot-report.html").write_text("<html><body>fixture</body></html>\\n")
(args.o / "fixture_environment.json").write_text(json.dumps({"PYTHONNOUSERSITE": os.environ["PYTHONNOUSERSITE"], "no_user_site": sys.flags.no_user_site}))
''')
            nanoplot.chmod(0o755)
            environment = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ.get("PATH", ""))
            environment.pop("PYTHONNOUSERSITE", None)
            output = project / "Assemblies/TEST/QC/nanoplot"
            result = subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", str(SNAKEFILE),
                 "--configfile", str(config), "--cores", "1", "--nocolor", "--latency-wait", "1",
                 str(output / "nanoplot.done")],
                cwd=root, env=environment, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=90,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertTrue((output / "nanoplot.done").is_file(), result.stdout)
            self.assertGreater((output / "NanoStats.txt").stat().st_size, 0)
            self.assertGreater((output / "NanoPlot-report.html").stat().st_size, 0)
            receipt = json.loads((output / "fixture_environment.json").read_text())
            self.assertEqual(receipt, {"PYTHONNOUSERSITE": "1", "no_user_site": 1})


if __name__ == "__main__":
    unittest.main()
