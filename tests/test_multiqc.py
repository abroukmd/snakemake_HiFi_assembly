"""Ensure BUSCO alias directories cannot duplicate MultiQC lineage reports."""
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


RULE = Path(__file__).resolve().parents[1] / "workflow/rules/multiqc.smk"
MULTIQC = os.environ.get("MULTIQC_TEST_EXECUTABLE") or shutil.which("multiqc")


@unittest.skipUnless(MULTIQC, "Set MULTIQC_TEST_EXECUTABLE to run native BUSCO discovery regression")
class MultiqcBuscoDiscoveryTests(unittest.TestCase):
    def test_top_level_lineages_survive_without_real_or_symlink_run_duplicates(self):
        with tempfile.TemporaryDirectory(prefix="hifi-multiqc-test-") as temporary:
            root = Path(temporary)
            qc = root / "Assemblies/TEST/QC"
            result = qc / "BUSCO-hifiasm/TEST-hifiasm"
            result.mkdir(parents=True)

            def summary(lineage):
                return ("# BUSCO version is: 5.8.3\n"
                        f"# The lineage dataset is: {lineage} (Creation date: 2025-07-01)\n"
                        "10\tComplete BUSCOs (C)\n9\tComplete and single-copy BUSCOs (S)\n"
                        "1\tComplete and duplicated BUSCOs (D)\n0\tFragmented BUSCOs (F)\n"
                        "0\tMissing BUSCOs (M)\n10\tTotal BUSCO groups searched\n")

            expected = {}
            for kind, lineage in [("generic", "eukaryota_odb12"), ("specific", "fungi_odb12")]:
                report = result / f"short_summary.{kind}.{lineage}.TEST-hifiasm.txt"
                report.write_text(summary(lineage))
                expected[report.resolve()] = lineage
            auto = result / "auto_lineage/run_eukaryota_odb12"
            auto.mkdir(parents=True)
            (auto / "short_summary.txt").write_text(summary("eukaryota_odb12"))
            (result / "run_eukaryota_odb12").symlink_to(auto, target_is_directory=True)
            native = result / "run_fungi_odb12"
            native.mkdir()
            (native / "short_summary.txt").write_text(summary("fungi_odb12"))

            patterns = re.findall(r"--ignore '([^']+)'", RULE.read_text())
            self.assertTrue(patterns, "The rule's scan exclusion arguments were not found")
            environment = dict(os.environ, MPLCONFIGDIR=str(root / "matplotlib"))

            def scan(name, ignored):
                output = root / name
                command = [MULTIQC, str(qc), "--module", "busco", "--dirs", "--dirs-depth", "3",
                           "--no-report", "--no-version-check", "--force", "-o", str(output)]
                for pattern in ignored:
                    command.extend(["--ignore", pattern])
                completed = subprocess.run(command, env=environment, text=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
                self.assertEqual(completed.returncode, 0, completed.stdout)
                with (output / "multiqc_data/multiqc_sources.txt").open() as handle:
                    sources = [row for row in csv.DictReader(handle, delimiter="\t") if row["Module"] == "BUSCO"]
                with (output / "multiqc_data/multiqc_busco.txt").open() as handle:
                    data = list(csv.DictReader(handle, delimiter="\t"))
                return sources, data

            # Establish that this fixture reproduces the old alias duplication.
            old_sources, _ = scan("control", ["*/auto_lineage/*"])
            self.assertEqual(len(old_sources), 4)
            sources, data = scan("fixed", patterns)
            self.assertEqual(len(sources), 2)
            self.assertEqual({Path(row["Source"]).resolve() for row in sources}, set(expected))
            self.assertTrue(all(Path(row["Source"]).parent == result for row in sources))
            self.assertEqual({row["lineage_dataset"] for row in data}, set(expected.values()))
            self.assertEqual(len(data), 2)


if __name__ == "__main__":
    unittest.main()
