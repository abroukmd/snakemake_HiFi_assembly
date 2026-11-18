"""Exercise assembler selection through the complete downstream DAG."""

import csv
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


SNAKEFILE = Path(__file__).resolve().parents[1] / "workflow" / "Snakefile"


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake runtime required")
class WorkflowDagTests(unittest.TestCase):
    def check_selection(self, hifiasm, lja, ragtag, mga="N", optional_qc=False):
        with tempfile.TemporaryDirectory(prefix="hifi-dag-test-") as tmp:
            root = Path(tmp)
            source = root / "input.fa"
            source.write_text(">fixture\nACGTACGT\n")
            sheet = root / "samples.tsv"
            with sheet.open("w") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "HiFi", "Ref", "fq1", "fq2", "ragtag", "run_hifiasm", "run_lja", "run_mga"])
                writer.writerow(["TEST", *([str(source)] * 4), ragtag, hifiasm, lja, mga])
            config = root / "config.json"
            config.write_text(json.dumps({
                "project_dir": str(root / "project"),
                "sample_sheet": str(sheet),
                "oatk_db": str(source),
                "run_tidk": optional_qc,
                "run_heterozygosity": optional_qc,
                "run_qc_summary": optional_qc,
            }))
            result = subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", str(SNAKEFILE),
                 "--configfile", str(config), "--dry-run", "--cores", "4", "--nocolor"],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=90,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            scheduled = re.findall(r"^(?:local)?rule ([A-Za-z0-9_]+):", result.stdout, re.MULTILINE)
            rules = set(scheduled)
            self.assertEqual("hifiasm" in rules, hifiasm == "Y", result.stdout)
            self.assertEqual("lja" in rules, lja == "Y", result.stdout)
            self.assertEqual("mga" in rules, mga == "Y", result.stdout)
            self.assertEqual("install_mga" in rules, mga == "Y", result.stdout)
            self.assertEqual("ragtag_purged" in rules, ragtag == "Y" and hifiasm == "Y", result.stdout)
            self.assertEqual("ragtag_purged_lja" in rules, ragtag == "Y" and lja == "Y", result.stdout)
            self.assertEqual("barrnap" in rules, hifiasm == "Y", result.stdout)
            for name in ("tidk", "heterozygosity"):
                self.assertEqual(scheduled.count(name), (int(hifiasm == "Y") + int(lja == "Y")) if optional_qc else 0, result.stdout)
            self.assertEqual("qc_summary" in rules, optional_qc, result.stdout)
            for name in ("quast", "nquire_ploidy", "merqury", "all"):
                self.assertIn(name, rules, result.stdout)

    def test_both_without_ragtag(self):
        self.check_selection("Y", "Y", "N")

    def test_both_with_ragtag(self):
        self.check_selection("Y", "Y", "Y")

    def test_hifiasm_only_without_ragtag(self):
        self.check_selection("Y", "N", "N")

    def test_hifiasm_only_with_ragtag(self):
        self.check_selection("Y", "N", "Y")

    def test_lja_only_without_ragtag(self):
        self.check_selection("N", "Y", "N")

    def test_lja_only_with_ragtag(self):
        self.check_selection("N", "Y", "Y")

    def test_optional_mga_first_phase(self):
        self.check_selection("Y", "N", "N", mga="Y")

    def test_optional_qc_with_both_assemblers(self):
        self.check_selection("Y", "Y", "Y", optional_qc=True)

    def test_optional_qc_hifiasm_only(self):
        self.check_selection("Y", "N", "N", optional_qc=True)

    def test_optional_qc_lja_only(self):
        self.check_selection("N", "Y", "Y", optional_qc=True)

    def test_retired_legacy_mitochondrial_filters_are_not_registered(self):
        with tempfile.TemporaryDirectory(prefix="hifi-retired-filter-test-") as tmp:
            root = Path(tmp)
            source = root / "input.fa"
            source.write_text(">fixture\nACGTACGT\n")
            sheet = root / "samples.tsv"
            with sheet.open("w") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "HiFi", "Ref", "fq1", "fq2"])
                writer.writerow(["TEST", *([str(source)] * 4)])
            config = root / "config.json"
            config.write_text(json.dumps({
                "project_dir": str(root / "project"),
                "sample_sheet": str(sheet),
                "oatk_db": str(source),
            }))
            for retired in ("filter_hifiasm", "filter_lja"):
                with self.subTest(rule=retired):
                    result = subprocess.run(
                        [sys.executable, "-m", "snakemake", retired,
                         "--snakefile", str(SNAKEFILE), "--configfile", str(config),
                         "--dry-run", "--cores", "4", "--nocolor"],
                        cwd=root, text=True, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, timeout=90,
                    )
                    self.assertNotEqual(result.returncode, 0, result.stdout)
                    self.assertIn("MissingRuleException", result.stdout)
                    self.assertIn(f"No rule to produce {retired}", result.stdout)


if __name__ == "__main__":
    unittest.main()
