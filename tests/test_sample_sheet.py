"""Exercise sample-table values through isolated complete workflow dry-runs."""

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
REQUIRED_COLUMNS = ["sample", "HiFi", "Ref", "fq1", "fq2"]


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake runtime required")
class SampleSheetDagTests(unittest.TestCase):
    def run_dag(self, rows, optional_columns=(), pad_required=False):
        with tempfile.TemporaryDirectory(prefix="hifi-sample-sheet-test-") as temporary:
            root = Path(temporary)
            source = root / "input.fa"
            source.write_text(">fixture\nACGTACGT\n")
            sheet = root / "samples.tsv"
            with sheet.open("w") as handle:
                writer = csv.DictWriter(handle, REQUIRED_COLUMNS + list(optional_columns), delimiter="\t")
                writer.writeheader()
                for values in rows:
                    row = {column: str(source) for column in REQUIRED_COLUMNS if column != "sample"}
                    row.update({column: "" for column in optional_columns})
                    row.update(values)
                    if pad_required:
                        row = {key: f"  {value}  " if key in REQUIRED_COLUMNS else value for key, value in row.items()}
                    writer.writerow(row)
            config = root / "config.json"
            config.write_text(json.dumps({
                "project_dir": str(root / "project"),
                "sample_sheet": str(sheet),
                "oatk_db": str(source),
            }))
            return subprocess.run(
                [sys.executable, "-m", "snakemake", "--snakefile", str(SNAKEFILE),
                 "--configfile", str(config), "--dry-run", "--cores", "4", "--nocolor"],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=90,
            )

    def scheduled_samples(self, result, selected_rule):
        self.assertEqual(result.returncode, 0, result.stdout)
        current_rule = None
        selected = []
        for line in result.stdout.splitlines():
            match = re.match(r"^(?:local)?rule ([A-Za-z0-9_]+):", line)
            if match:
                current_rule = match[1]
            if current_rule == selected_rule and line.strip().startswith("wildcards:"):
                match = re.search(r"(?:wildcards: |, )sample=([^,]+)", line)
                if match:
                    selected.append(match[1])
        return sorted(selected)

    def test_header_only_table_fails_instead_of_an_empty_successful_dag(self):
        result = self.run_dag([])
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("Sample sheet must contain at least one sample row", result.stdout)

    def test_zero_padded_numeric_ids_remain_distinct_from_unpadded_ids(self):
        result = self.run_dag([{"sample": "001"}, {"sample": "1"}])
        self.assertEqual(self.scheduled_samples(result, "hifiasm"), ["001", "1"], result.stdout)
        self.assertEqual(self.scheduled_samples(result, "lja"), ["001", "1"], result.stdout)
        for sample in ("001", "1"):
            self.assertIn(f"/{sample}_FINAL_ASSEMBLIES/hifiasm/{sample}.purged.final.fasta", result.stdout)

    def test_literal_na_and_null_are_valid_sample_ids(self):
        result = self.run_dag([{"sample": "NA"}, {"sample": "NULL"}])
        self.assertEqual(self.scheduled_samples(result, "hifiasm"), ["NA", "NULL"], result.stdout)
        self.assertEqual(self.scheduled_samples(result, "lja"), ["NA", "NULL"], result.stdout)

    def test_numeric_one_flags_with_blank_cells_select_the_intended_jobs(self):
        result = self.run_dag([
            {"sample": "FIRST", "run_hifiasm": "1", "run_lja": "", "ragtag": "1", "run_mga": "1"},
            {"sample": "SECOND", "run_hifiasm": "", "run_lja": "1", "ragtag": "", "run_mga": ""},
        ], ["run_hifiasm", "run_lja", "ragtag", "run_mga"])
        for rule, expected in [("hifiasm", ["FIRST"]), ("lja", ["SECOND"]),
                               ("ragtag_purged", ["FIRST"]), ("ragtag_purged_lja", []),
                               ("mga", ["FIRST"])]:
            self.assertEqual(self.scheduled_samples(result, rule), expected, result.stdout)

    def test_blank_sample_ids_are_rejected(self):
        for value in ("", "   "):
            with self.subTest(sample=value):
                result = self.run_dag([{"sample": value}])
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn("Sample IDs must not be blank", result.stdout)

    def test_blank_required_input_paths_are_rejected(self):
        for column, value in [("HiFi", ""), ("Ref", "   "), ("fq1", ""), ("fq2", "   ")]:
            with self.subTest(column=column):
                result = self.run_dag([{"sample": "TEST", column: value}])
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(f"column {column!r} contains blank paths", result.stdout)

    def test_surrounding_whitespace_is_trimmed_from_required_values(self):
        result = self.run_dag([{"sample": "PADDED", "run_hifiasm": " 1 ", "run_lja": "0", "ragtag": " yes "}],
                              ["run_hifiasm", "run_lja", "ragtag"], pad_required=True)
        self.assertEqual(self.scheduled_samples(result, "hifiasm"), ["PADDED"], result.stdout)
        self.assertEqual(self.scheduled_samples(result, "lja"), [], result.stdout)
        self.assertEqual(self.scheduled_samples(result, "ragtag_purged"), ["PADDED"], result.stdout)


if __name__ == "__main__":
    unittest.main()
