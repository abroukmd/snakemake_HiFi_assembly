"""Exercise actual RagTag shell fragments for sequence retention and name maps."""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import textwrap
import unittest


RULES = Path(__file__).resolve().parents[1] / "workflow/rules/ragtag.smk"


def rule_shell(name):
    section = RULES.read_text().split(f"rule {name}:", 1)[1]
    return re.search(r"shell:\s+r'''(.*?)'''", section, re.DOTALL).group(1)


def render(fragment, substitutions):
    for name, value in substitutions.items():
        fragment = fragment.replace("{" + name + "}", str(value))
    return fragment.replace("{{", "{").replace("}}", "}")


class RagTagHelperTests(unittest.TestCase):
    def test_scaffold_conservation_guard_accepts_reordering_and_rejects_loss(self):
        for rule in ("ragtag_purged", "ragtag_purged_lja"):
            with self.subTest(rule=rule):
                shell = rule_shell(rule)
                guard = re.search(r"^        if ! cmp -s .*?^        fi$", shell,
                                  re.MULTILINE | re.DOTALL)
                self.assertIsNotNone(guard, "A completed reorder must verify all scaffold IDs and lengths")
                for case, ordered, succeeds in (
                    ("reordered", "ref2\t4\t0\t4\t5\nref1\t4\t0\t4\t5\n", True),
                    ("missing_scaffold", "ref1\t4\t0\t4\t5\n", False),
                    ("changed_length", "ref1\t4\t0\t4\t5\nref2\t3\t0\t3\t4\n", False),
                ):
                    with self.subTest(case=case), tempfile.TemporaryDirectory(prefix="ragtag-conservation-") as temporary:
                        root = Path(temporary)
                        native_fai, ordered_fai = root / "native.fai", root / "ordered.fai"
                        native_fai.write_text("ref1\t4\t0\t4\t5\nref2\t4\t0\t4\t5\n")
                        ordered_fai.write_text(ordered)
                        done = root / "complete.done"
                        command = render("set -euo pipefail\n" + guard.group(0), {
                            "output.fai": native_fai, "output.ordered_fai": ordered_fai,
                            "log.err": root / "stderr.log"})
                        command += "\ntouch " + shlex.quote(str(done)) + "\n"
                        result = subprocess.run(["bash", "-c", command], capture_output=True,
                                                text=True, timeout=10)
                        self.assertEqual(result.returncode == 0, succeeds, result.stdout + result.stderr)
                        self.assertEqual(done.exists(), succeeds)

    def test_failed_scaffold_extraction_aborts_before_completion(self):
        for rule in ("ragtag_purged", "ragtag_purged_lja"):
            with self.subTest(rule=rule), tempfile.TemporaryDirectory(prefix="ragtag-faidx-failure-") as temporary:
                root = Path(temporary)
                tools = root / "bin"
                tools.mkdir()
                samtools = tools / "samtools"
                samtools.write_text("#!/usr/bin/env python3\n" + '''import pathlib, sys
assert sys.argv[1] == "faidx"
fasta = pathlib.Path(sys.argv[2])
if len(sys.argv) == 3:
    fasta.with_name(fasta.name + ".fai").write_text("ref1_RagTag\\t4\\t0\\t4\\t5\\nref2_RagTag\\t4\\t0\\t4\\t5\\n")
elif sys.argv[3] == "ref1_RagTag":
    print(">ref1_RagTag\\nACGT")
else:
    print("simulated failure extracting second scaffold", file=sys.stderr)
    sys.exit(1)
''')
                samtools.chmod(0o755)
                reference, agp, fasta = root / "reference.fa", root / "scaffolds.agp", root / "scaffolds.fa"
                reference.write_text(">ref1\nACGT\n>ref2\nTGCA\n")
                agp.write_text("ref1_RagTag\t1\t4\t1\tW\tcontig1\t1\t4\t+\n"
                               "ref2_RagTag\t1\t4\t1\tW\tcontig2\t1\t4\t+\n")
                fasta.write_text(">ref1_RagTag\nACGT\n>ref2_RagTag\nTGCA\n")
                shell = rule_shell(rule)
                fragment = shell.split("        reorder_scaffolds_by_ref_keep_unmatched()", 1)[1]
                fragment = "reorder_scaffolds_by_ref_keep_unmatched()" + fragment.split(
                    "\n        rename_fasta_to_sample_scaffs()", 1)[0]
                fragment = render(textwrap.dedent(fragment),
                                  {"log.out": root / "stdout.log", "log.err": root / "stderr.log"})
                output, done = root / "ordered.fa", root / "complete.done"
                command = "set -euo pipefail\n" + fragment + "\n" + \
                    "reorder_scaffolds_by_ref_keep_unmatched " + " ".join(
                        shlex.quote(str(path)) for path in (reference, agp, fasta, output)) + \
                    "\ntouch " + shlex.quote(str(done)) + "\n"
                result = subprocess.run(["bash", "-c", command],
                                        env=dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"]),
                                        capture_output=True, text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertFalse(done.exists(), "A partial scaffold extraction must not be marked complete")
                self.assertIn("ref1_RagTag", output.read_text(), "First scaffold must succeed before the injected failure")

    def test_copied_haplotigs_have_identity_name_maps(self):
        for rule in ("ragtag_purged", "ragtag_purged_lja"):
            with self.subTest(rule=rule), tempfile.TemporaryDirectory(prefix="ragtag-hap-map-") as temporary:
                root = Path(temporary)
                original, copied, mapping = root / "hap.fa", root / "final_hap.fa", root / "hap.map.tsv"
                original.write_text(">hap_alpha description\nACGT\n>hap_beta\nTGCA\n")
                shell = rule_shell(rule)
                fragments = [line.strip() for line in shell.splitlines()
                             if (line.strip().startswith("cp -f ") and "{output.final_hap}" in line)
                             or (line.strip().startswith("awk ") and "{output.final_hap_map}" in line)]
                self.assertEqual(len(fragments), 2)
                command = render("set -euo pipefail\n" + "\n".join(fragments), {
                    "input.hap": original, "output.final_hap": copied,
                    "output.final_hap_map": mapping, "wildcards.sample": "TEST"})
                result = subprocess.run(["bash", "-c", command], capture_output=True, text=True, timeout=10)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(copied.read_bytes(), original.read_bytes())
                identifiers = [line[1:].split()[0] for line in copied.read_text().splitlines() if line.startswith(">")]
                self.assertEqual([line.split("\t") for line in mapping.read_text().splitlines()],
                                 [[identifier, identifier] for identifier in identifiers])


if __name__ == "__main__":
    unittest.main()
