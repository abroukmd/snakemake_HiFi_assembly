"""Regression coverage for the callable-site diploid-model estimator."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("estimate_heterozygosity", ROOT / "workflow/scripts/estimate_heterozygosity.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


class CallableLogicTests(unittest.TestCase):
    def test_explicit_organelle_ids_do_not_hide_nuclear_name_matches(self):
        reference = {"nuclear_mito_homolog": "ACGTNACGT", "S_mito": "AACCGG"}
        excluded = module.explicit_exclusions(reference, {"ctg1": "AACCGG"}, "S", [])
        with tempfile.TemporaryDirectory() as temporary:
            intervals, total = module.nuclear_intervals(reference, excluded, Path(temporary) / "nuclear.bed")
            self.assertEqual(total, 8)
            self.assertEqual(intervals["nuclear_mito_homolog"], [(0, 4), (5, 9)])
            self.assertNotIn("S_mito", intervals)
        with self.assertRaisesRegex(ValueError, "does not match"):
            module.explicit_exclusions(reference, {"ctg1": "TTTTTT"}, "S", [])

    def test_callable_mask_counts_exact_intervals_including_depth_bounds(self):
        rows = [("chr", 1, 10), ("chr", 2, 11), ("chr", 3, 0),
                ("chr", 4, 40), ("chr", 5, 41), ("chr", 6, 9)]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "callable.bed"
            self.assertEqual(module.callable_mask(rows, 10, 40, path), 3)
            mask, total = module.mask_from_bed(path)
            self.assertEqual(total, 3)
            index = module.mask_index(mask)
            self.assertTrue(module.in_mask(index, "chr", 1))
            self.assertFalse(module.in_mask(index, "chr", 3))
            self.assertTrue(module.in_mask(index, "chr", 4))
            self.assertFalse(module.in_mask(index, "chr", 5))

    def test_bam_header_rejects_wrong_reference_and_documents_identity_limit(self):
        reference = {"chr": "ACGT"}
        identity = module.check_bam_header("@SQ\tSN:chr\tLN:4\n", reference)
        self.assertTrue(identity["SN_LN_and_order_match"])
        self.assertIn("do not prove", identity["identity_boundary"])
        with self.assertRaisesRegex(ValueError, "names/lengths/order"):
            module.check_bam_header("@SQ\tSN:chr\tLN:5\n", reference)
        with self.assertRaisesRegex(ValueError, "M5 differs"):
            module.check_bam_header("@SQ\tSN:chr\tLN:4\tM5:wrong\n", reference)


@unittest.skipUnless(shutil.which("bcftools") and shutil.which("samtools"),
                     "Real caller fixture needs bcftools and samtools on PATH")
class CallerFixtureTests(unittest.TestCase):
    def fixture(self, folder, low_depth=False):
        rng = random.Random(20261022)
        sequence = "".join(rng.choice("ACGT") for _ in range(300))
        reference = folder / "input.fasta"
        reference.write_text(f">chr1\n{sequence}\n>S_mito\n{sequence}\n")
        mito = folder / "oatk.fasta"
        mito.write_text(f">oatk_ctg\n{sequence}\n")
        sam = folder / "reads.sam"
        with sam.open("w") as out:
            out.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n@SQ\tSN:S_mito\tLN:300\n@RG\tID:test\tSM:test\n")
            for contig in ["chr1", "S_mito"]:
                for j in range(2 if low_depth else 20):
                    bases = list(sequence[80:280])
                    alternatives = [base for base in "ACGT" if base != sequence[119]]
                    if j >= 10:
                        bases[39] = alternatives[0]  # True biallelic SNP at120.
                    alternatives = [base for base in "ACGT" if base != sequence[149]]
                    if j < 7:
                        bases[69] = alternatives[0]
                    elif j < 14:
                        bases[69] = alternatives[1]  # Original triallelic site at150.
                    qualities = list("I" * 200)
                    qualities[49] = "!"  # No BQ-passing evidence at130.
                    cigar = "200M"
                    if j >= 10:
                        bases.insert(100, "A")
                        qualities.insert(100, "I")
                        cigar = "100M1I100M"  # Insertion at180 must not duplicate DP positions.
                    out.write(f"{contig}_{j}\t0\t{contig}\t81\t60\t{cigar}\t*\t0\t0\t{''.join(bases)}\t{''.join(qualities)}\tRG:Z:test\n")
        bam = folder / "reads.bam"
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True, capture_output=True)
        subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
        return argparse.Namespace(reference=reference, bam=bam, bam_index=Path(str(bam) + ".bai"),
                                  mito_source=mito, sample="S", assembler="hifiasm", outdir=folder / "outputs",
                                  extra_excluded_ids=[], min_callable_fraction=.5, threads=4,
                                  bcftools="bcftools", samtools="samtools")

    def test_true_snp_multiallelic_indel_low_bq_and_organelle(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            args = self.fixture(folder)
            args.min_callable_fraction = .9
            before = {str(path): module.sha256(path) for path in [args.reference, args.bam, args.bam_index, args.mito_source]}
            metrics = module.estimate(args)
            self.assertEqual(metrics["status"], "low_callable_fraction")
            self.assertEqual(metrics["heterozygous_snps"], 1)
            self.assertEqual(metrics["callable_bases"], 199)
            self.assertAlmostEqual(metrics["heterozygosity_percent"], 100 / 199)
            with gzip.open(args.outdir / "heterozygous_snps.tsv.gz", "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(row["Contig"], row["Position1"]) for row in rows], [("chr1", "120")])
            mask, count = module.mask_from_bed(args.outdir / "callable.bed")
            self.assertEqual(count, 199)
            self.assertNotIn("S_mito", mask)
            self.assertFalse(module.in_mask(module.mask_index(mask), "chr1", 130))
            self.assertFalse(Path(str(args.reference) + ".fai").exists())
            self.assertTrue((args.outdir / "heterozygosity.done").is_file())
            receipt = json.loads((args.outdir / "receipt.json").read_text())
            self.assertEqual(receipt["metrics"]["genotype_model"], module.MODEL)
            self.assertEqual(receipt["exclusions"]["known_appended_OATK_ids"], ["S_mito"])
            for path, digest in before.items():
                self.assertEqual(module.sha256(path), digest)
            for name, record in receipt["outputs"].items():
                self.assertEqual(module.sha256(args.outdir / name), record["sha256"])

    def test_zero_callable_sites_are_untestable_not_zero_heterozygosity(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.fixture(Path(temporary), low_depth=True)
            args.threads = 1  # Sequential pileup and call honor a single-core allocation.
            metrics = module.estimate(args)
            self.assertEqual(metrics["status"], "untestable_no_callable_bases")
            self.assertEqual(metrics["callable_bases"], 0)
            self.assertIsNone(metrics["heterozygosity_percent"])
            self.assertEqual((args.outdir / "callable.bed").stat().st_size, 0)
            self.assertTrue((args.outdir / "heterozygosity.done").is_file())


    def test_explicit_all_excluded_reference_is_untestable_with_valid_empty_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.fixture(Path(temporary))
            args.extra_excluded_ids = ["chr1"]
            metrics = module.estimate(args)
            self.assertEqual(metrics["status"], "untestable_no_nuclear_acgt")
            self.assertEqual(metrics["nuclear_acgt_bases"], 0)
            self.assertIsNone(metrics["callable_fraction"])
            self.assertIsNone(metrics["heterozygosity_percent"])
            self.assertTrue((args.outdir / "heterozygous_snps.vcf.gz.tbi").is_file())

    def test_failed_reference_check_keeps_diagnostics_and_has_no_success_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.fixture(Path(temporary))
            text = args.reference.read_text()
            args.reference.write_text(text.replace(">chr1\n", ">chr1\nN", 1))
            with self.assertRaisesRegex(ValueError, "names/lengths/order"):
                module.estimate(args)
            status = json.loads((args.outdir / "STATUS.json").read_text())
            self.assertEqual(status["status"], "failed_incomplete")
            self.assertTrue(Path(status["scratch"]).is_dir())
            self.assertTrue((args.outdir / "commands.jsonl").is_file())
            self.assertFalse((args.outdir / "heterozygosity.done").exists())

if __name__ == "__main__":
    unittest.main()
