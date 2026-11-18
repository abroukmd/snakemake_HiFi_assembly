#!/usr/bin/env python3
"""Callable-site SNP heterozygosity under an explicit, unconfirmed diploid model.

Portable adaptation of the corrected EPIC BAM estimator (22 October 2026).
The same all-site bcftools caller FORMAT/DP supplies denominator and SNP filters.
Originally multiallelic records are never decomposed. Exclusions are explicit
IDs, backed by exact OATK-to-final sequence matches, rather than name regexes.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import tempfile
import traceback

PARAMETERS = {"min_mq": 30, "min_bq": 20, "min_dp": 10, "max_dp_multiplier": 2.0,
              "min_qual": 30, "min_ref_ad": 3, "min_alt_ad": 3, "ploidy": 2,
              "keep_alts": True, "pileup_max_depth": 10000,
              "excluded_read_flags": "UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP",
              "depth_mean": "positive_nuclear_ACGT_all_site_caller_FORMAT_DP",
              "BAQ": "bcftools_default", "overlap_handling": "bcftools_default",
              "anomalous_pairs": "excluded_by_bcftools_default", "GQ_filter": None,
              "duplicate_marking": "none; preflagged duplicate reads excluded"}
NON_NUCLEOTIDE_TYPES = 'TYPE="indel" || TYPE="mnp" || TYPE="other" || TYPE="bnd"'
SNP_FIELDS = ["Contig", "Position1", "REF", "ALT", "QUAL", "GT", "DP", "REF_AD", "ALT_AD", "ALT_fraction", "GQ"]
MODEL = "STANDARDIZED_DIPLOID_MODEL_NOT_CONFIRMED_BIOLOGICAL_PLOIDY"


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta(path):
    records, identifier, sequence = {}, None, []
    def append():
        if identifier is not None:
            value = "".join(sequence).upper()
            require(value and identifier not in records, f"Empty sequence or duplicate FASTA ID {identifier}: {path}")
            require(not (set(value) - set("ACGTURYSWKMBDHVN")), f"Invalid nucleotide sequence {identifier}: {path}")
            records[identifier] = value
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                append()
                require(line[1:].strip(), f"Empty FASTA header: {path}")
                identifier, sequence = line[1:].split()[0], []
            else:
                require(identifier is not None, f"Sequence before header: {path}")
                sequence.append("".join(line.split()))
    append()
    require(records, f"No FASTA records: {path}")
    return records


def appended_mito_ids(sample, count):
    return ([f"{sample}_mito"] if count == 1 else
            [f"{sample}_mito_{index}" for index in range(1, count + 1)])


def explicit_exclusions(reference, mitochondria, sample, extra):
    require(isinstance(extra, list) and all(isinstance(x, str) and x for x in extra),
            "Extra excluded IDs must be a JSON list of nonempty strings")
    added = appended_mito_ids(sample, len(mitochondria))
    for identifier, sequence in zip(added, mitochondria.values()):
        require(reference.get(identifier) == sequence,
                f"Known appended mitochondrial ID/sequence does not match final FASTA: {identifier}")
    excluded = set(added) | set(extra)
    require(excluded <= set(reference), f"Excluded IDs absent from reference: {sorted(excluded - set(reference))}")
    return sorted(excluded)


def nuclear_intervals(reference, excluded, output):
    intervals, total = {}, 0
    with Path(output).open("w") as handle:
        for identifier, sequence in reference.items():
            if identifier in excluded:
                continue
            intervals[identifier] = [(match.start(), match.end()) for match in re.finditer("[ACGT]+", sequence)]
            for start, end in intervals[identifier]:
                handle.write(f"{identifier}\t{start}\t{end}\n")
                total += end - start
    return intervals, total


def mask_from_bed(path):
    intervals = defaultdict(list)
    total = 0
    with Path(path).open() as handle:
        for line in handle:
            chrom, start, end = line.split()
            start, end = int(start), int(end)
            require(0 <= start < end, f"Invalid BED interval: {line.strip()}")
            require(not intervals[chrom] or start >= intervals[chrom][-1][1], "Overlapping or unsorted BED intervals")
            intervals[chrom].append((start, end))
            total += end - start
    return intervals, total


def mask_index(intervals):
    return {chrom: ([start for start, _ in records], [end for _, end in records])
            for chrom, records in intervals.items()}


def in_mask(index, chrom, position1):
    if chrom not in index:
        return False
    starts, ends = index[chrom]
    item = bisect_right(starts, position1 - 1) - 1
    return item >= 0 and position1 <= ends[item]


class Tools:
    def __init__(self, folder, bcftools, samtools):
        self.folder, self.bcftools, self.samtools = folder, bcftools, samtools

    def record(self, **command):
        with (self.folder / "commands.jsonl").open("a") as handle:
            handle.write(json.dumps({"time_utc": utcnow(), **command}) + "\n")

    def run(self, argv, stdout=None):
        self.record(argv=argv)
        with (self.folder / "tools.log").open("ab") as log:
            subprocess.run(argv, stdout=stdout if stdout is not None else log, stderr=log, check=True)

    def capture(self, argv):
        self.record(argv=argv)
        with (self.folder / "tools.log").open("ab") as log:
            return subprocess.check_output(argv, stderr=log, text=True)

    def stream(self, argv):
        self.record(argv=argv)
        with (self.folder / "tools.log").open("ab") as log:
            process = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=log, text=True)
            try:
                for line in process.stdout:
                    yield line.rstrip("\n")
            finally:
                process.stdout.close()
                code = process.wait()
                if code:
                    raise subprocess.CalledProcessError(code, argv)

    def pipeline(self, first, second):
        self.record(pipeline=[first, second])
        with (self.folder / "tools.log").open("ab") as log:
            producer = subprocess.Popen(first, stdout=subprocess.PIPE, stderr=log)
            consumer = subprocess.Popen(second, stdin=producer.stdout, stdout=log, stderr=log)
            producer.stdout.close()
            codes = [consumer.wait(), producer.wait()]
            require(codes == [0, 0], f"Pileup/caller pipeline failed: exit codes {codes}")


def check_bam_header(header, reference):
    sq, samples, md5_checked = [], set(), []
    for line in header.splitlines():
        fields = line.split("\t")
        tags = dict(field.split(":", 1) for field in fields[1:] if ":" in field)
        if fields[0] == "@SQ":
            require("SN" in tags and "LN" in tags, "BAM @SQ is missing SN/LN")
            sq.append((tags["SN"], int(tags["LN"])))
            if "M5" in tags:
                require(tags["SN"] in reference and
                        hashlib.md5(reference[tags["SN"]].encode()).hexdigest().lower() == tags["M5"].lower(),
                        f"BAM M5 differs from final reference for {tags['SN']}")
                md5_checked.append(tags["SN"])
        if fields[0] == "@RG" and "SM" in tags:
            samples.add(tags["SM"])
    require(sq == [(name, len(sequence)) for name, sequence in reference.items()],
            "BAM reference names/lengths/order differ from exact final FASTA")
    require(len(samples) <= 1, f"Multiple BAM samples found: {sorted(samples)}")
    return {"SN_LN_and_order_match": True, "M5_checked_contigs": md5_checked,
            "read_group_samples": sorted(samples),
            "identity_boundary": "SN/LN do not prove sequence identity where M5 is absent; exact final accessor association is inherited from the mapping DAG"}


def depth_rows(tools, calls, nuclear_index):
    previous = {}
    command = [tools.bcftools, "query", "-e", NON_NUCLEOTIDE_TYPES, "-f", "%CHROM\t%POS\t%REF[\t%DP]\n", str(calls)]
    for line in tools.stream(command):
        fields = line.split("\t")
        require(len(fields) == 4, f"Expected one called sample, found {len(fields) - 3}")
        chrom, position, refbase, depth = fields
        position = int(position)
        require(position > previous.get(chrom, 0), f"Repeated or unsorted nucleotide position: {chrom}:{position}")
        previous[chrom] = position
        require(len(refbase) == 1 and refbase in "ACGT" and in_mask(nuclear_index, chrom, position),
                f"Non-nuclear/non-ACGT position in callable evidence: {chrom}:{position}")
        if depth != ".":
            depth = int(depth)
            require(depth >= 0, "Negative FORMAT/DP")
            yield chrom, position, depth


def callable_mask(rows, minimum, maximum, output):
    count, current = 0, None
    with Path(output).open("w") as handle:
        for chrom, position, depth in rows:
            if minimum <= depth <= maximum:
                count += 1
                if current and chrom == current[0] and position - 1 == current[2]:
                    current = (chrom, current[1], position)
                else:
                    if current:
                        handle.write("\t".join(map(str, current)) + "\n")
                    current = (chrom, position - 1, position)
            elif current:
                handle.write("\t".join(map(str, current)) + "\n")
                current = None
        if current:
            handle.write("\t".join(map(str, current)) + "\n")
    return count


def extract_snps(tools, passing, callable_index, maximum, output, balances_output):
    seen, balances, bins = set(), [], Counter()
    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SNP_FIELDS, delimiter="\t")
        writer.writeheader()
        command = [tools.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AD\t%GQ]\n", str(passing)]
        for line in tools.stream(command):
            fields = line.split("\t")
            require(len(fields) == 9, "Expected exactly one genotype sample in SNP query")
            chrom, position, refbase, alt, qual, gt, depth, ad, gq = fields
            position, depth, quality = int(position), int(depth), float(qual)
            alleles = ad.split(",")
            require(len(alleles) == 2, "Originally multiallelic site reached the biallelic SNP set")
            ref_ad, alt_ad = map(int, alleles)
            key = (chrom, position)
            require(key not in seen and in_mask(callable_index, chrom, position),
                    f"Duplicate SNP or SNP outside denominator mask: {key}")
            require(gt in {"0/1", "1/0", "0|1", "1|0"} and len(refbase) == len(alt) == 1 and
                    refbase in "ACGT" and alt in "ACGT" and alt != refbase, "Non-biallelic heterozygous SNP")
            require(math.isfinite(quality) and quality >= PARAMETERS["min_qual"] and
                    PARAMETERS["min_dp"] <= depth <= maximum and
                    ref_ad >= PARAMETERS["min_ref_ad"] and alt_ad >= PARAMETERS["min_alt_ad"], "SNP filter mismatch")
            seen.add(key)
            balance = alt_ad / (ref_ad + alt_ad)
            balances.append(balance)
            bins[min(19, int(balance * 20))] += 1
            writer.writerow(dict(zip(SNP_FIELDS, [chrom, position, refbase, alt, qual, gt,
                                                  depth, ref_ad, alt_ad, balance, gq])))
    with Path(balances_output).open("w") as handle:
        handle.write("ALT_fraction_start\tALT_fraction_end\tSNP_sites\n")
        for index in range(20):
            handle.write(f"{index / 20:.2f}\t{(index + 1) / 20:.2f}\t{bins[index]}\n")
    return len(seen), statistics.median(balances) if balances else None


def estimate(args):
    folder = args.outdir.absolute()
    folder.mkdir(parents=True, exist_ok=True)
    done = folder / "heterozygosity.done"
    done.unlink(missing_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix="work.", dir=folder))
    started = utcnow()
    save_json(folder / "STATUS.json", {"status": "running", "started_utc": started,
                                       "pid": os.getpid(), "scratch": str(scratch)})
    try:
        reference_path = args.reference.resolve(strict=True)
        bam = args.bam.resolve(strict=True)
        bam_index = args.bam_index.resolve(strict=True)
        mito_path = args.mito_source.resolve(strict=True)
        before_stats = {str(path): [path.stat().st_size, path.stat().st_mtime_ns] for path in [reference_path, bam, bam_index, mito_path]}
        tools = Tools(folder, args.bcftools, args.samtools)
        versions = {name: tools.capture([binary, "--version"]).splitlines()[0]
                    for name, binary in [("bcftools", args.bcftools), ("samtools", args.samtools)]}
        original_reference = fasta(reference_path)
        mitochondrial = fasta(mito_path)
        excluded = explicit_exclusions(original_reference, mitochondrial, args.sample, args.extra_excluded_ids)
        identity = check_bam_header(tools.capture([args.samtools, "view", "-H", str(bam)]), original_reference)
        tools.run([args.samtools, "quickcheck", "-v", str(bam)])
        inputs = {"reference": {"path": str(reference_path), "sha256": sha256(reference_path)},
                  "bam": {"path": str(bam), "sha256": sha256(bam)},
                  "bam_index": {"path": str(bam_index), "sha256": sha256(bam_index)},
                  "mitochondrial_source": {"path": str(mito_path), "sha256": sha256(mito_path)}}
        private_reference = scratch / "reference.fasta"
        shutil.copyfile(reference_path, private_reference)
        require(sha256(private_reference) == inputs["reference"]["sha256"], "Private reference copy hash mismatch")
        tools.run([args.samtools, "faidx", str(private_reference)])
        # A private BAM link and explicit index link avoid source-adjacent writes.
        private_bam = scratch / "reads.bam"
        private_bam.symlink_to(bam)
        (scratch / "reads.bam.bai").symlink_to(bam_index)
        nuclear_bed = scratch / "nuclear_acgt.bed"
        intervals, nuclear_total = nuclear_intervals(original_reference, excluded, nuclear_bed)
        nuclear_index = mask_index(intervals)
        (scratch / "excluded_contigs.txt").write_text("".join(identifier + "\n" for identifier in excluded))
        calls = scratch / "all_sites.bcf"
        compression_threads = max(0, args.threads - 2)
        if nuclear_total:
            mpileup = [args.bcftools, "mpileup", "--threads", "0", "--fasta-ref", str(private_reference),
                       "--regions-file", str(nuclear_bed), "--annotate", "FORMAT/AD,FORMAT/DP",
                       "--min-MQ", str(PARAMETERS["min_mq"]), "--min-BQ", str(PARAMETERS["min_bq"]),
                       "--max-depth", str(PARAMETERS["pileup_max_depth"]), "--skip-any-set", PARAMETERS["excluded_read_flags"],
                       "--seed", "0", "--output-type", "u", str(private_bam)]
            caller = [args.bcftools, "call", "--threads", str(compression_threads), "--multiallelic-caller", "--keep-alts",
                      "--ploidy", "2", "--annotate", "FORMAT/GQ", "--output-type", "b", "--output", str(calls)]
            if args.threads == 1:
                # Honor a single-core allocation by avoiding concurrent processes.
                pileup_file = scratch / "pileup.bcf"
                tools.run(mpileup[:-1] + ["--output", str(pileup_file), mpileup[-1]])
                tools.run(caller + [str(pileup_file)])
            else:
                tools.pipeline(mpileup, caller)
        else:
            empty = scratch / "empty.vcf"
            with empty.open("w") as handle:
                handle.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                for identifier, sequence in original_reference.items():
                    handle.write(f"##contig=<ID={identifier},length={len(sequence)}>\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            tools.run([args.bcftools, "view", "-Ob", "-o", str(calls), str(empty)])
        histogram, visited = Counter(), 0
        if nuclear_total:
            for _, _, depth in depth_rows(tools, calls, nuclear_index):
                visited += 1
                histogram[depth] += 1
        require(visited <= nuclear_total, "Nucleotide evidence contains duplicate genomic positions")
        covered = sum(count for depth, count in histogram.items() if depth > 0)
        positive_mean = sum(depth * count for depth, count in histogram.items() if depth > 0) / covered if covered else 0
        maximum = max(PARAMETERS["min_dp"], int(positive_mean * PARAMETERS["max_dp_multiplier"] + .5))
        mask_path = scratch / "callable.bed"
        callable_bases = callable_mask(depth_rows(tools, calls, nuclear_index) if nuclear_total else [], PARAMETERS["min_dp"], maximum, mask_path)
        expected = sum(count for depth, count in histogram.items() if PARAMETERS["min_dp"] <= depth <= maximum)
        mask, mask_bases = mask_from_bed(mask_path)
        require(callable_bases == expected == mask_bases, "Callable BED and depth histogram disagree")
        with (scratch / "depth_histogram.tsv").open("w") as handle:
            handle.write("DP\tPositions\n")
            for depth, count in sorted(histogram.items()):
                handle.write(f"{depth}\t{count}\n")
        raw_snps, passing = scratch / "biallelic_snps.vcf.gz", scratch / "heterozygous_snps.vcf.gz"
        if callable_bases:
            tools.run([args.bcftools, "view", "--threads", str(max(0, args.threads - 1)), "--targets-file", str(mask_path),
                       "--targets-overlap", "0", "--types", "snps", "--min-alleles", "2", "--max-alleles", "2",
                       "--include", f'QUAL>={PARAMETERS["min_qual"]}', "-Oz", "-o", str(raw_snps), str(calls)])
            expression = (f'GT="het" && FMT/DP>={PARAMETERS["min_dp"]} && FMT/DP<={maximum} && '
                          f'FMT/AD[0:0]>={PARAMETERS["min_ref_ad"]} && FMT/AD[0:1]>={PARAMETERS["min_alt_ad"]}')
            tools.run([args.bcftools, "view", "--threads", str(max(0, args.threads - 1)), "--include", expression,
                       "-Oz", "-o", str(passing), str(raw_snps)])
        else:
            tools.run([args.bcftools, "view", "--header-only", "-Oz", "-o", str(raw_snps), str(calls)])
            shutil.copyfile(raw_snps, passing)
        for variants in [raw_snps, passing]:
            tools.run([args.bcftools, "index", "--tbi", "--force", str(variants)])
        het_sites, median_balance = extract_snps(tools, passing, mask_index(mask), maximum,
            scratch / "heterozygous_snps.tsv.gz", scratch / "allele_balance.tsv") if callable_bases else empty_snps(scratch)
        fraction = callable_bases / nuclear_total if nuclear_total else None
        status = ("untestable_no_nuclear_acgt" if not nuclear_total else "untestable_no_callable_bases" if not callable_bases
                  else "low_callable_fraction" if fraction < args.min_callable_fraction else "complete")
        metrics = {"sample": args.sample, "assembler": args.assembler, "status": status, "genotype_model": MODEL,
                   "nuclear_acgt_bases": nuclear_total, "positions_with_nucleotide_evidence": visited,
                   "mean_positive_nuclear_evidence_dp": positive_mean, "minimum_callable_dp": PARAMETERS["min_dp"],
                   "maximum_callable_dp": maximum, "callable_bases": callable_bases, "callable_fraction": fraction,
                   "heterozygous_snps": het_sites, "heterozygosity_proportion": het_sites / callable_bases if callable_bases else None,
                   "heterozygosity_percent": 100 * het_sites / callable_bases if callable_bases else None,
                   "het_snps_per_kb": 1000 * het_sites / callable_bases if callable_bases else None,
                   "median_alt_fraction": median_balance, "callability_review_threshold": args.min_callable_fraction,
                   "excluded_contig_ids": excluded, "pileup_depth_cap": PARAMETERS["pileup_max_depth"],
                   "positions_at_depth_cap": sum(count for depth, count in histogram.items() if depth >= PARAMETERS["pileup_max_depth"]),
                   "interpretation": "Conditional depth-callable diploid-model estimate; mapping, collapsed paralogs, haplotype separation, contamination, and unknown organellar contigs can affect it"}
        save_json(scratch / "summary.json", metrics)
        with (scratch / "summary.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(metrics), delimiter="\t")
            writer.writeheader()
            writer.writerow({key: json.dumps(value) if isinstance(value, list) else value for key, value in metrics.items()})
        for path in [reference_path, bam, bam_index, mito_path]:
            require([path.stat().st_size, path.stat().st_mtime_ns] == before_stats[str(path)], f"Source changed during analysis: {path}")
        outputs = ["summary.json", "summary.tsv", "callable.bed", "nuclear_acgt.bed", "excluded_contigs.txt",
                   "depth_histogram.tsv", "biallelic_snps.vcf.gz", "biallelic_snps.vcf.gz.tbi",
                   "heterozygous_snps.vcf.gz", "heterozygous_snps.vcf.gz.tbi", "heterozygous_snps.tsv.gz", "allele_balance.tsv"]
        for name in outputs:
            os.replace(scratch / name, folder / name)
        receipt = {"state": "completed", "started_utc": started, "finished_utc": utcnow(), "inputs": inputs,
                   "implementation": {"script": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
                   "parameters": {**PARAMETERS, "threads": args.threads, "callability_review_threshold": args.min_callable_fraction},
                   "reference_consistency": identity, "versions": versions, "metrics": metrics,
                   "exclusions": {"known_appended_OATK_ids": appended_mito_ids(args.sample, len(mitochondrial)), "additional_configured_ids": args.extra_excluded_ids,
                                  "limitation": "Unknown organellar sequences are not detected by broad name matching"},
                   "outputs": {name: {"sha256": sha256(folder / name), "bytes": (folder / name).stat().st_size} for name in outputs},
                   "private_reference": str(private_reference), "scratch": str(scratch),
                   "large_intermediate_policy": "all_sites.bcf deleted only after validated outputs and receipt; failed scratch is retained"}
        save_json(folder / "receipt.json", receipt)
        calls.unlink()
        (scratch / "pileup.bcf").unlink(missing_ok=True)
        save_json(folder / "STATUS.json", {"status": status, "finished_utc": utcnow(), "receipt": str(folder / "receipt.json")})
        done.touch()
        print(json.dumps(metrics, allow_nan=False), flush=True)
        return metrics
    except BaseException as error:
        save_json(folder / "STATUS.json", {"status": "failed_incomplete", "error": str(error),
            "traceback": traceback.format_exc(), "finished_utc": utcnow(), "scratch": str(scratch)})
        raise


def empty_snps(scratch):
    with gzip.open(scratch / "heterozygous_snps.tsv.gz", "wt") as handle:
        handle.write("\t".join(SNP_FIELDS) + "\n")
    with (scratch / "allele_balance.tsv").open("w") as handle:
        handle.write("ALT_fraction_start\tALT_fraction_end\tSNP_sites\n")
        for index in range(20):
            handle.write(f"{index / 20:.2f}\t{(index + 1) / 20:.2f}\t0\n")
    return 0, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["reference", "bam", "bam-index", "mito-source", "outdir"]:
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembler", choices=["hifiasm", "lja"], required=True)
    parser.add_argument("--extra-excluded-ids-json", default="[]")
    parser.add_argument("--min-callable-fraction", type=float, default=.5)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    args.extra_excluded_ids = json.loads(args.extra_excluded_ids_json)
    require(1 <= args.threads <= 4, "This module supports 1–4 allocated threads")
    require(0 <= args.min_callable_fraction <= 1, "Callability review threshold must be within [0,1]")
    estimate(args)


if __name__ == "__main__":
    main()
