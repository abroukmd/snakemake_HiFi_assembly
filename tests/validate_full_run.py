#!/usr/bin/env python3
"""Validate completed pipeline outputs without rerunning scientific analyses.

Usage: python tests/validate_full_run.py --config /run/config.json --output /run/validation.json
A partial run produces a JSON report and a nonzero exit status. Source reads and
assemblies are never modified. Only the explicitly requested JSON is written.
"""
import argparse
from bisect import bisect_right
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import statistics
import subprocess
import sys
import zlib
import xml.etree.ElementTree as ET

TRUTHY = {"1", "true", "t", "yes", "y"}
COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")


def truthy(value):
    return str(value).strip().lower() in TRUTHY


def selected_assemblers(row):
    selected = []
    for assembler, aliases in [("hifiasm", ["run_hifiasm", "hifiasm", "use_hifiasm"]),
                               ("lja", ["run_lja", "lja", "use_lja"])]:
        flag = next((name for name in aliases if name in row), None)
        if flag is None or truthy(row[flag]):
            selected.append(assembler)
    return selected


def require(condition, message):
    if not condition:
        raise ValueError(message)


def nonempty(path):
    path = Path(path)
    require(path.is_file() and path.stat().st_size > 0, f"Missing or empty file: {path}")
    return path


def marker(path, artifacts=()):
    path = Path(path)
    require(path.is_file(), f"Missing completion marker: {path}")
    for artifact in artifacts:
        artifact = Path(artifact)
        if artifact.name.endswith(".done"):
            require(artifact.is_file(), f"Missing prerequisite marker: {artifact}")
        else:
            artifact = nonempty(artifact)
        require(path.stat().st_mtime + 2 >= artifact.stat().st_mtime,
                f"Marker predates its artifact: {path} / {artifact}")


def fasta(path):
    path = nonempty(path)
    records, identifier, sequence = {}, None, []
    def append():
        if identifier is not None:
            value = "".join(sequence).upper()
            require(value, f"Empty FASTA record {identifier} in {path}")
            require(not (set(value) - set("ACGTURYSWKMBDHVN")),
                    f"Invalid nucleotide characters in {identifier}: {path}")
            require(identifier not in records, f"Duplicate FASTA identifier {identifier}: {path}")
            records[identifier] = value
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                append()
                require(line[1:].strip(), f"Empty FASTA header: {path}")
                identifier = line[1:].split()[0]
                sequence = []
            else:
                require(identifier is not None, f"Sequence before FASTA header: {path}")
                sequence.append("".join(line.split()))
    append()
    require(records, f"No FASTA records: {path}")
    return records


def fasta_stats(path, records, k):
    lengths = sorted(map(len, records.values()), reverse=True)
    total, cumulative, n50 = sum(lengths), 0, 0
    for length in lengths:
        cumulative += length
        if cumulative * 2 >= total:
            n50 = length
            break
    counts = Counter("".join(records.values()))
    unambiguous = sum(counts[base] for base in "ACGT")
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest(), "sequence_count": len(lengths),
            "total_bp": total, "n50_bp": n50, "longest_bp": lengths[0],
            "n_bp": counts["N"], "gc_percent_excluding_ambiguity":
            100 * (counts["G"] + counts["C"]) / unambiguous if unambiguous else None,
            "maximum_possible_kmer_windows": sum(max(0, length - k + 1) for length in lengths)}


def png_info(path):
    path = nonempty(path)
    with path.open("rb") as handle:
        require(handle.read(8) == b"\x89PNG\r\n\x1a\n", f"Invalid PNG signature: {path}")
        seen = set()
        dimensions = None
        while True:
            header = handle.read(8)
            require(len(header) == 8, f"Truncated PNG chunk header: {path}")
            length, kind = struct.unpack(">I4s", header)
            require(length <= path.stat().st_size, f"Invalid PNG chunk length: {path}")
            data, checksum = handle.read(length), handle.read(4)
            require(len(data) == length and len(checksum) == 4, f"Truncated PNG: {path}")
            require(zlib.crc32(kind + data) & 0xffffffff == struct.unpack(">I", checksum)[0],
                    f"PNG CRC mismatch: {path}")
            seen.add(kind)
            if kind == b"IHDR":
                require(length == 13 and dimensions is None, f"Invalid PNG IHDR: {path}")
                dimensions = struct.unpack(">II", data[:8])
                require(min(dimensions) > 0, f"Empty PNG dimensions: {path}")
            if kind == b"IEND":
                require(length == 0 and b"IHDR" in seen and b"IDAT" in seen,
                        f"Incomplete PNG: {path}")
                break
        require(not handle.read(1), f"Unexpected bytes after PNG IEND: {path}")
    return {"path": str(path), "width": dimensions[0], "height": dimensions[1],
            "size_bytes": path.stat().st_size}


def circular_equal(first, second):
    if len(first) != len(second):
        return False
    return first in second + second or first.translate(COMPLEMENT)[::-1] in second + second


def integration_check(project, sample, assembler, raw_final, mito, ragtag, config):
    assembly_root = project / "Assemblies" / sample
    if ragtag:
        branch = "RAGTAG" if assembler == "hifiasm" else "RAGTAG_LJA"
        before_path = assembly_root / branch / f"FINAL_{sample}" / f"{sample}.purged.final.fasta"
        done_name = f"{sample}_ragtag_purged.done" if assembler == "hifiasm" else f"{sample}_ragtag_lja.done"
        marker(assembly_root / branch / done_name, [before_path])
    else:
        branch = "PURGE_DUPS" if assembler == "hifiasm" else "PURGE_DUPS_LJA"
        before_path = assembly_root / branch / f"purged-{assembler}.sorted.fa"
    before = fasta(before_path)
    integration = assembly_root / "OATK/mito_integration" / assembler
    removed_path = integration / "contigs_to_remove.txt"
    require(removed_path.is_file(), f"Missing removed-ID list: {removed_path}")
    removed = set(removed_path.read_text().split())
    require(removed <= set(before), f"Removed IDs are absent from assembly input: {removed - set(before)}")
    require(not (removed & set(raw_final)), f"Removed IDs survive in final FASTA: {removed & set(raw_final)}")
    report_path = nonempty(integration / "integration.tsv")
    with report_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    recorded_removed = set()
    for row in rows:
        require(row["assembly_id"] in before and row["mito_id"] in mito,
                f"Integration report references unknown IDs: {report_path}")
        qcov, tcov = float(row["query_coverage"]), float(row["target_coverage"])
        require(all(math.isfinite(x) and 0 <= x <= 1 for x in (qcov, tcov)),
                f"Invalid integration coverage: {report_path}")
        require(row["decision"] in {"keep", "remove"}, f"Invalid integration decision: {report_path}")
        if row["decision"] == "remove":
            require(qcov + 1e-8 >= float(config.get("mito_min_query_coverage", .8)) and
                    tcov + 1e-8 >= float(config.get("mito_min_target_coverage", .95)),
                    f"Removal violates coverage thresholds: {report_path}")
            recorded_removed.add(row["assembly_id"])
    require(recorded_removed == removed, f"Integration report/list disagree: {report_path}")
    added = ([f"{sample}_mito"] if len(mito) == 1 else
             [f"{sample}_mito_{number}" for number in range(1, len(mito) + 1)])
    expected = {name: seq for name, seq in before.items() if name not in removed}
    for name, sequence in zip(added, mito.values()):
        require(name not in expected, f"Mitochondrial header collides with retained input ID: {name}")
        expected[name] = sequence
        require(raw_final.get(name) == sequence, f"Added mitochondrial record differs from OATK: {name}")
    require(raw_final == expected, f"Final FASTA differs from retained input plus OATK: {sample}/{assembler}")
    # Compare to every final record, in both orientations, allowing a shifted origin.
    duplicate_details = []
    for name in added:
        same = [other for other, sequence in raw_final.items() if circular_equal(raw_final[name], sequence)]
        require(same == [name], f"Circular mitochondrial duplicate of {name}: {same}")
        duplicate_details.append({"mito_id": name, "length_bp": len(raw_final[name]), "matching_final_ids": same})
    return {"input_assembly": str(before_path), "removed_ids": sorted(removed),
            "added_ids": added, "circular_duplicate_checks": duplicate_details}


def merqury_check(directory, assembly_path, stats, sample, stage, assembler, project, k):
    label = f"{sample}.{stage}.{assembler}"
    qv_path, comp_path = nonempty(directory / "merqury.qv"), nonempty(directory / "completeness.stats")
    qv = [line.split() for line in qv_path.read_text().splitlines() if line.strip()]
    comp = [line.split() for line in comp_path.read_text().splitlines() if line.strip()]
    require(len(qv) == len(comp) == 1 and len(qv[0]) == len(comp[0]) == 5,
            f"Expected one assembly per Merqury report: {directory}")
    require(qv[0][0] == comp[0][0] == label and comp[0][1] == "all",
            f"Merqury input label mismatch: {directory}")
    errors, total = int(qv[0][1]), int(qv[0][2])
    quality, rate = float(qv[0][3]), float(qv[0][4])
    require(0 <= errors <= total and 0 < total <= stats["maximum_possible_kmer_windows"],
            f"Merqury k-mer total exceeds assembly windows or counts are invalid: {directory}")
    require(not math.isnan(quality) and quality >= 0 and math.isfinite(rate) and 0 <= rate <= 1,
            f"Invalid Merqury quality/error rate: {directory}")
    require(math.isfinite(quality) or errors == 0, f"Infinite Merqury QV with nonzero errors: {directory}")
    recovered, solid = int(comp[0][2]), int(comp[0][3])
    percent = float(comp[0][4])
    require(solid > 0 and 0 <= recovered <= solid and math.isfinite(percent) and 0 <= percent <= 100,
            f"Invalid Merqury completeness: {directory}")
    require(abs(percent - 100 * recovered / solid) < .01, f"Inconsistent Merqury completeness: {directory}")
    manifest = json.loads(nonempty(directory / "inputs.json").read_text())
    expected_reads = project / "Assemblies" / sample / "QC/merqury" / f"{sample}.reads.meryl"
    require(Path(manifest["assembly"]).resolve() == assembly_path.resolve(), f"Wrong Merqury assembly provenance: {directory}")
    if "assembly_sha256" in manifest:
        require(manifest["assembly_sha256"] == stats["sha256"],
                f"Merqury describes changed FASTA contents: {directory}")
    if manifest.get("schema_version", 1) >= 2:
        require(manifest.get("assembly_sha256") == stats["sha256"],
                f"Merqury receipt lacks its assembly hash: {directory}")
        tools = manifest.get("tools", {})
        require(Path(tools["merqury_executable"]).parent == Path(tools["merqury_root"]),
                f"Merqury executable does not match its recorded source root: {directory}")
        for source, expected_hash in tools["source_sha256"].items():
            require(hashlib.sha256(nonempty(source).read_bytes()).hexdigest() == expected_hash,
                    f"Merqury source/binary changed: {source}")
        adaptation = tools.get("plot_adaptation")
        if adaptation:
            original = nonempty(adaptation["original_path"])
            patched = nonempty(directory / adaptation["patched_path"])
            require(hashlib.sha256(original.read_bytes()).hexdigest() == adaptation["original_sha256"],
                    f"Merqury original plotting source changed: {directory}")
            require(hashlib.sha256(patched.read_bytes()).hexdigest() == adaptation["patched_sha256"],
                    f"Merqury retained plotting adaptation changed: {directory}")
    require(Path(manifest["read_db"]).resolve() == expected_reads.resolve(), f"Wrong Merqury read provenance: {directory}")
    require(manifest["label"] == label and int(manifest["k"]) == k and
            manifest["mode"] == "one assembly per Merqury call", f"Wrong Merqury mode/settings: {directory}")
    require((directory / f"{label}.fasta").samefile(assembly_path), f"Wrong Merqury assembly link: {directory}")
    require((directory / "illumina.reads.meryl").samefile(expected_reads), f"Wrong Merqury reads link: {directory}")
    plot = png_info(directory / "merqury.spectra-asm.ln.png")
    return {"label": label, "errors": errors, "total_assembly_kmers": total,
            "qv": quality if math.isfinite(quality) else "inf", "error_rate": rate,
            "recovered_solid_kmers": recovered, "total_solid_kmers": solid,
            "completeness_percent": percent, "plot": plot, "inputs": manifest}


def nquire_check(directory, sample, expected_bam):
    """Independently reconcile native bins, fixed-model ranks and status receipts."""
    results = nonempty(directory / f"{sample}_nquire.results")
    ploidy = nonempty(directory / f"{sample}_ploidy.tsv")
    receipt_path = directory / f"{sample}_nquire.summary.json"
    receipt = json.loads(nonempty(receipt_path).read_text()) if receipt_path.exists() else None

    def positions(path):
        with nonempty(path).open("rb") as stream:
            header = stream.read(9)
            require(len(header) == 9, f"Truncated nQuire header: {path}")
            count, extended, minimum = struct.unpack(">IBI", header)
            require(extended in (0, 1) and 0 <= minimum < 500, f"Invalid nQuire header: {path}")
            if extended:
                require(len(stream.read(4)) == 4, f"Truncated nQuire hash: {path}")
            for _ in range(count):
                if extended:
                    require(len(stream.read(8)) == 8, f"Truncated nQuire coordinates: {path}")
                record = stream.read(3)
                require(len(record) == 3, f"Truncated nQuire record: {path}")
                coverage, alleles = struct.unpack(">HB", record)
                counts = stream.read(2 * alleles)
                require(coverage > 0 and 1 <= alleles <= 6 and len(counts) == 2 * alleles,
                        f"Invalid nQuire allele record: {path}")
                support = struct.unpack(">" + "H" * alleles, counts)
                require(all(0 < value <= coverage for value in support) and sum(support) <= coverage,
                        f"Invalid nQuire support: {path}")
            require(not stream.read(1), f"Trailing nQuire bin data: {path}")
        return count

    raw_path, denoised_path = directory / f"{sample}.bin", directory / f"{sample}-denoised.bin"
    raw_count = positions(raw_path)
    status = receipt["status"] if receipt else "legacy_ranked_fixed_models"
    denoised_count = positions(denoised_path) if denoised_path.exists() else None
    if status.startswith("untestable_"):
        require(receipt and receipt["state"] == "completed", "Untestable nQuire result lacks a completed receipt")
        require((status == "untestable_zero_raw_positions" and raw_count == 0 and denoised_count is None)
                or (status == "untestable_zero_denoised_positions" and raw_count > 0 and denoised_count == 0),
                "nQuire untestable status disagrees with native position counts")
        require(results.read_text().startswith(f"# NOT_RUN: {status}\n") and
                not receipt["lrdmodel_executed"] and receipt["best_supported_fixed_model"] is None,
                "nQuire no-site result masquerades as an executed model comparison")
        best, deltas = None, None
    else:
        require(raw_count > 0 and denoised_count is not None and 0 < denoised_count <= raw_count,
                "nQuire model comparison has no valid native site support")
        rows = [line.split("\t") for line in results.read_text().splitlines() if line.strip()]
        header = ["file", "free", "dip", "tri", "tet", "d_dip", "d_tri", "d_tet"]
        require(len(rows) == 2 and rows[0] == header and len(rows[1]) == 8 and
                rows[1][0] == sample + "-denoised.bin", "nQuire model report identity/schema mismatch")
        values = dict(zip(header[1:], map(float, rows[1][1:])))
        require(all(math.isfinite(value) for value in values.values()), "nQuire likelihoods are not finite")
        for model in ["dip", "tri", "tet"]:
            require(math.isclose(values["free"] - values[model], values["d_" + model], rel_tol=1e-10, abs_tol=2e-6),
                    "nQuire likelihood delta is inconsistent")
        deltas = {"2": values["d_dip"], "3": values["d_tri"], "4": values["d_tet"]}
        tied = [model for model, value in deltas.items() if value == min(deltas.values())]
        best = tied[0] if len(tied) == 1 else None
        if receipt:
            require(status == ("ranked_fixed_models" if best else "ambiguous_tie") and
                    receipt["lrdmodel_executed"] and receipt["best_supported_fixed_model"] == best and
                    receipt["fixed_model_deltas"] == deltas, "nQuire receipt contradicts native fixed-model evidence")
    require(ploidy.read_text().split() == [sample, best or "NA"], "nQuire fixed-model TSV mismatch")
    artifacts = [results, ploidy]
    if receipt:
        require(receipt["state"] == "completed" and receipt["sample"] == sample and
                receipt["interpretation"] == "exploratory_fixed_model_comparison_not_validated_biological_ploidy",
                "nQuire receipt state, identity or interpretation mismatch")
        require(Path(receipt["bam"]["path"]).resolve() == expected_bam.resolve() and
                receipt["bam"]["sha256"] == file_sha256(expected_bam), "nQuire BAM provenance mismatch")
        require(receipt["raw_bin"]["positions"] == raw_count and
                (receipt["denoised_bin"]["positions"] if receipt["denoised_bin"] else None) == denoised_count,
                "nQuire receipt position counts mismatch")
        script = Path(__file__).resolve().parents[1] / "workflow/scripts/run_nquire.py"
        require(Path(receipt["implementation"]["script"]).name == script.name and
                receipt["implementation"]["sha256"] == file_sha256(script), "nQuire implementation mismatch")
        require({results.name, ploidy.name, raw_path.name} <= set(receipt["outputs"]),
                "nQuire receipt omits required native artifacts")
        for name, record in receipt["outputs"].items():
            path = safe_receipt_path(directory, name)
            require(path.stat().st_size == record["bytes"] and file_sha256(path) == record["sha256"],
                    f"nQuire artifact hash/size mismatch: {name}")
        artifacts.append(receipt_path)
    marker(directory / "nquire.done", artifacts)
    return {"status": status, "raw_positions": raw_count, "denoised_positions": denoised_count,
            "best_supported_fixed_model": best, "fixed_model_deltas": deltas,
            "execution_provenance": "receipt_checked" if receipt else "legacy_no_receipt"}


def genomescope_check(path):
    with nonempty(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(rows, f"GenomeScope summary has no model rows: {path}")
    result = []
    scores = []
    for row in rows:
        all_fit, full_fit = float(row["model_fit_all_percent"]), float(row["model_fit_full_percent"])
        require(all(math.isfinite(value) and value <= 100 for value in (all_fit, full_fit)),
                f"Invalid GenomeScope fit values: {path}")
        require(row["selection_status"] == "exploratory_fit_ranking_not_validated_ploidy",
                f"GenomeScope ranking presented as confirmed ploidy: {path}")
        require("best" not in row, f"Obsolete GenomeScope definitive best flag: {path}")
        source = nonempty(row["source_summary"])
        source_fit = re.search(r"(?m)^Model Fit\s+([-+0-9.eE]+)%\s+([-+0-9.eE]+)%", source.read_text())
        require(source_fit and [float(x) for x in source_fit.groups()] == [all_fit, full_fit],
                f"GenomeScope summary disagrees with source: {path}")
        result.append({"ploidy_model": int(row["ploidy"]), "all_fit_percent": all_fit,
                       "full_fit_percent": full_fit, "fit_rank": int(row["fit_rank"]),
                       "best_by_model_fit": row["best_by_model_fit"], "source_summary": str(source)})
        scores.append(all_fit)
    ranked = sorted(set(scores), reverse=True)
    for row in result:
        rank = ranked.index(row["all_fit_percent"]) + 1
        require(row["fit_rank"] == rank and row["best_by_model_fit"] == ("YES" if rank == 1 else ""),
                f"GenomeScope fit ranking is inconsistent: {path}")
    return result


def gff_check(path):
    content = nonempty(path).read_text()
    lines = content.splitlines()
    require(lines and lines[0].strip() == "##gff-version 3", f"Missing GFF3 header or leading diagnostics: {path}")
    count = 0
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        require(len(fields) == 9, f"Non-GFF diagnostic or wrong field count: {path}")
        require(int(fields[3]) > 0 and int(fields[4]) >= int(fields[3]), f"Invalid GFF coordinates: {path}")
        count += 1
    return {"path": str(path), "feature_count": count}


def quast_check(path, assemblers):
    with nonempty(path).open() as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    expected = [f"{stage}_{assembler}" for stage in ["raw", "final"] for assembler in assemblers]
    require(rows and rows[0][0] == "Assembly" and rows[0][1:] == expected,
            f"QUAST labels do not identify the requested assemblies: {path}")
    require(len(rows) > 1 and all(len(row) == len(expected) + 1 for row in rows), f"Malformed QUAST report: {path}")
    return {"path": str(path), "labels": expected,
            "metrics": {row[0]: dict(zip(expected, row[1:])) for row in rows[1:]}}


def busco_check(directory):
    summaries = sorted(path for path in directory.glob("short_summary*.txt") if path.is_file())
    require(summaries, f"No direct BUSCO short summary: {directory}")
    parsed = []
    expression = r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:(\d+)"
    for source in summaries:
        content = nonempty(source).read_text()
        match = re.search(expression, content)
        require(match is not None, f"Unparseable BUSCO summary: {source}")
        complete, single, duplicate, fragmented, missing = map(float, match.groups()[:5])
        total = int(match.group(6))
        require(total > 0 and all(math.isfinite(x) and 0 <= x <= 100 for x in
                (complete, single, duplicate, fragmented, missing)), f"Invalid BUSCO scores: {source}")
        require(abs(complete + fragmented + missing - 100) <= .3 and
                abs(single + duplicate - complete) <= .3, f"Inconsistent BUSCO scores: {source}")
        lineage = re.search(r"lineage dataset is:\s*(\S+)", content)
        parsed.append({"path": str(source), "lineage": lineage.group(1) if lineage else None,
                       "complete_percent": complete, "single_percent": single,
                       "duplicated_percent": duplicate, "fragmented_percent": fragmented,
                       "missing_percent": missing, "markers": total})
    return parsed


def mapping_check(directory, assembler, samtools):
    bam = nonempty(directory / f"purged-{assembler}.bam")
    nonempty(directory / f"purged-{assembler}.bam.bai")
    flagstat = nonempty(directory / f"purged-{assembler}.flagstat.txt")
    counts = {}
    for line in flagstat.read_text().splitlines():
        match = re.match(r"(\d+) \+ (\d+) (.+)", line)
        require(match is not None, f"Malformed flagstat line: {flagstat}")
        description = match.group(3).split(" (")[0]
        counts[description] = {"qc_passed": int(match.group(1)), "qc_failed": int(match.group(2))}
    require("paired in sequencing" in counts and "properly paired" in counts and "mapped" in counts,
            f"Incomplete paired-read flagstat: {flagstat}")
    denominator = counts["paired in sequencing"]["qc_passed"]
    numerator = counts["properly paired"]["qc_passed"]
    require(denominator > 0 and 0 <= numerator <= denominator, f"Invalid paired-read counts: {flagstat}")
    status = "not_run_samtools_unavailable"
    if samtools:
        checked = subprocess.run([samtools, "quickcheck", "-v", str(bam)], capture_output=True, text=True, timeout=60)
        require(checked.returncode == 0, f"BAM quickcheck failed: {bam}: {checked.stderr} {checked.stdout}")
        status = "passed"
    return {"bam": str(bam), "flagstat": str(flagstat), "counts": counts,
            "properly_paired_percent": 100 * numerator / denominator, "bam_quickcheck": status}


def multiqc_data_check(data, project, samples):
    """Require every requested report identity and its exact source provenance.

    Generic and specific BUSCO lineage summaries remain distinct rows. Their
    count can exceed the logical sample/stage/assembler bundle count.
    """
    with nonempty(data / "multiqc_sources.txt").open() as handle:
        sources = list(csv.DictReader(handle, delimiter="\t"))
    expected_quast, expected_mapping, expected_busco = {}, {}, {}
    for sample in samples:
        name = sample["sample"].strip()
        qc = project / "Assemblies" / name / "QC"
        for assembler in selected_assemblers(sample):
            expected_mapping[(name, assembler)] = qc / "Mapping" / f"purged-{assembler}.flagstat.txt"
            for stage in ["raw", "final"]:
                expected_quast[(name, f"{stage}_{assembler}")] = qc / "QUAST/report.tsv"
                busco_name = assembler if stage == "raw" else f"purged-{assembler}"
                expected_busco[(name, stage, assembler)] = qc / f"BUSCO-{busco_name}" / f"{name}-{busco_name}"
    result = {}
    for module, filename, expected in [
            ("quast", "multiqc_quast.txt", expected_quast),
            ("samtools", "multiqc_samtools_flagstat.txt", expected_mapping),
            ("busco", "multiqc_busco.txt", expected_busco)]:
        with nonempty(data / filename).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        require(rows, f"Empty MultiQC parsed {module} data: {data}")
        names = [row["Sample"] for row in rows]
        require(len(set(names)) == len(names), f"Duplicate MultiQC {module} sample names: {data}")
        matched, records = set(), []
        for row in rows:
            display = row["Sample"]
            provenance = [entry for entry in sources if entry["Module"].lower() == module and
                          entry["Sample Name"] == display and
                          (module != "samtools" or entry["Section"].lower() == "flagstat")]
            paths = {Path(entry["Source"]).resolve() for entry in provenance}
            require(len(paths) == 1, f"Missing or ambiguous MultiQC {module} source for {display}")
            source = paths.pop()
            require("auto_lineage" not in source.parts,
                    f"MultiQC includes an automatic-lineage intermediate: {source}")
            if module == "busco":
                candidates = [key for key, directory in expected.items() if source.parent == directory.resolve()]
                require(source.name.startswith("short_summary"), f"Unexpected BUSCO source: {source}")
            else:
                suffix = display.split(" | ")[-1]
                candidates = [key for key, path in expected.items() if source == path.resolve() and
                              suffix == (key[1] if module == "quast" else f"purged-{key[1]}")]
            require(len(candidates) == 1, f"Unrecognized or merged MultiQC {module} identity: {display}: {source}")
            identity = candidates[0]
            require(identity[0] in display, f"MultiQC display name lost strain identity: {display}")
            matched.add(identity)
            if module == "quast":
                with nonempty(source).open() as handle:
                    original = list(csv.reader(handle, delimiter="\t"))
                index = original[0].index(identity[1])
                original_values = {values[0]: values[index] for values in original[1:] if values}
                for metric in ["N50", "Total length", "# contigs"]:
                    require(metric in row and metric in original_values and
                            float(row[metric]) == float(original_values[metric]),
                            f"MultiQC QUAST {metric} disagrees with {source} / {identity[1]}")
            elif module == "samtools":
                original = {}
                for line in nonempty(source).read_text().splitlines():
                    match = re.match(r"(\d+) \+ (\d+) (.+)", line)
                    if match:
                        original[match.group(3).split(" (")[0]] = int(match.group(1))
                for exported, source_key in [("total_passed", "in total"),
                                              ("mapped_passed", "mapped"),
                                              ("properly paired_passed", "properly paired")]:
                    require(exported in row and source_key in original and
                            int(float(row[exported])) == original[source_key],
                            f"MultiQC mapping {exported} disagrees with {source}")
            else:
                content = nonempty(source).read_text()
                complete = re.search(r"(?m)^\s*(\d+)\s+Complete BUSCOs", content)
                lineage = re.search(r"lineage dataset is:\s*(\S+)", content)
                require(complete and float(row["complete"]) == int(complete.group(1)),
                        f"MultiQC BUSCO complete count disagrees with {source}")
                require(lineage and row["lineage_dataset"] == lineage.group(1),
                        f"MultiQC BUSCO lineage identity disagrees with {source}")
            records.append({"display_name": display, "identity": list(identity),
                            "source": str(source), "lineage": row.get("lineage_dataset")})
        require(matched == set(expected),
                f"MultiQC {module} lost report identities: missing {sorted(set(expected) - matched)}")
        if module != "busco":
            require(len(rows) == len(expected), f"MultiQC {module} row count differs from requested assemblies")
        result[module] = {"row_count": len(rows), "expected_bundle_count": len(expected),
                          "observed_bundle_count": len(matched), "reports": records}
    return result


# Inserted into tests/validate_full_run.py; independent of runtime helpers.
_HASH_CACHE = {}


def file_sha256(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns)
    if key not in _HASH_CACHE:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(chunk)
        _HASH_CACHE[key] = digest.hexdigest()
    return _HASH_CACHE[key]


def tsv_rows(path):
    with nonempty(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
                f"Missing or duplicate TSV fields: {path}")
        rows = list(reader)
    require(all(None not in row and all(value is not None for value in row.values()) for row in rows),
            f"Malformed TSV rows: {path}")
    return rows


def close_number(actual, expected, name, allow_infinite=False):
    if expected is None:
        require(actual is None or actual == "NA", f"Expected explicit NA for {name}, found {actual}")
        return
    require(actual is not None and actual not in {"", "NA"}, f"Missing numerical value for {name}")
    value = float(actual)
    if math.isinf(float(expected)):
        require(allow_infinite and value == float(expected), f"Unexpected infinity for {name}")
    else:
        require(math.isfinite(value) and math.isclose(value, float(expected), rel_tol=1e-8, abs_tol=1e-8),
                f"Metric mismatch for {name}: {actual} versus {expected}")


def safe_receipt_path(folder, name):
    path = (folder / name).resolve()
    require(path.is_relative_to(folder.resolve()), f"Receipt path escapes output directory: {name}")
    require(path.is_file(), f"Missing receipt output: {path}")
    return path


def interval_union_length(intervals):
    total, right = 0, -1
    for start, end in sorted(intervals):
        total += max(0, end - max(start, right))
        right = max(right, end)
    return total


def motif_forms(motif):
    require(re.fullmatch("[ACGT]+", motif) is not None, f"Invalid motif: {motif}")
    reverse = motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return {value[index:] + value[:index] for value in [motif, reverse] for index in range(len(value))}


def tidk_check(directory, assembly, sequences, config):
    run_config = json.loads(nonempty(directory / "run_config.json").read_text())
    receipt = json.loads(nonempty(directory / "complete.json").read_text())
    require(Path(run_config["input_path"]).resolve() == assembly.resolve() and
            run_config["input_sha256"] == file_sha256(assembly), f"TIDK input identity mismatch: {directory}")
    require(re.search(r"(?:^|\s)0\.2\.7(?:\s|$)", run_config["tidk_version"]), "Unexpected TIDK version")
    require(file_sha256(run_config["tidk"]) == run_config["tidk_sha256"], "TIDK executable checksum mismatch")
    script = Path(__file__).resolve().parents[1] / "workflow/scripts/tidk_telomeres.py"
    require(file_sha256(script) == run_config["script_sha256"], "TIDK implementation checksum mismatch")
    parameters = {"minimum": 5, "maximum": 30, "top_n": 10, "end_size": 10000,
                  "window": 10000, "copies": 5, "threshold": 20, "distance": .05,
                  "min_sequence_length": 0}
    for key, default in parameters.items():
        close_number(run_config[key], config.get("tidk_" + key, default), f"TIDK {key}")
    required = {"assembly_summary.tsv", "candidate_motifs.tsv", "sequence_inclusion.tsv",
                "telomere_counts_by_sequence.tsv", "telomere_validation_summary.tsv",
                "telomere_validation_by_sequence.tsv", "repeat_arrays.tsv", "run_config.json"}
    require(required <= set(receipt["output_sha256"]), "TIDK receipt omits required outputs")
    for name, digest in receipt["output_sha256"].items():
        require(file_sha256(safe_receipt_path(directory, name)) == digest, f"TIDK output checksum mismatch: {name}")
    require(file_sha256(directory / "screened_sequences.fasta") == run_config["screened_fasta_sha256"],
            "TIDK screened FASTA checksum mismatch")
    included = {name: seq for name, seq in sequences.items() if len(seq) > int(run_config["min_sequence_length"])}
    screened = fasta(directory / "screened_sequences.fasta") if included else {}
    require(screened == included, "TIDK screened sequences differ from exact retained input records")
    inclusion = tsv_rows(directory / "sequence_inclusion.tsv")
    require(len(inclusion) == len(sequences) and {row["sequence"] for row in inclusion} == set(sequences),
            "TIDK sequence inclusion identities disagree with input")
    for row in inclusion:
        name = row["sequence"]
        require(int(row["sequence_length"]) == len(sequences[name]) and
                int(row["included"]) == int(name in included), "TIDK inclusion length/flag mismatch")
        require(int(row["effective_end_window_bp"]) == min(int(run_config["end_size"]), len(sequences[name]) // 2),
                "TIDK terminal windows overlap or have incorrect size")
    candidates = tsv_rows(directory / "candidate_motifs.tsv")
    motifs = {row["motif"] for row in candidates}
    require(len(motifs) == len(candidates), "Duplicate TIDK candidate motifs")
    for row in candidates:
        motif = row["motif"]
        require(int(row["repeat_length"]) == len(motif) and motif == min(motif_forms(motif)) and
                int(row["explore_score"]) >= 0, "Invalid TIDK motif identity or score")
        require(not any(len(motif) % size == 0 and motif[:size] * (len(motif) // size) == motif
                        for size in range(1, len(motif))), "TIDK motif is not primitive")
    spans, array_seen = {}, set()
    arrays = tsv_rows(directory / "repeat_arrays.tsv")
    for row in arrays:
        name, motif = row["sequence"], row["motif"]
        require(name in included and motif in motifs, "TIDK array references an unknown retained sequence/motif")
        start, end = int(row["start_1based"]) - 1, int(row["end_1based"])
        require(0 <= start < end <= len(included[name]), "TIDK repeat coordinates exceed the input sequence")
        identity = (name, motif, start, end)
        require(identity not in array_seen, "Duplicate TIDK repeat array")
        array_seen.add(identity)
        piece = included[name][start:end]
        unit = piece[:len(motif)]
        require(unit in motif_forms(motif) and
                piece == (unit * ((len(piece) + len(unit) - 1) // len(unit)))[:len(piece)],
                "TIDK claimed array is not an exact phase/strand-equivalent tandem repeat")
        copies = (end - start) // len(motif)
        require(int(row["array_bp"]) == end - start and int(row["complete_copies"]) == copies and
                copies >= int(run_config["copies"]), "TIDK array copy/length mismatch")
        require(int(row["distance_left_bp"]) == start and
                int(row["distance_right_bp"]) == len(included[name]) - end, "TIDK array end-distance mismatch")
        window = min(int(run_config["end_size"]), len(included[name]) // 2)
        expected_region = ("BOTH_WINDOWS" if start < window and end > len(included[name]) - window else
                           "LEFT" if start < window else "RIGHT" if end > len(included[name]) - window else "INTERNAL")
        require(row["region"] == expected_region, "TIDK repeat region mismatch")
        spans.setdefault((name, motif), []).append((start, end))
    def validate_details(rows, expected_motifs):
        by_key = {(row["sequence"], row["motif"]): row for row in rows}
        require(len(by_key) == len(rows) and set(by_key) == {(name, motif) for name in included for motif in expected_motifs},
                "TIDK per-sequence/motif identities are incomplete or duplicated")
        for (name, motif), row in by_key.items():
            length = len(included[name])
            window = min(int(run_config["end_size"]), length // 2)
            require(int(row["sequence_length"]) == length and int(row["effective_end_window_bp"]) == window,
                    "TIDK per-sequence length/window mismatch")
            detected = 0
            for side in ["left", "right"]:
                eligible = [span for span in spans.get((name, motif), []) if
                            (span[0] < window if side == "left" else span[1] > length - window)]
                closest = min(eligible, key=lambda span: (span[0] if side == "left" else length - span[1],
                                                        -(span[1] - span[0]))) if eligible else None
                flag = int(row[side + "_detected"])
                require(flag == int(closest is not None), "TIDK terminal-array detection flag mismatch")
                detected += flag
                if closest:
                    start, end = closest
                    for key, value in [("start_1based", start + 1), ("end_1based", end),
                                       ("array_bp", end - start), ("complete_copies", (end - start) // len(motif)),
                                       ("distance_bp", start if side == "left" else length - end)]:
                        require(int(row[side + "_" + key]) == value, "TIDK selected terminal array mismatch")
                else:
                    require(all(row[side + "_" + key] == "" for key in
                                ["start_1based", "end_1based", "array_bp", "complete_copies", "distance_bp"]),
                            "TIDK absent end has unexplained coordinates")
            require(int(row["detected_ends"]) == detected <= 2, "TIDK detected-end count mismatch")
        return by_key
    details = validate_details(tsv_rows(directory / "telomere_validation_by_sequence.tsv"), motifs)
    validation = tsv_rows(directory / "telomere_validation_summary.tsv")
    require(len(validation) == len(motifs) and {row["motif"] for row in validation} == motifs,
            "TIDK motif-level summaries are incomplete")
    total_bp = sum(map(len, included.values()))
    terminal_space = sum(2 * min(int(run_config["end_size"]), len(seq) // 2) for seq in included.values())
    for row in validation:
        motif = row["motif"]
        detail = [details[(name, motif)] for name in included]
        ends = sum(int(value["detected_ends"]) for value in detail)
        terminal_bp = internal_bp = 0
        for name, sequence in included.items():
            window = min(int(run_config["end_size"]), len(sequence) // 2)
            regions = spans.get((name, motif), [])
            terminal = [(a, min(b, window)) for a, b in regions if a < window]
            terminal += [(max(a, len(sequence) - window), b) for a, b in regions if b > len(sequence) - window]
            internal = [(max(a, window), min(b, len(sequence) - window)) for a, b in regions
                        if min(b, len(sequence) - window) > max(a, window)]
            terminal_bp += interval_union_length(terminal)
            internal_bp += interval_union_length(internal)
        terminal_density = terminal_bp / terminal_space
        internal_density = internal_bp / (total_bp - terminal_space) if total_bp > terminal_space else None
        for field, expected in [("terminal_ends_detected", ends), ("terminal_repeat_bp", terminal_bp),
                                ("internal_repeat_bp", internal_bp),
                                ("both_ends", sum(int(value["detected_ends"]) == 2 for value in detail)),
                                ("one_end", sum(int(value["detected_ends"]) == 1 for value in detail)),
                                ("zero_ends", sum(int(value["detected_ends"]) == 0 for value in detail))]:
            close_number(row[field], expected, "TIDK " + field)
        close_number(row["terminal_density"], terminal_density, "TIDK terminal density")
        require(row["internal_density"] == "" if internal_density is None else
                math.isclose(float(row["internal_density"]), internal_density, rel_tol=1e-8, abs_tol=1e-8),
                "TIDK internal density mismatch")
        require(int(row["supported"]) == int(ends > 0 and internal_density is not None and terminal_density > internal_density),
                "TIDK exploratory support rule mismatch")
    summary_rows = tsv_rows(directory / "assembly_summary.tsv")
    require(len(summary_rows) == 1, "TIDK must have one assembly summary")
    summary = summary_rows[0]
    best, status = summary["best_motif"], summary["screen_status"]
    supported = {row["motif"] for row in validation if row["supported"] == "1"}
    require((best in supported and status == "CANDIDATE_IDENTIFIED") if supported else
            (not best and status == ("NO_TELOMERIC_REPEAT_IDENTIFIED" if included else "NO_SEQUENCES_ABOVE_THRESHOLD")),
            "TIDK best-motif status does not describe its evidence")
    final_details = tsv_rows(directory / "telomere_counts_by_sequence.tsv")
    if best:
        validate_details(final_details, {best})
    else:
        require(len(final_details) == len(included) and {row["sequence"] for row in final_details} == set(included),
                "TIDK nondetection omits retained sequences")
        require(all(row["motif"] == "" and int(row["detected_ends"]) == int(row["left_detected"]) ==
                    int(row["right_detected"]) == 0 for row in final_details), "TIDK nondetection contains unexplained detections")
    require(all(row["screen_status"] == status for row in final_details), "TIDK sequence statuses disagree with assembly status")
    expected_counts = {"input_sequences": len(sequences), "retained_sequences": len(included),
                       "excluded_sequences": len(sequences) - len(included), "candidate_motifs_tested": len(motifs),
                       "detected_ends": sum(int(row["detected_ends"]) for row in final_details),
                       "possible_sequence_ends": 2 * len(included),
                       "both_ends": sum(int(row["detected_ends"]) == 2 for row in final_details),
                       "one_end": sum(int(row["detected_ends"]) == 1 for row in final_details),
                       "zero_ends": sum(int(row["detected_ends"]) == 0 for row in final_details)}
    for field, expected in expected_counts.items():
        close_number(summary[field], expected, "TIDK assembly " + field)
        close_number(receipt[field], expected, "TIDK receipt " + field)
    require(receipt["screen_status"] == status and receipt["best_motif"] == best and
            receipt["interpretation"] == "EXPLORATORY_EXACT_REPEAT_CANDIDATES_NOT_VALIDATED_TELOMERES",
            "TIDK receipt status/interpretation mismatch")
    commands = run_config["commands"]
    require(len(commands) == (1 + len(motifs) if included else 0) and
            all(Path(command[0]).resolve() == Path(run_config["tidk"]).resolve() for command in commands),
            "TIDK command provenance does not match retained candidate searches")
    marker(directory / "tidk.done", [directory / "assembly_summary.tsv", directory / "complete.json"])
    return {"summary": summary, "input_sha256": run_config["input_sha256"], "tool_version": run_config["tidk_version"],
            "arrays_checked": len(arrays), "output_hashes_checked": len(receipt["output_sha256"]),
            "interpretation": receipt["interpretation"]}


def checked_bed(path, sequences):
    path = Path(path)
    require(path.is_file(), f"Missing BED: {path}")
    intervals, total, last_order = {}, 0, -1
    order = {name: index for index, name in enumerate(sequences)}
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            require(len(fields) == 3 and fields[0] in sequences, f"Unknown/malformed BED interval: {path}")
            name, start, end = fields[0], int(fields[1]), int(fields[2])
            require(order[name] >= last_order and 0 <= start < end <= len(sequences[name]), f"Unsorted/out-of-range BED: {path}")
            require(not intervals.get(name) or start >= intervals[name][-1][1], f"Overlapping BED intervals: {path}")
            intervals.setdefault(name, []).append((start, end))
            total += end - start
            last_order = order[name]
    return intervals, total


def interval_lookup(intervals):
    return {name: ([start for start, end in spans], [end for start, end in spans]) for name, spans in intervals.items()}


def position_in_mask(index, name, position1):
    if name not in index:
        return False
    starts, ends = index[name]
    item = bisect_right(starts, position1 - 1) - 1
    return item >= 0 and position1 <= ends[item]


def vcf_snp_records(path, sequences, callable_index, minimum_dp, maximum_dp, passing):
    records, sample_header = {}, False
    with gzip.open(nonempty(path), "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                require(len(line.rstrip().split("\t")) == 10, f"Expected exactly one VCF sample: {path}")
                sample_header = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            require(sample_header, f"VCF is missing its sample header: {path}")
            fields = line.rstrip().split("\t")
            require(len(fields) == 10, f"Invalid single-sample VCF record: {path}")
            chrom, position, _, refbase, alt, qual, _, _, fmt, values = fields
            position, quality = int(position), float(qual)
            require(chrom in sequences and 1 <= position <= len(sequences[chrom]) and
                    refbase == sequences[chrom][position - 1] and refbase in "ACGT" and
                    len(alt) == 1 and alt in "ACGT" and alt != refbase,
                    f"VCF contains a non-originally-biallelic SNP or wrong reference allele: {path}")
            require(math.isfinite(quality) and quality >= 30 and position_in_mask(callable_index, chrom, position),
                    f"SNP QUAL or common callable-mask membership failed: {path}")
            key = (chrom, position)
            require(key not in records, f"Duplicate SNP position: {path}: {key}")
            genotype = dict(zip(fmt.split(":"), values.split(":")))
            if passing:
                require(genotype.get("GT") in {"0/1", "1/0", "0|1", "1|0"}, f"Non-diploid-heterozygous GT: {path}")
                depth = int(genotype["DP"])
                allele_depth = list(map(int, genotype["AD"].split(",")))
                require(len(allele_depth) == 2 and min(allele_depth) >= 3 and minimum_dp <= depth <= maximum_dp,
                        f"SNP DP/allele-depth filters failed: {path}")
            records[key] = {"REF": refbase, "ALT": alt, "QUAL": quality, **genotype}
    require(sample_header, f"Missing VCF header: {path}")
    return records


def heterozygosity_check(directory, project, sample, assembler, assembly, sequences, config, samtools):
    summary = json.loads(nonempty(directory / "summary.json").read_text())
    receipt = json.loads(nonempty(directory / "receipt.json").read_text())
    require(receipt["state"] == "completed" and receipt["metrics"] == summary,
            "Heterozygosity summary differs from completed receipt")
    require(summary["sample"] == sample and summary["assembler"] == assembler and
            summary["genotype_model"] == "STANDARDIZED_DIPLOID_MODEL_NOT_CONFIRMED_BIOLOGICAL_PLOIDY",
            "Heterozygosity sample/assembler/model identity mismatch")
    expected_inputs = {"reference": assembly,
        "bam": project / "Assemblies" / sample / "QC/Mapping" / f"purged-{assembler}.bam",
        "bam_index": project / "Assemblies" / sample / "QC/Mapping" / f"purged-{assembler}.bam.bai",
        "mitochondrial_source": project / "Assemblies" / sample / "OATK" / f"{sample}.mito.ctg.fasta"}
    for name, path in expected_inputs.items():
        record = receipt["inputs"][name]
        require(Path(record["path"]).resolve() == path.resolve() and file_sha256(path) == record["sha256"],
                f"Heterozygosity {name} input identity/checksum mismatch")
    required = {"summary.json", "summary.tsv", "callable.bed", "nuclear_acgt.bed", "excluded_contigs.txt",
                "depth_histogram.tsv", "biallelic_snps.vcf.gz", "biallelic_snps.vcf.gz.tbi",
                "heterozygous_snps.vcf.gz", "heterozygous_snps.vcf.gz.tbi", "heterozygous_snps.tsv.gz", "allele_balance.tsv"}
    require(required <= set(receipt["outputs"]), "Heterozygosity receipt omits required artifacts")
    for name, record in receipt["outputs"].items():
        path = safe_receipt_path(directory, name)
        require(path.stat().st_size == record["bytes"] and file_sha256(path) == record["sha256"],
                f"Heterozygosity output checksum/size mismatch: {name}")
    script = Path(__file__).resolve().parents[1] / "workflow/scripts/estimate_heterozygosity.py"
    # workflow.source_path executes a byte-identical cached helper. Its runtime
    # pathname is historical provenance, not an assertion of repository location.
    require(Path(receipt["implementation"]["script"]).name == script.name and
            file_sha256(script) == receipt["implementation"]["sha256"], "Heterozygosity implementation checksum mismatch")
    parameters = receipt["parameters"]
    for key, value in {"min_mq": 30, "min_bq": 20, "min_dp": 10, "max_dp_multiplier": 2,
                       "min_qual": 30, "min_ref_ad": 3, "min_alt_ad": 3, "ploidy": 2,
                       "pileup_max_depth": 10000}.items():
        close_number(parameters[key], value, "Heterozygosity " + key)
    require(parameters["keep_alts"] is True and parameters["GQ_filter"] is None and
            parameters["excluded_read_flags"] == "UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP",
            "Heterozygosity calling/read-filter method differs from declared model")
    threshold = float(config.get("heterozygosity_min_callable_fraction", .5))
    close_number(parameters["callability_review_threshold"], threshold, "callability review threshold")
    require(receipt["reference_consistency"]["SN_LN_and_order_match"] is True and
            "do not prove" in receipt["reference_consistency"]["identity_boundary"],
            "Heterozygosity lacks a qualified reference consistency check")
    if samtools:
        checked = subprocess.run([samtools, "view", "-H", str(expected_inputs["bam"])], capture_output=True, text=True, timeout=60)
        require(checked.returncode == 0, f"Cannot independently inspect heterozygosity BAM header: {checked.stderr}")
        header_lengths = []
        for line in checked.stdout.splitlines():
            if line.startswith("@SQ\t"):
                tags = dict(value.split(":", 1) for value in line.split("\t")[1:] if ":" in value)
                header_lengths.append((tags["SN"], int(tags["LN"])))
                if "M5" in tags:
                    require(hashlib.md5(sequences[tags["SN"]].encode()).hexdigest().lower() == tags["M5"].lower(),
                            "Heterozygosity BAM reference M5 mismatch")
        require(header_lengths == [(name, len(seq)) for name, seq in sequences.items()],
                "Heterozygosity BAM reference SN/LN/order mismatch")
    mito = fasta(expected_inputs["mitochondrial_source"])
    known_ids = ([f"{sample}_mito"] if len(mito) == 1 else [f"{sample}_mito_{i}" for i in range(1, len(mito) + 1)])
    for name, sequence in zip(known_ids, mito.values()):
        require(sequences.get(name) == sequence, "Heterozygosity excluded OATK sequence identity mismatch")
    extra = config.get("heterozygosity_exclude_contigs", {}).get(sample, {}).get(assembler, [])
    excluded = set(known_ids) | set(extra)
    require(excluded <= set(sequences) and set((directory / "excluded_contigs.txt").read_text().split()) == excluded and
            set(summary["excluded_contig_ids"]) == excluded and
            receipt["exclusions"]["known_appended_OATK_ids"] == known_ids and
            receipt["exclusions"]["additional_configured_ids"] == extra,
            "Heterozygosity explicit exclusions disagree with exact inputs/configuration")
    nuclear, nuclear_total = checked_bed(directory / "nuclear_acgt.bed", sequences)
    expected_nuclear = {name: [(match.start(), match.end()) for match in re.finditer("[ACGT]+", seq)]
                        for name, seq in sequences.items() if name not in excluded}
    expected_nuclear = {name: spans for name, spans in expected_nuclear.items() if spans}
    require(nuclear == expected_nuclear, "Heterozygosity nuclear ACGT mask does not exactly match reference/exclusions")
    callable_mask, callable_bases = checked_bed(directory / "callable.bed", sequences)
    nuclear_index = interval_lookup(nuclear)
    for name, spans in callable_mask.items():
        for start, end in spans:
            starts, ends = nuclear_index.get(name, ([], []))
            item = bisect_right(starts, start) - 1
            require(item >= 0 and end <= ends[item], "Callable BED includes excluded or ambiguous reference bases")
    hist = {}
    for row in tsv_rows(directory / "depth_histogram.tsv"):
        depth, positions = int(row["DP"]), int(row["Positions"])
        require(depth >= 0 and positions >= 0 and depth not in hist, "Invalid duplicate/negative depth histogram bin")
        hist[depth] = positions
    positive_bases = sum(count for depth, count in hist.items() if depth > 0)
    mean = sum(depth * count for depth, count in hist.items() if depth > 0) / positive_bases if positive_bases else 0
    maximum = max(10, int(2 * mean + .5))
    require(sum(count for depth, count in hist.items() if 10 <= depth <= maximum) == callable_bases and
            sum(hist.values()) <= nuclear_total, "Callable BED and evidence-depth histogram disagree")
    for field, expected in [("nuclear_acgt_bases", nuclear_total), ("callable_bases", callable_bases),
                            ("positions_with_nucleotide_evidence", sum(hist.values())),
                            ("mean_positive_nuclear_evidence_dp", mean), ("minimum_callable_dp", 10),
                            ("maximum_callable_dp", maximum), ("pileup_depth_cap", 10000),
                            ("positions_at_depth_cap", sum(count for depth, count in hist.items() if depth >= 10000))]:
        close_number(summary[field], expected, "heterozygosity " + field)
    index = interval_lookup(callable_mask)
    original = vcf_snp_records(directory / "biallelic_snps.vcf.gz", sequences, index, 10, maximum, False)
    passed = vcf_snp_records(directory / "heterozygous_snps.vcf.gz", sequences, index, 10, maximum, True)
    require(set(passed) <= set(original) and all(passed[key]["REF"] == original[key]["REF"] and
            passed[key]["ALT"] == original[key]["ALT"] for key in passed), "Heterozygous VCF differs from original biallelic set")
    with gzip.open(directory / "heterozygous_snps.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and {"Contig", "Position1", "REF", "ALT", "QUAL", "GT", "DP", "REF_AD", "ALT_AD", "ALT_fraction", "GQ"} <= set(reader.fieldnames),
                "Heterozygosity SNP TSV omits required fields")
        rows = list(reader)
    seen, balances, bins = set(), [], Counter()
    for row in rows:
        key = (row["Contig"], int(row["Position1"]))
        require(key in passed and key not in seen, "Heterozygosity TSV positions differ from authoritative VCF")
        seen.add(key)
        variant = passed[key]
        for field in ["REF", "ALT", "GT", "GQ"]:
            require(row[field] == variant[field], "Heterozygosity TSV allele/genotype differs from VCF")
        for field in ["QUAL", "DP"]:
            close_number(row[field], variant[field], "SNP TSV " + field)
        ref_ad, alt_ad = map(int, variant["AD"].split(","))
        require(int(row["REF_AD"]) == ref_ad and int(row["ALT_AD"]) == alt_ad, "SNP TSV allele depth mismatch")
        balance = alt_ad / (ref_ad + alt_ad)
        close_number(row["ALT_fraction"], balance, "SNP allele balance")
        balances.append(balance)
        bins[min(19, int(balance * 20))] += 1
    require(seen == set(passed), "Heterozygosity TSV omits called SNPs")
    balance_rows = tsv_rows(directory / "allele_balance.tsv")
    require(len(balance_rows) == 20, "Incomplete allele-balance histogram")
    for item, row in enumerate(balance_rows):
        close_number(row["ALT_fraction_start"], item / 20, "allele-balance bin start")
        close_number(row["ALT_fraction_end"], (item + 1) / 20, "allele-balance bin end")
        require(int(row["SNP_sites"]) == bins[item], "Allele-balance histogram and SNPs disagree")
    sites = len(passed)
    fraction = callable_bases / nuclear_total if nuclear_total else None
    expected_status = ("untestable_no_nuclear_acgt" if not nuclear_total else "untestable_no_callable_bases" if not callable_bases
                       else "low_callable_fraction" if fraction < threshold else "complete")
    require(summary["status"] == expected_status and summary["heterozygous_snps"] == sites,
            "Heterozygosity numerator or explicit callability status mismatch")
    for field, expected in [("callable_fraction", fraction), ("callability_review_threshold", threshold),
                            ("heterozygosity_percent", 100 * sites / callable_bases if callable_bases else None),
                            ("heterozygosity_proportion", sites / callable_bases if callable_bases else None),
                            ("het_snps_per_kb", 1000 * sites / callable_bases if callable_bases else None),
                            ("median_alt_fraction", statistics.median(balances) if balances else None)]:
        close_number(summary[field], expected, "heterozygosity " + field)
    table = tsv_rows(directory / "summary.tsv")
    require(len(table) == 1 and table[0]["sample"] == sample and table[0]["assembler"] == assembler,
            "Heterozygosity summary TSV identity mismatch")
    require(table[0]["status"] == summary["status"] and table[0]["genotype_model"] == summary["genotype_model"],
            "Heterozygosity summary TSV status/model mismatch")
    for field in ["callable_bases", "nuclear_acgt_bases", "heterozygous_snps", "heterozygosity_percent",
                  "callable_fraction", "maximum_callable_dp", "minimum_callable_dp", "median_alt_fraction"]:
        if summary[field] is None:
            require(table[0][field] == "", "Untestable summary TSV metric must remain empty under its explicit status")
        else:
            close_number(table[0][field], summary[field], "heterozygosity summary TSV " + field)
    marker(directory / "heterozygosity.done", [directory / "summary.json", directory / "receipt.json"])
    return {"summary": summary, "input_hashes_checked": len(expected_inputs), "output_hashes_checked": len(receipt["outputs"]),
            "unique_snp_positions_checked": sites, "header_independently_checked": bool(samtools),
            "validation_boundary": "Checks retained masks and variants, not a fresh genotype call or biological ploidy calibration"}


def qc_summary_check(directory, project, samples, native, config):
    table = directory / "assembly_qc.tsv"
    rows = tsv_rows(table)
    identities = [(row["sample"], row["assembler"], row["stage"]) for row in rows]
    expected = {(sample["sample"].strip(), assembler, stage) for sample in samples
                for assembler in selected_assemblers(sample) for stage in ["raw", "final"]}
    require(len(identities) == len(set(identities)) and set(identities) == expected,
            "QC summary drops, duplicates, or adds sample/assembler/stage identities")
    provenance = json.loads(nonempty(directory / "provenance.json").read_text())
    require(len(provenance["identities"]) == len(expected) and {tuple(identity) for identity in provenance["identities"]} == expected,
            "QC summary provenance identities disagree with table")
    script = Path(__file__).resolve().parents[1] / "workflow/scripts/qc_summary.py"
    require(provenance["script_sha256"] == file_sha256(script), "QC summary collector checksum mismatch")
    expected_sources, finite_qv, infinite_qv, paired_rows, tidk_testable, het_testable = set(), 0, 0, 0, 0, 0
    for row in rows:
        sample, assembler, stage = row["sample"], row["assembler"], row["stage"]
        key = f"{stage}/{assembler}"
        evidence = native[sample]
        stats = evidence["fasta"][key]
        require(stats is not None and Path(row["assembly"]).resolve() == Path(stats["path"]).resolve() and
                row["assembly_sha256"] == stats["sha256"], "QC summary assembly identity mismatch")
        expected_sources.add(Path(stats["path"]).resolve())
        fields = {"sequence_count": "sequence_count", "total_bp": "total_bp", "sequence_n50_bp": "n50_bp",
                  "longest_sequence_bp": "longest_bp", "gap_n_bp": "n_bp", "gc_percent_acgt": "gc_percent_excluding_ambiguity"}
        for field, source in fields.items():
            close_number(row[field], stats[source], "QC summary " + field)
        busco = evidence["busco"][key]
        require(busco is not None, "Cannot validate QC summary without native BUSCO evidence")
        specific = [report for report in busco if Path(report["path"]).name.startswith("short_summary.specific.")]
        selected = specific or [report for report in busco if Path(report["path"]).name.startswith("short_summary.generic.")]
        require(len(selected) == 1, "QC summary selected BUSCO lineage is ambiguous")
        selected = selected[0]
        require(Path(row["busco_report"]).resolve() == Path(selected["path"]).resolve() and
                row["busco_lineage"] == selected["lineage"], "QC summary BUSCO lineage/report identity mismatch")
        expected_sources.add(Path(selected["path"]).resolve())
        described = re.search(r"Summarized benchmarking in BUSCO notation for file (.+)", Path(selected["path"]).read_text())
        require(described and Path(described.group(1).strip()).resolve() == Path(stats["path"]).resolve(),
                "Selected BUSCO report describes a different assembly")
        for field, source in [("busco_markers", "markers"), ("busco_complete_percent", "complete_percent"),
                              ("busco_single_percent", "single_percent"), ("busco_duplicated_percent", "duplicated_percent"),
                              ("busco_fragmented_percent", "fragmented_percent"), ("busco_missing_percent", "missing_percent")]:
            close_number(row[field], selected[source], "QC summary " + field)
        merqury = evidence["merqury"][key]
        require(merqury is not None, "Cannot validate QC summary without native Merqury evidence")
        qv = float(merqury["qv"])
        close_number(row["merqury_qv"], qv, "QC summary Merqury QV", True)
        close_number(row["merqury_completeness_percent"], merqury["completeness_percent"], "QC summary Merqury completeness")
        close_number(row["merqury_k"], merqury["inputs"]["k"], "QC summary Merqury k")
        require(row["merqury_status"] == ("INFINITE_QV_ZERO_ERROR_KMERS" if math.isinf(qv) else "MEASURED"),
                "QC summary hides infinite/missing Merqury QV status")
        finite_qv += int(math.isfinite(qv))
        infinite_qv += int(math.isinf(qv))
        mq = project / "Assemblies" / sample / "QC/merqury" / stage / assembler
        expected_sources.update((mq / name).resolve() for name in ["merqury.qv", "completeness.stats", "inputs.json"])
        if stage == "raw":
            require(row["mapping_status"] == row["tidk_status"] == row["heterozygosity_status"] == "NOT_APPLICABLE_RAW",
                    "QC summary raw rows silently imply final-only QC")
            for field in ["properly_paired_reads_percent", "mapping_report", "tidk_best_motif", "tidk_detected_ends", "tidk_possible_ends",
                          "genotype_model", "callable_bases", "nuclear_acgt_bases", "callable_acgt_percent", "heterozygous_snps", "heterozygosity_percent"]:
                require(row[field] == "NA", f"Unexpected/missingness-hidden raw QC field: {field}")
            continue
        mapping = evidence["mapping"][assembler]
        require(mapping is not None and row["mapping_status"] == "MEASURED" and
                Path(row["mapping_report"]).resolve() == Path(mapping["flagstat"]).resolve(), "QC summary mapping identity/status mismatch")
        close_number(row["properly_paired_reads_percent"], mapping["properly_paired_percent"], "QC summary proper pairing")
        paired_rows += 1
        expected_sources.add(Path(mapping["flagstat"]).resolve())
        if truthy(config.get("run_tidk", False)):
            tidk = evidence["tidk"][assembler]
            require(tidk is not None, "Cannot validate QC summary without native TIDK evidence")
            summary = tidk["summary"]
            require(row["tidk_status"] == summary["screen_status"] and row["tidk_best_motif"] == summary["best_motif"],
                    "QC summary TIDK status/motif mismatch")
            close_number(row["tidk_detected_ends"], summary["detected_ends"], "QC summary TIDK detections")
            close_number(row["tidk_possible_ends"], summary["possible_sequence_ends"], "QC summary TIDK possible ends")
            tidk_testable += int(int(summary["possible_sequence_ends"]) > 0)
            folder = project / "Assemblies" / sample / "QC/tidk" / assembler
            expected_sources.update((folder / name).resolve() for name in ["assembly_summary.tsv", "run_config.json", "complete.json"])
        else:
            require(row["tidk_status"] == "NOT_REQUESTED" and all(row[field] == "NA" for field in
                    ["tidk_best_motif", "tidk_detected_ends", "tidk_possible_ends"]), "QC summary hides unrequested TIDK missingness")
        if truthy(config.get("run_heterozygosity", False)):
            het = evidence["heterozygosity"][assembler]
            require(het is not None, "Cannot validate QC summary without native heterozygosity evidence")
            summary = het["summary"]
            require(row["heterozygosity_status"] == summary["status"] and row["genotype_model"] == summary["genotype_model"],
                    "QC summary heterozygosity status/model mismatch")
            for field in ["callable_bases", "nuclear_acgt_bases", "heterozygous_snps", "heterozygosity_percent"]:
                close_number(row[field], summary[field], "QC summary " + field)
            close_number(row["callable_acgt_percent"], 100 * summary["callable_bases"] / summary["nuclear_acgt_bases"]
                         if summary["nuclear_acgt_bases"] else None, "QC summary callable fraction")
            het_testable += int(summary["heterozygosity_percent"] is not None)
            folder = project / "Assemblies" / sample / "QC/heterozygosity" / assembler
            expected_sources.update((folder / name).resolve() for name in ["summary.json", "receipt.json"])
        else:
            require(row["heterozygosity_status"] == "NOT_REQUESTED" and all(row[field] == "NA" for field in
                    ["genotype_model", "callable_bases", "nuclear_acgt_bases", "callable_acgt_percent", "heterozygous_snps", "heterozygosity_percent"]),
                    "QC summary hides unrequested heterozygosity missingness")
    recorded_sources = {Path(path).resolve(): digest for path, digest in provenance["source_sha256"].items()}
    require(set(recorded_sources) == expected_sources, "QC summary source manifest does not cover exactly its native evidence")
    for path, digest in recorded_sources.items():
        require(file_sha256(path) == digest, f"QC summary source checksum mismatch: {path}")
    plot_counts = {"unique_samples": len({row["sample"] for row in rows}),
                   "all_assembly_rows": len(rows), "final_assembly_rows": len(rows) // 2,
                   "finite_qv_rows": finite_qv, "infinite_qv_rows": infinite_qv,
                   "proper_pairing_rows": paired_rows, "tidk_testable_rows": tidk_testable,
                   "heterozygosity_testable_rows": het_testable}
    require(provenance["plot_counts"] == plot_counts and provenance["plotting_versions"]["backend"] == "Agg",
            "QC summary plot counts hide missing/infinite/untestable evidence")
    images = {}
    png = nonempty(directory / "assembly_qc.png")
    require(png.stat().st_size <= 100 * 1024 * 1024, "QC summary PNG has unbounded file size")
    images["png"] = png_info(png)
    require(100 <= images["png"]["width"] <= 10000 and 100 <= images["png"]["height"] <= 10000,
            "QC summary PNG canvas exceeds bounded dimensions")
    pdf = nonempty(directory / "assembly_qc.pdf")
    require(pdf.stat().st_size <= 50 * 1024 * 1024, "QC summary PDF has unbounded file size")
    content = pdf.read_bytes()
    require(content.startswith(b"%PDF-") and b"%%EOF" in content[-1024:], "Invalid QC summary PDF")
    boxes = re.findall(rb"/MediaBox\s*\[\s*([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s*\]", content)
    require(boxes and len(re.findall(rb"/Type\s*/Page\b", content)) == 1,
            "QC summary PDF must have one bounded page")
    for box in boxes:
        x1, y1, x2, y2 = map(float, box)
        require(all(math.isfinite(x) for x in [x1, y1, x2, y2]) and 100 <= x2 - x1 <= 2000 and 100 <= y2 - y1 <= 2000,
                "QC summary PDF canvas is out of bounds")
    images["pdf"] = {"path": str(pdf), "sha256": file_sha256(pdf), "bytes": pdf.stat().st_size}
    svg = nonempty(directory / "assembly_qc.svg")
    require(svg.stat().st_size <= 50 * 1024 * 1024, "QC summary SVG has unbounded file size")
    root = ET.parse(svg).getroot()
    require(root.tag.split("}")[-1] == "svg", "QC summary SVG has an invalid root")
    viewbox = list(map(float, root.attrib["viewBox"].split()))
    require(len(viewbox) == 4 and all(math.isfinite(x) for x in viewbox) and
            100 <= viewbox[2] <= 2000 and 100 <= viewbox[3] <= 2000, "QC summary SVG canvas is out of bounds")
    images["svg"] = {"path": str(svg), "sha256": file_sha256(svg), "viewBox": viewbox}
    return {"row_count": len(rows), "identities": [list(identity) for identity in identities],
            "native_sources_checked": len(recorded_sources), "plot_counts": plot_counts, "images": images}


class Validator:
    def __init__(self):
        self.checks = []
        self.errors = []

    def attempt(self, name, operation):
        try:
            value = operation()
        except Exception as error:
            self.errors.append({"check": name, "error": str(error), "type": type(error).__name__})
            self.checks.append({"check": name, "status": "failed"})
            return None
        self.checks.append({"check": name, "status": "passed"})
        return value


def validate(config, samtools):
    project = Path(config["project_dir"]).resolve()
    with Path(config["sample_sheet"]).open() as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    require(samples, "Sample sheet has no samples")
    names = [row["sample"].strip() for row in samples]
    require(len(names) == len(set(names)), "Duplicate sample identifiers in validation input")
    k = int(config.get("merqury_k", 21))
    require(k > 0, "Invalid k-mer size")
    validator, results = Validator(), {}
    for row in samples:
        sample, assemblers = row["sample"].strip(), selected_assemblers(row)
        require(assemblers, f"No enabled assembler for {sample}")
        branch = project / "Assemblies" / sample
        qc = branch / "QC"
        result = {"assemblers": assemblers, "ragtag": truthy(row.get("ragtag", False)),
                  "fasta": {}, "mitochondrial_integration": {}, "merqury": {},
                  "chromeister": {}, "busco": {}, "mapping": {}, "tidk": {}, "heterozygosity": {}}
        results[sample] = result
        mitochondrial = validator.attempt(f"{sample}/OATK_FASTA", lambda: fasta(branch / "OATK" / f"{sample}.mito.ctg.fasta"))
        validator.attempt(f"{sample}/OATK_marker", lambda: marker(branch / "OATK/.done", [branch / "OATK" / f"{sample}.mito.ctg.fasta"]))
        validator.attempt(f"{sample}/NanoPlot", lambda: marker(qc / "nanoplot/nanoplot.done", [qc / "nanoplot/NanoStats.txt"]))
        ranked = branch / "TELOSEARCH" / f"TeloSearchLR.k4.K20.t1000.{sample}.hifireads.rankedRepeatsTable.txt"
        validator.attempt(f"{sample}/TeloSearchLR", lambda: (
            nonempty(branch / "TELOSEARCH" / f"{sample}.hifireads.fasta"), nonempty(ranked),
            marker(branch / "TELOSEARCH/telosearch.done", [ranked]),
            require("Script finished." in nonempty(branch / "logs/telosearchlr.log").read_text(),
                    "TeloSearchLR lacks its completion message")))
        summary = qc / "GenomeScope2" / f"{sample}.genomescope2.summary.tsv"
        result["genomescope"] = validator.attempt(f"{sample}/GenomeScope", lambda: genomescope_check(summary))
        validator.attempt(f"{sample}/GenomeScope_marker", lambda: marker(qc / "GenomeScope2/ploidy.done", [summary, qc / "GenomeScope2" / f"{sample}.histo"]))
        preferred_assembler = assemblers[0]
        result["nquire"] = validator.attempt(f"{sample}/nQuire", lambda: nquire_check(
            qc / "nQuire", sample, qc / "Mapping" / f"purged-{preferred_assembler}.bam"))
        if "hifiasm" in assemblers:
            gff = branch / "FILTERED" / f"{sample}-barrnap-hifiasm.gff3"
            result["barrnap"] = validator.attempt(f"{sample}/Barrnap", lambda: gff_check(gff))
            validator.attempt(f"{sample}/Barrnap_marker", lambda: marker(branch / "FILTERED/barrnap.done", [gff]))
        for assembler in assemblers:
            raw_path = branch / (f"HIFIASM/{sample}.fasta" if assembler == "hifiasm" else "LJA/assembly.fasta")
            final_path = project / f"{sample}_FINAL_ASSEMBLIES" / assembler / f"{sample}.purged.final.fasta"
            for stage, path in [("raw", raw_path), ("final", final_path)]:
                key = f"{stage}/{assembler}"
                sequences = validator.attempt(f"{sample}/{key}/FASTA", lambda path=path: fasta(path))
                stats = validator.attempt(f"{sample}/{key}/FASTA_stats", lambda path=path, seq=sequences: fasta_stats(path, seq, k)) if sequences else None
                result["fasta"][key] = stats
                if stage == "raw":
                    done = branch / ("HIFIASM/hifiasm.done" if assembler == "hifiasm" else "LJA/lja.done")
                    validator.attempt(f"{sample}/{key}/assembly_marker", lambda done=done, path=path: marker(done, [path]))
                if stage == "final" and sequences and mitochondrial:
                    result["mitochondrial_integration"][assembler] = validator.attempt(
                        f"{sample}/{assembler}/mitochondrial_integration",
                        lambda seq=sequences: integration_check(project, sample, assembler, seq, mitochondrial, result["ragtag"], config))
                if stage == "final" and truthy(config.get("run_tidk", False)):
                    result["tidk"][assembler] = validator.attempt(f"{sample}/{assembler}/TIDK",
                        lambda path=path, sequences=sequences: tidk_check(qc / "tidk" / assembler, path, sequences or fasta(path), config))
                if stage == "final" and truthy(config.get("run_heterozygosity", False)):
                    result["heterozygosity"][assembler] = validator.attempt(f"{sample}/{assembler}/heterozygosity",
                        lambda path=path, sequences=sequences: heterozygosity_check(qc / "heterozygosity" / assembler,
                            project, sample, assembler, path, sequences or fasta(path), config, samtools))
                directory = qc / "merqury" / stage / assembler
                if stats:
                    result["merqury"][key] = validator.attempt(f"{sample}/{key}/Merqury",
                        lambda directory=directory, path=path, stats=stats: merqury_check(directory, path, stats, sample, stage, assembler, project, k))
                else:
                    validator.attempt(f"{sample}/{key}/Merqury", lambda directory=directory: nonempty(directory / "merqury.qv"))
                prefix = f"{sample}-{assembler}" if stage == "raw" else f"{sample}-final_{assembler}"
                chrome = branch / "CHROMEISTER" if stage == "raw" else final_path.parent / "CHROMEISTER"
                png = chrome / f"{prefix}.mat.filt.png"
                result["chromeister"][key] = validator.attempt(f"{sample}/{key}/Chromeister_PNG", lambda png=png: png_info(png))
                done = chrome / (f".{assembler}.done" if stage == "raw" else f".final_{assembler}.done")
                validator.attempt(f"{sample}/{key}/Chromeister_marker", lambda done=done, chrome=chrome, prefix=prefix, png=png: marker(done, [png, chrome / f"{prefix}.mat", chrome / f"{prefix}_score.txt"]))
                busco_name = assembler if stage == "raw" else f"purged-{assembler}"
                directory = qc / f"BUSCO-{busco_name}" / f"{sample}-{busco_name}"
                result["busco"][key] = validator.attempt(f"{sample}/{key}/BUSCO", lambda directory=directory: busco_check(directory))
                done = qc / f"BUSCO-{busco_name}" / (f"busco_{assembler}.done" if stage == "raw" else f"busco_purged_{assembler}.done")
                validator.attempt(f"{sample}/{key}/BUSCO_marker", lambda done=done: marker(done))
            purge = branch / ("PURGE_DUPS" if assembler == "hifiasm" else "PURGE_DUPS_LJA")
            done = purge / ("purge_hifiasm.done" if assembler == "hifiasm" else "purge_lja.done")
            validator.attempt(f"{sample}/{assembler}/purge_marker", lambda done=done, purge=purge, assembler=assembler: marker(done, [purge / f"purged.{assembler}.fa", purge / f"purged-{assembler}.sorted.fa"]))
            result["mapping"][assembler] = validator.attempt(f"{sample}/{assembler}/mapping", lambda assembler=assembler: mapping_check(qc / "Mapping", assembler, samtools))
        result["quast"] = validator.attempt(f"{sample}/QUAST", lambda: quast_check(qc / "QUAST/report.tsv", assemblers))
        validator.attempt(f"{sample}/QUAST_marker", lambda: marker(qc / "QUAST/quast.done", [qc / "QUAST/report.tsv"]))
        validator.attempt(f"{sample}/mapping_marker", lambda: marker(qc / "Mapping/mapping_qc.done", [qc / "Mapping" / f"purged-{a}.flagstat.txt" for a in assemblers]))
        for stage in ["raw", "final"]:
            validator.attempt(f"{sample}/Merqury_{stage}_marker", lambda stage=stage: marker(qc / "merqury" / stage / f"{stage}.done", [qc / "merqury" / stage / a / "merqury.qv" for a in assemblers]))
        validator.attempt(f"{sample}/Merqury_marker", lambda: marker(qc / "merqury/merqury.done", [qc / "merqury/raw/raw.done", qc / "merqury/final/final.done"]))
    multiqc = None
    if truthy(config.get("run_multiqc", False)):
        html = project / "Assemblies/multiqc/multiqc_report.html"
        validator.attempt("MultiQC", lambda: (
            marker(html.parent / "multiqc.done", [html]),
            require("<html" in html.read_text(errors="replace").lower(), "MultiQC output is not HTML")))
        multiqc = validator.attempt("MultiQC_report_identity_and_provenance",
                                   lambda: multiqc_data_check(html.parent / "multiqc_data", project, samples))
    qc_summary = None
    if truthy(config.get("run_qc_summary", False)):
        qc_summary = validator.attempt("QC_summary_native_metrics_provenance_and_plots",
            lambda: qc_summary_check(project / "Assemblies/qc_summary", project, samples, results, config))
    return {"status": "passed" if not validator.errors else "failed", "project_dir": str(project),
            "sample_count": len(samples), "sample_ids": names, "kmer_k": k,
            "samtools": samtools, "check_count": len(validator.checks),
            "passed_check_count": sum(check["status"] == "passed" for check in validator.checks),
            "checks": validator.checks, "errors": validator.errors, "samples": results, "multiqc": multiqc, "qc_summary": qc_summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samtools", default=shutil.which("samtools"), help="Optional samtools executable for BAM quickcheck")
    args = parser.parse_args()
    try:
        config = json.loads(args.config.read_text())
        result = validate(config, args.samtools)
    except Exception as error:
        result = {"status": "failed", "errors": [{"check": "configuration", "error": str(error), "type": type(error).__name__}]}
    result.update(validated_at_utc=datetime.now(timezone.utc).isoformat(), config=str(args.config.resolve()),
                  validator=str(Path(__file__).resolve()), analyses_rerun=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Validation {result['status']}: {result.get('passed_check_count', 0)}/{result.get('check_count', 0)} checks passed; {len(result['errors'])} issue(s). Report: {args.output}")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
