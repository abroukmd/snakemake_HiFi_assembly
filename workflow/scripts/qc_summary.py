#!/usr/bin/env python3
"""Collect traceable sample/assembler/stage QC and draw a descriptive overview."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys


FIELDS = [
    "sample", "assembler", "stage", "assembly", "assembly_sha256",
    "sequence_count", "total_bp", "sequence_n50_bp", "longest_sequence_bp",
    "gap_n_bp", "gc_percent_acgt", "busco_lineage", "busco_markers",
    "busco_complete_percent", "busco_single_percent", "busco_duplicated_percent",
    "busco_fragmented_percent", "busco_missing_percent", "busco_report",
    "merqury_qv", "merqury_completeness_percent", "merqury_k", "merqury_status",
    "properly_paired_reads_percent", "mapping_status", "mapping_report",
    "tidk_status", "tidk_best_motif", "tidk_detected_ends", "tidk_possible_ends",
    "heterozygosity_status", "genotype_model", "callable_bases",
    "nuclear_acgt_bases", "callable_acgt_percent", "heterozygous_snps",
    "heterozygosity_percent",
]


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def fasta_statistics(path, k=21):
    """Count FASTA records, including short and organellar records, directly."""
    lengths, seen, counts = [], set(), Counter()
    current = None
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    if current == 0:
                        raise ValueError(f"Empty sequence in {path}")
                    lengths.append(current)
                identifier = line[1:].split()[0] if line[1:].split() else ""
                if not identifier or identifier in seen:
                    raise ValueError(f"Empty/duplicate FASTA identifier in {path}:{number}")
                seen.add(identifier)
                current = 0
            else:
                sequence = line.upper()
                if current is None or set(sequence) - set("ACGTRYSWKMBDHVNU"):
                    raise ValueError(f"Invalid nucleotide FASTA in {path}:{number}")
                current += len(sequence)
                counts.update(sequence)
    if current is None or current == 0:
        raise ValueError(f"Missing or empty FASTA records in {path}")
    lengths.append(current)
    total, cumulative, n50 = sum(lengths), 0, None
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= total:
            n50 = length
            break
    acgt = sum(counts[base] for base in "ACGT")
    return {"assembly_sha256": digest(path), "sequence_count": len(lengths),
            "total_bp": total, "sequence_n50_bp": n50,
            "longest_sequence_bp": max(lengths), "gap_n_bp": counts["N"],
            "gc_percent_acgt": 100 * (counts["G"] + counts["C"]) / acgt if acgt else None,
            "maximum_kmer_windows": sum(max(0, length - k + 1) for length in lengths)}


def busco_summary(directory, assembly):
    directory = Path(directory)
    specific = sorted(directory.glob("short_summary.specific.*.txt"))
    candidates = specific or sorted(directory.glob("short_summary.generic.*.txt"))
    if len(candidates) != 1:
        raise ValueError(f"Expected one selected BUSCO summary in {directory}: {candidates}")
    report = candidates[0]
    content = report.read_text()
    lineage = re.search(r"The lineage dataset is:\s*(\S+)", content)
    scores = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", content)
    source = re.search(r"Summarized benchmarking in BUSCO notation for file (.+)", content)
    if not lineage or not scores or not source:
        raise ValueError(f"Incomplete BUSCO summary: {report}")
    if Path(source.group(1).strip()).resolve() != Path(assembly).resolve():
        raise ValueError(f"BUSCO describes a different FASTA: {report}")
    values = list(map(float, scores.groups()[:5]))
    complete, single, duplicated, fragmented, missing = values
    if (any(not math.isfinite(v) or not 0 <= v <= 100 for v in values)
            or abs(complete - single - duplicated) > .31
            or abs(complete + fragmented + missing - 100) > .31):
        raise ValueError(f"Inconsistent BUSCO percentages: {report}")
    return {"busco_lineage": lineage.group(1), "busco_markers": int(scores.group(6)),
            **dict(zip(["busco_complete_percent", "busco_single_percent", "busco_duplicated_percent",
                        "busco_fragmented_percent", "busco_missing_percent"], values)),
            "busco_report": str(report)}


def merqury_summary(directory, expected_label, assembly, maximum_windows):
    directory = Path(directory)
    qv = [line.split() for line in (directory / "merqury.qv").read_text().splitlines() if line.strip()]
    complete = [line.split() for line in (directory / "completeness.stats").read_text().splitlines() if line.strip()]
    if (len(qv) != 1 or len(qv[0]) != 5 or qv[0][0] != expected_label
            or len(complete) != 1 or len(complete[0]) != 5 or complete[0][0] != expected_label):
        raise ValueError(f"Ambiguous or mislabelled single-assembly Merqury reports: {directory}")
    quality, percent = float(qv[0][3]), float(complete[0][4])
    errors, total = int(qv[0][1]), int(qv[0][2])
    recovered, solid = int(complete[0][2]), int(complete[0][3])
    if (not 0 <= errors <= total <= maximum_windows or total == 0 or math.isnan(quality)
            or quality < 0 or (math.isinf(quality) and errors != 0)
            or not 0 <= recovered <= solid or solid == 0
            or not math.isfinite(percent) or not 0 <= percent <= 100):
        raise ValueError(f"Invalid Merqury metrics: {directory}")
    provenance = json.loads((directory / "inputs.json").read_text())
    if (provenance["label"] != expected_label
            or Path(provenance["assembly"]).resolve() != Path(assembly).resolve()
            or provenance.get("assembly_sha256") != digest(assembly)):
        raise ValueError(f"Merqury provenance does not match its assembly: {directory}")
    return {"merqury_qv": quality, "merqury_completeness_percent": percent,
            "merqury_k": int(provenance["k"]),
            "merqury_status": "INFINITE_QV_ZERO_ERROR_KMERS" if math.isinf(quality) else "MEASURED"}


def mapping_summary(path):
    content = Path(path).read_text()
    proper = re.search(r"^(\d+) \+ \d+ properly paired", content, re.MULTILINE)
    paired = re.search(r"^(\d+) \+ \d+ paired in sequencing", content, re.MULTILINE)
    if not proper or not paired:
        raise ValueError(f"Missing proper-pair denominator in {path}")
    numerator, denominator = int(proper.group(1)), int(paired.group(1))
    if not 0 <= numerator <= denominator:
        raise ValueError(f"Inconsistent proper-pair counts in {path}")
    return {"properly_paired_reads_percent": 100 * numerator / denominator if denominator else None,
            "mapping_status": "MEASURED" if denominator else "NO_PAIRED_READS",
            "mapping_report": str(path)}


def collect(records, table, provenance, source_script=None):
    if not records:
        raise ValueError("No selected assembly records to summarize")
    rows, sources, identities = [], {}, set()
    for record in records:
        identity = tuple(record[field] for field in ("sample", "assembler", "stage"))
        if identity in identities:
            raise ValueError(f"Duplicate QC identity: {identity}")
        identities.add(identity)
        assembly = Path(record["assembly"])
        mq_inputs = json.loads((Path(record["merqury"]) / "inputs.json").read_text())
        statistics = fasta_statistics(assembly, int(mq_inputs["k"]))
        row = {field: None for field in FIELDS}
        row.update({field: record[field] for field in ("sample", "assembler", "stage")})
        row.update(assembly=str(assembly), **{k: v for k, v in statistics.items() if k in FIELDS})
        row.update(busco_summary(record["busco"], assembly))
        label = f"{record['sample']}.{record['stage']}.{record['assembler']}"
        row.update(merqury_summary(record["merqury"], label, assembly,
                                   statistics["maximum_kmer_windows"]))
        # Native Merqury labels are sample.stage.assembler, not table-key order.
        row["mapping_status"] = "NOT_APPLICABLE_RAW"
        row["tidk_status"] = row["heterozygosity_status"] = "NOT_APPLICABLE_RAW"
        paths = [assembly, Path(row["busco_report"]), Path(record["merqury"]) / "merqury.qv",
                 Path(record["merqury"]) / "completeness.stats", Path(record["merqury"]) / "inputs.json"]
        if record["stage"] == "final":
            row.update(mapping_summary(record["mapping"]))
            paths.append(Path(record["mapping"]))
            row["tidk_status"] = row["heterozygosity_status"] = "NOT_REQUESTED"
            if record.get("tidk"):
                tidk_path = Path(record["tidk"])
                config_path, receipt_path = tidk_path.parent / "run_config.json", tidk_path.parent / "complete.json"
                config, receipt = json.loads(config_path.read_text()), json.loads(receipt_path.read_text())
                if (Path(config["input_path"]).resolve() != assembly.resolve()
                        or config["input_sha256"] != row["assembly_sha256"]
                        or receipt["output_sha256"][tidk_path.name] != digest(tidk_path)
                        or receipt["output_sha256"][config_path.name] != digest(config_path)):
                    raise ValueError(f"TIDK provenance does not match its assembly/artifacts: {tidk_path}")
                with tidk_path.open() as handle:
                    tidk_rows = list(csv.DictReader(handle, delimiter="\t"))
                if len(tidk_rows) != 1:
                    raise ValueError(f"Expected one assembly TIDK summary: {record['tidk']}")
                tidk = tidk_rows[0]
                row.update(tidk_status=tidk["screen_status"], tidk_best_motif=tidk["best_motif"],
                           tidk_detected_ends=int(tidk["detected_ends"]),
                           tidk_possible_ends=int(tidk["possible_sequence_ends"]))
                if (int(tidk["input_sequences"]) != row["sequence_count"]
                        or not 0 <= row["tidk_detected_ends"] <= row["tidk_possible_ends"] <= 2 * row["sequence_count"]):
                    raise ValueError(f"Inconsistent TIDK sequence/end counts: {tidk_path}")
                paths.extend([tidk_path, config_path, receipt_path])
            if record.get("heterozygosity"):
                het_path = Path(record["heterozygosity"])
                receipt_path = het_path.parent / "receipt.json"
                het, receipt = json.loads(het_path.read_text()), json.loads(receipt_path.read_text())
                source = receipt["inputs"]["reference"]
                if (het["sample"] != record["sample"] or het["assembler"] != record["assembler"]
                        or receipt["state"] != "completed"
                        or Path(source["path"]).resolve() != assembly.resolve()
                        or source["sha256"] != row["assembly_sha256"]
                        or receipt["outputs"][het_path.name]["sha256"] != digest(het_path)):
                    raise ValueError(f"Heterozygosity provenance does not match its assembly/artifacts: {het_path}")
                row.update(heterozygosity_status=het["status"], genotype_model=het["genotype_model"])
                for field in ("callable_bases", "nuclear_acgt_bases", "heterozygous_snps", "heterozygosity_percent"):
                    row[field] = het.get(field)
                total, callable_bases = row["nuclear_acgt_bases"], row["callable_bases"]
                sites, percent = row["heterozygous_snps"], row["heterozygosity_percent"]
                if not 0 <= sites <= callable_bases <= total:
                    raise ValueError(f"Invalid callable-site/SNP counts: {het_path}")
                if ((callable_bases == 0 and percent is not None)
                        or (callable_bases > 0 and (percent is None or not math.isfinite(percent)
                            or not math.isclose(percent, 100 * sites / callable_bases, rel_tol=1e-8, abs_tol=1e-10)))):
                    raise ValueError(f"Inconsistent callable-site heterozygosity: {het_path}")
                row["callable_acgt_percent"] = 100 * callable_bases / total if total and callable_bases is not None else None
                paths.extend([het_path, receipt_path])
        rows.append(row)
        sources.update({str(path): digest(path) for path in paths})
    Path(table).parent.mkdir(parents=True, exist_ok=True)
    with Path(table).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if v is None else "+inf" if isinstance(v, float) and math.isinf(v) else v
                             for k, v in row.items()})
    Path(provenance).write_text(json.dumps({
        "identities": [list(identity) for identity in sorted(identities)], "source_sha256": sources,
        "script_sha256": digest(source_script or __file__), "python": sys.version,
        "definitions": {"unit": "sample x assembler x stage; assemblies of a strain are not biological replicates",
                        "sequence_count": "all FASTA records, including short and organellar sequences",
                        "proper_pairing": "QC-passed properly-paired / QC-passed paired-in-sequencing reads",
                        "heterozygosity": "depth-callable diploid-model estimate; denominator has no missing-GT/GQ gate; see per-assembly filters and receipt",
                        "telomeres": "exact-repeat candidate screening; not a T2T or chromosome-count claim",
                        "busco": "selected specific summary, or generic summary when no specific summary exists"}},
        indent=2) + "\n")
    return rows


def plot(rows, pdf, svg, png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    final = [row for row in rows if row["stage"] == "final"]
    extra = any(row["tidk_status"] != "NOT_REQUESTED" or row["heterozygosity_status"] != "NOT_REQUESTED" for row in final)
    fig, axes = plt.subplots(2, 3 if extra else 2, figsize=(17 if extra else 13, 10), squeeze=False)
    axes = list(axes.flat)
    colors = {"hifiasm": "#2468a2", "lja": "#c16b23"}
    markers = {"raw": "o", "final": "^"}
    for row in rows:
        style = dict(color=colors.get(row["assembler"], "#666666"), marker=markers[row["stage"]], alpha=.8)
        axes[0].scatter(row["total_bp"] / 1e6, row["sequence_n50_bp"] / 1e6, **style)
        qv = row["merqury_qv"]
        if math.isfinite(qv):
            axes[1].scatter(row["merqury_completeness_percent"], qv, **style)
        else:
            axes[1].scatter(row["merqury_completeness_percent"], 1.03,
                            transform=axes[1].get_xaxis_transform(), clip_on=False, **style)
            axes[1].annotate("∞", (row["merqury_completeness_percent"], 1.04),
                             xycoords=axes[1].get_xaxis_transform(), ha="center")
    axes[0].set(xlabel="Assembly length (Mbp)", ylabel="FASTA sequence N50 (Mbp)", title="A  Contiguity")
    axes[1].set(xlabel="Merqury completeness (%)", ylabel="Consensus QV", title="B  K-mer evidence")
    # Mixed data/axes transforms do not initialize data limits. Without an
    # explicit physical range, an all-infinite cohort can expand tight exports
    # to hundreds of thousands of pixels instead of showing a small top strip.
    axes[1].set_xlim(0, 100)
    finite_qvs = [row["merqury_qv"] for row in rows if math.isfinite(row["merqury_qv"])]
    if finite_qvs:
        lower, upper = min(finite_qvs), max(finite_qvs)
        padding = max(1, .08 * (upper - lower))
        axes[1].set_ylim(max(0, lower - padding), upper + padding)
    else:
        axes[1].set_ylim(0, 1)
        axes[1].set_yticks([])
        axes[1].text(.5, .5, "No finite QV values", transform=axes[1].transAxes, ha="center")
    if any(math.isinf(row["merqury_qv"]) for row in rows):
        axes[1].text(.02, .02, "∞ shown above axes; no finite QV substituted", transform=axes[1].transAxes, fontsize=8)
    ax = axes[2]
    if len(final) <= 24:
        left = [0.] * len(final)
        for key, name, color in [("single", "Single", "#228833"), ("duplicated", "Duplicated", "#77aa66"),
                                 ("fragmented", "Fragmented", "#ddaa33"), ("missing", "Missing", "#bb5566")]:
            values = [row[f"busco_{key}_percent"] for row in final]
            ax.barh(range(len(final)), values, left=left, label=name, color=color)
            left = [a + b for a, b in zip(left, values)]
        ax.set_yticks(range(len(final)), [f"{r['sample']} {r['assembler']}\n{r['busco_lineage']} n={r['busco_markers']}" for r in final], fontsize=7)
        ax.set(xlabel="BUSCO components (%)", xlim=(0, 100))
        ax.legend(fontsize=7, ncol=2)
    else:
        lineages = sorted({row["busco_lineage"] for row in final})
        for index, lineage in enumerate(lineages):
            selected = [row for row in final if row["busco_lineage"] == lineage]
            for asm in sorted({row["assembler"] for row in selected}):
                values = [row["busco_complete_percent"] for row in selected if row["assembler"] == asm]
                shift = -.12 if asm == "hifiasm" else .12
                ax.scatter([index + shift] * len(values), values, color=colors.get(asm, "#666666"), alpha=.5, s=15)
        ax.set_xticks(range(len(lineages)), lineages, rotation=60, ha="right", fontsize=7)
        ax.set(ylabel="Complete BUSCO (%)")
    ax.set_title("C  Final BUSCO (retain lineage context)")
    for row in final:
        color = colors.get(row["assembler"], "#666666")
        if row["properly_paired_reads_percent"] is not None:
            axes[3].scatter(row["total_bp"] / 1e6, row["properly_paired_reads_percent"], color=color, alpha=.8)
        if extra:
            if row["tidk_possible_ends"]:
                axes[4].scatter(row["tidk_possible_ends"], row["tidk_detected_ends"], color=color, alpha=.8)
            if row["callable_acgt_percent"] is not None and row["heterozygosity_percent"] is not None:
                axes[5].scatter(row["callable_acgt_percent"], row["heterozygosity_percent"], color=color, alpha=.8)
    axes[3].set(xlabel="Final assembly length (Mbp)", ylabel="Properly paired reads (%)", title="D  Illumina alignment")
    if extra:
        axes[4].set(xlabel="Possible retained sequence ends", ylabel="Candidate-positive ends", title="E  Exact-repeat end screening")
        axes[5].set(xlabel="Depth-callable nuclear ACGT (%)", ylabel="Heterozygous SNPs / callable bp (%)", title="F  Conditional diploid-model estimate")
        if not any(row["tidk_possible_ends"] for row in final):
            message = ("TIDK not requested" if all(row["tidk_status"] == "NOT_REQUESTED" for row in final)
                       else "No retained sequences to screen")
            axes[4].text(.5, .5, message, transform=axes[4].transAxes, ha="center")
        if not any(row["heterozygosity_percent"] is not None for row in final):
            message = ("Heterozygosity not requested" if all(row["heterozygosity_status"] == "NOT_REQUESTED" for row in final)
                       else "Untestable with the current callable mask")
            axes[5].text(.5, .5, message, transform=axes[5].transAxes, ha="center")
    if not any(row["properly_paired_reads_percent"] is not None for row in final):
        axes[3].text(.5, .5, "No paired-read denominator", transform=axes[3].transAxes, ha="center")
    legend = [Line2D([], [], marker="o", ls="", color=colors.get(asm, "#666666"), label=asm)
              for asm in sorted({row["assembler"] for row in rows})]
    legend += [Line2D([], [], marker=markers[stage], ls="", color="#555555", label=stage)
               for stage in sorted({row["stage"] for row in rows})]
    fig.legend(handles=legend, loc="upper center", ncol=4, frameon=False)
    fig.text(.5, .01, "Descriptive QC; lineage sets differ. Telomere candidates and diploid-model SNP estimates do not establish T2T status or cellular ploidy.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .95), h_pad=3, w_pad=3)
    for target in [pdf, svg, png]:
        # Keep a fixed physical canvas even for off-axis infinity annotations.
        fig.savefig(target, dpi=250)
    plt.close(fig)


def main():
    source_script = None
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        records = list(sm.params.records)
        outputs = sm.output
        source_script = sm.input.source_script
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--manifest", type=Path, required=True)
        parser.add_argument("--outdir", type=Path, required=True)
        args = parser.parse_args()
        records = json.loads(args.manifest.read_text())
        args.outdir.mkdir(parents=True, exist_ok=True)
        from types import SimpleNamespace
        outputs = SimpleNamespace(table=args.outdir / "assembly_qc.tsv", provenance=args.outdir / "provenance.json",
                                  pdf=args.outdir / "assembly_qc.pdf", svg=args.outdir / "assembly_qc.svg", png=args.outdir / "assembly_qc.png")
    # Keep the input label definition independent of display/table-key order.
    def execute():
        rows = collect(records, outputs.table, outputs.provenance, source_script)
        plot(rows, outputs.pdf, outputs.svg, outputs.png)
        import matplotlib
        provenance = json.loads(Path(outputs.provenance).read_text())
        provenance["plotting_versions"] = {"matplotlib": matplotlib.__version__, "backend": "Agg"}
        final = [row for row in rows if row["stage"] == "final"]
        provenance["plot_counts"] = {
            "unique_samples": len({row["sample"] for row in rows}),
            "all_assembly_rows": len(rows), "final_assembly_rows": len(final),
            "finite_qv_rows": sum(math.isfinite(row["merqury_qv"]) for row in rows),
            "infinite_qv_rows": sum(math.isinf(row["merqury_qv"]) for row in rows),
            "proper_pairing_rows": sum(row["properly_paired_reads_percent"] is not None for row in final),
            "tidk_testable_rows": sum(bool(row["tidk_possible_ends"]) for row in final),
            "heterozygosity_testable_rows": sum(row["heterozygosity_percent"] is not None for row in final),
        }
        Path(outputs.provenance).write_text(json.dumps(provenance, indent=2) + "\n")
    if "snakemake" in globals():
        import contextlib
        import traceback
        with Path(sm.log[0]).open("w") as handle, contextlib.redirect_stdout(handle), contextlib.redirect_stderr(handle):
            try:
                execute()
            except Exception:
                traceback.print_exc()
                raise
    else:
        execute()


if __name__ == "__main__":
    main()
