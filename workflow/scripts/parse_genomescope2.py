#!/usr/bin/env python3
"""Extract GenomeScope summary fit scores without asserting biological ploidy.

The console field called 'model fit' is adups in GenomeScope 2.0, not a fit
score. Only summary.txt's explicit Model Fit percentages are used here.
"""
import argparse
import csv
import math
from pathlib import Path
import re

NUMBER = r"[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?"
HEADER = re.compile(r"(?:GenomeScope analyzing .*?p=|\[INFO\] GenomeScope: p=)(\d+)\s+k=(\d+).*?outdir=(\S+)")
FIT = re.compile(rf"^Model Fit\s+({NUMBER})%\s+({NUMBER})%\s*$", re.MULTILINE)
LENGTH = re.compile(r"^Genome Haploid Length\s+([0-9,]+) bp\s+([0-9,]+) bp", re.MULTILINE)
FIELDS = ["ploidy", "k", "outdir", "genome_len_min", "genome_len_max",
          "model_fit_all_percent", "model_fit_full_percent", "source_summary",
          "source_log", "fit_rank", "best_by_model_fit", "selection_status"]


def parse_one(log_path):
    log_path = Path(log_path).resolve()
    rows = []
    seen = set()
    for line in log_path.read_text(errors="replace").splitlines():
        match = HEADER.search(line)
        if not match:
            continue
        ploidy, k, folder = match.groups()
        directory = Path(folder)
        if not directory.is_absolute():
            directory = log_path.parent / directory
        summary = directory / "summary.txt"
        if summary in seen or not summary.is_file():
            continue
        seen.add(summary)
        content = summary.read_text(errors="replace")
        fit, length = FIT.search(content), LENGTH.search(content)
        if not fit or not length:
            continue
        all_fit, full_fit = map(float, fit.groups())
        lengths = [int(x.replace(",", "")) for x in length.groups()]
        if not all(math.isfinite(x) and x <= 100 for x in (all_fit, full_fit)) or min(lengths) <= 0:
            continue
        rows.append({"ploidy": int(ploidy), "k": int(k), "outdir": str(directory),
                     "genome_len_min": min(lengths), "genome_len_max": max(lengths),
                     "model_fit_all_percent": all_fit, "model_fit_full_percent": full_fit,
                     "source_summary": str(summary), "source_log": str(log_path)})
    return rows


def summarize(logs):
    rows = [row for log in logs for row in parse_one(log)]
    if not rows:
        raise ValueError("No converged GenomeScope summary with valid Model Fit and genome length")
    unique = {row["source_summary"]: row for row in rows}
    rows = list(unique.values())
    scores = sorted({row["model_fit_all_percent"] for row in rows}, reverse=True)
    for row in rows:
        rank = scores.index(row["model_fit_all_percent"]) + 1
        row.update(fit_rank=rank, best_by_model_fit="YES" if rank == 1 else "",
                   selection_status="exploratory_fit_ranking_not_validated_ploidy")
    return sorted(rows, key=lambda row: (row["ploidy"], row["source_summary"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("logs", type=Path, nargs="+")
    args = parser.parse_args()
    rows = summarize(args.logs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
