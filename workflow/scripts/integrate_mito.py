#!/usr/bin/env python3
"""Replace whole mitochondrial contigs without deleting nuclear flanks.

PAF queries are OATK mitochondrial records and targets are assembly records.
Union coverage combines split alignments, including circular origin rotations;
both query and target coverage must pass before an assembly record is removed.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    records, sequence, seen = [], [], set()
    header = None
    with Path(path).open() as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                if not header or not header.split():
                    raise ValueError(f"Empty FASTA header in {path}:{lineno}")
                name = header.split()[0]
                if name in seen:
                    raise ValueError(f"Duplicate FASTA ID {name!r} in {path}")
                seen.add(name)
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}:{lineno}")
                sequence.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records or any(not seq for _, seq in records):
        raise ValueError(f"Missing or empty FASTA sequence in {path}")
    return records


def covered_length(intervals):
    """Count the union of half-open intervals without double-counting overlaps."""
    total, right = 0, 0
    for start, end in sorted(intervals):
        total += max(0, end - max(start, right))
        right = max(right, end)
    return total


def mitochondrial_matches(paf, assemblies, mitochondria, min_identity, min_qcov, min_tcov):
    hits = defaultdict(lambda: {"query": [], "target": []})
    with Path(paf).open() as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF record at {paf}:{lineno}")
            query, target = fields[0], fields[5]
            if query not in mitochondria or target not in assemblies:
                raise ValueError(f"PAF ID absent from input FASTAs at {paf}:{lineno}")
            qlen, qstart, qend = map(int, fields[1:4])
            tlen, tstart, tend = map(int, fields[6:9])
            matches, alnlen = map(int, fields[9:11])
            if (qlen != len(mitochondria[query]) or tlen != len(assemblies[target])
                    or not 0 <= qstart < qend <= qlen
                    or not 0 <= tstart < tend <= tlen
                    or not 0 <= matches <= alnlen or alnlen <= 0):
                raise ValueError(f"Invalid PAF coordinates or lengths at {paf}:{lineno}")
            if matches / alnlen >= min_identity:
                hits[(target, query)]["query"].append((qstart, qend))
                hits[(target, query)]["target"].append((tstart, tend))
    removed, report = set(), []
    for (target, query), intervals in sorted(hits.items()):
        qcov = covered_length(intervals["query"]) / len(mitochondria[query])
        tcov = covered_length(intervals["target"]) / len(assemblies[target])
        remove = qcov >= min_qcov and tcov >= min_tcov
        if remove:
            removed.add(target)
        report.append({"assembly_id": target, "mito_id": query,
                       "query_coverage": f"{qcov:.8f}",
                       "target_coverage": f"{tcov:.8f}",
                       "decision": "remove" if remove else "keep"})
    return removed, report


def integrate(assembly, mito, paf, sample, output, removed_path, report_path,
              min_identity=0.95, min_qcov=0.80, min_tcov=0.95):
    for name, value in [("identity", min_identity), ("query coverage", min_qcov),
                        ("target coverage", min_tcov)]:
        if not 0 < value <= 1:
            raise ValueError(f"Minimum {name} must be in (0, 1]")
    assembly_records = read_fasta(assembly)
    mito_records = read_fasta(mito)
    assemblies = {header.split()[0]: seq for header, seq in assembly_records}
    mitochondria = {header.split()[0]: seq for header, seq in mito_records}
    removed, report = mitochondrial_matches(
        paf, assemblies, mitochondria, min_identity, min_qcov, min_tcov)
    mito_names = ([f"{sample}_mito"] if len(mito_records) == 1 else
                  [f"{sample}_mito_{i}" for i in range(1, len(mito_records) + 1)])
    kept = [(header, seq) for header, seq in assembly_records
            if header.split()[0] not in removed]
    if set(mito_names) & {header.split()[0] for header, _ in kept}:
        raise ValueError("Renamed mitochondrial IDs collide with retained assembly IDs")
    for path in [output, removed_path, report_path]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(output).open("w") as handle:
        for header, seq in kept + list(zip(mito_names, [seq for _, seq in mito_records])):
            handle.write(f">{header}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")
    Path(removed_path).write_text("".join(name + "\n" for name in sorted(removed)))
    with Path(report_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
            "assembly_id", "mito_id", "query_coverage", "target_coverage", "decision"])
        writer.writeheader()
        writer.writerows(report)
    return removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["assembly", "mito", "paf", "sample", "output", "removed", "report"]:
        parser.add_argument("--" + option, required=True)
    parser.add_argument("--min-identity", type=float, default=0.95)
    parser.add_argument("--min-query-coverage", type=float, default=0.80)
    parser.add_argument("--min-target-coverage", type=float, default=0.95)
    args = parser.parse_args()
    removed = integrate(args.assembly, args.mito, args.paf, args.sample, args.output,
                        args.removed, args.report, args.min_identity,
                        args.min_query_coverage, args.min_target_coverage)
    print(f"Removed {len(removed)} whole mitochondrial assembly records; "
          f"wrote final FASTA to {args.output}")


if __name__ == "__main__":
    main()
