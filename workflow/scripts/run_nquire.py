#!/usr/bin/env python3
"""Run nQuire safely and report a fixed-model comparison with explicit status."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import traceback


HEADER = ["file", "free", "dip", "tri", "tet", "d_dip", "d_tri", "d_tet"]
INTERPRETATION = "exploratory_fixed_model_comparison_not_validated_biological_ploidy"


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_bin(path):
    """Validate the network-endian format in pinned nQuire dump_utils.c.

    Headers contain uint32 positions, uint8 extended flag, and uint32 minimum
    fraction times 1000, followed by a uint32 BAM header hash when extended.
    Each position optionally has two uint32 coordinates, then uint16 coverage,
    uint8 allele count, and that many uint16 allele-support counts.
    """
    path = Path(path)
    with path.open("rb") as handle:
        header = handle.read(9)
        if len(header) != 9:
            raise ValueError(f"Truncated nQuire bin header: {path}")
        positions, extended, minimum = struct.unpack(">IBI", header)
        if extended not in (0, 1) or not 0 <= minimum < 500:
            raise ValueError(f"Invalid nQuire bin header: {path}")
        bam_header_hash = None
        if extended:
            value = handle.read(4)
            if len(value) != 4:
                raise ValueError(f"Truncated extended nQuire bin header: {path}")
            bam_header_hash = struct.unpack(">I", value)[0]
        for index in range(positions):
            if extended and len(handle.read(8)) != 8:
                raise ValueError(f"Truncated nQuire position coordinates: {path}")
            record = handle.read(3)
            if len(record) != 3:
                raise ValueError(f"Truncated nQuire position {index + 1}: {path}")
            coverage, alleles = struct.unpack(">HB", record)
            if not coverage or not 1 <= alleles <= 6:
                raise ValueError(f"Invalid nQuire position {index + 1}: {path}")
            counts = handle.read(2 * alleles)
            if len(counts) != 2 * alleles:
                raise ValueError(f"Truncated nQuire allele counts: {path}")
            counts = struct.unpack(">" + "H" * alleles, counts)
            if any(count == 0 or count > coverage for count in counts) or sum(counts) > coverage:
                raise ValueError(f"Invalid nQuire allele counts: {path}")
        if handle.read(1):
            raise ValueError(f"Unexpected trailing nQuire bin data: {path}")
    return {"positions": positions, "extended": bool(extended),
            "minimum_allele_fraction": minimum / 1000,
            "bam_header_hash": bam_header_hash, "bytes": path.stat().st_size,
            "sha256": sha256(path)}


def parse_results(path, expected_label):
    """Require one complete finite row; preserve ties without choosing diploidy."""
    rows = [line.split("\t") for line in Path(path).read_text().splitlines() if line.strip()]
    if len(rows) != 2 or rows[0] != HEADER or len(rows[1]) != len(HEADER):
        raise ValueError("nQuire results require the exact header and one eight-column model row")
    if rows[1][0] != expected_label:
        raise ValueError("nQuire results name a different input bin")
    try:
        values = dict(zip(HEADER[1:], map(float, rows[1][1:])))
    except ValueError as error:
        raise ValueError("nQuire likelihoods must be numeric") from error
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("nQuire likelihoods and deltas must be finite")
    for model in ("dip", "tri", "tet"):
        if not math.isclose(values["free"] - values[model], values["d_" + model],
                            rel_tol=1e-10, abs_tol=2e-6):
            raise ValueError("nQuire delta does not equal free minus fixed-model likelihood")
    deltas = {"2": values["d_dip"], "3": values["d_tri"], "4": values["d_tet"]}
    minimum = min(deltas.values())
    tied = [model for model, delta in deltas.items() if delta == minimum]
    return {"status": "ranked_fixed_models" if len(tied) == 1 else "ambiguous_tie",
            "best_supported_fixed_model": tied[0] if len(tied) == 1 else None,
            "minimum_delta_models": tied, "fixed_model_deltas": deltas,
            "native_likelihoods": values}


def validate_analysis(results, sample, raw_bin, denoised_bin):
    """Validate saved native outputs without rerunning the scientific analysis."""
    raw, denoised = read_bin(raw_bin), read_bin(denoised_bin)
    if not raw["positions"] or not denoised["positions"]:
        raise ValueError("Model comparison requires nonzero raw and denoised position counts")
    validate_denoised(raw, denoised)
    return {**parse_results(results, sample + "-denoised.bin"),
            "raw_bin": raw, "denoised_bin": denoised}


def validate_denoised(raw, denoised):
    if denoised["positions"] > raw["positions"]:
        raise ValueError("nQuire denoising increased the position count")
    for key in ("extended", "minimum_allele_fraction", "bam_header_hash"):
        if raw[key] != denoised[key]:
            raise ValueError(f"nQuire denoising changed bin header field {key}")


def source_revision(executable):
    if not (executable.parent / ".git").exists() or not shutil.which("git"):
        return None
    result = subprocess.run(["git", "-C", str(executable.parent), "rev-parse", "HEAD"],
                            capture_output=True, text=True, timeout=10)
    revision = result.stdout.strip()
    return revision if result.returncode == 0 and re.fullmatch(r"[0-9a-f]{40}", revision) else None


def run(executable, bam, outdir, sample, threads, revision=""):
    if threads < 1 or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
        raise ValueError("A positive thread allocation and safe sample identifier are required")
    executable, bam = Path(executable).resolve(strict=True), Path(bam).resolve(strict=True)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("nQuire executable is not an executable file")
    if not bam.is_file() or not bam.stat().st_size:
        raise ValueError("Input BAM must be a nonempty file")
    outdir = Path(outdir).absolute()
    outdir.parent.mkdir(parents=True, exist_ok=True)
    if outdir.is_symlink() or (outdir.exists() and not outdir.is_dir()):
        raise ValueError("Refusing to replace a symlink or non-directory nQuire output")
    work = Path(tempfile.mkdtemp(prefix=f".{sample}.nquire.", dir=outdir.parent))
    previous = None
    summary = {"schema_version": 1, "state": "running", "sample": sample,
               "interpretation": INTERPRETATION, "started_utc": utcnow(),
               "threads": threads, "commands": [], "scratch": str(work)}
    try:
        # Archive earlier results before a rerun so a failed attempt cannot leave
        # an old completion marker beside stale native bins or model rows.
        if outdir.exists():
            if any(outdir.iterdir()):
                previous = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.previous.", dir=outdir.parent))
                os.replace(outdir, previous)
            else:
                outdir.rmdir()
        summary.update(previous_output_directory=str(previous) if previous else None,
                       bam={"path": str(bam), "bytes": bam.stat().st_size, "sha256": sha256(bam)},
                       tool={"executable": str(executable), "sha256": sha256(executable),
                             "configured_revision": revision or None,
                             "source_checkout_revision": source_revision(executable)},
                       implementation={"script": str(Path(__file__).resolve()), "sha256": sha256(__file__)})
        log = work / f"{sample}.nquire.log"
        results = work / f"{sample}_nquire.results"

        def execute(step, arguments, stdout_path=None):
            command = [str(executable), *arguments]
            record = {"step": step, "command": command, "started_utc": utcnow()}
            summary["commands"].append(record)
            with log.open("ab") as diagnostics:
                diagnostics.write(("\n[" + step + "] " + shlex.join(command) + "\n").encode())
                diagnostics.flush()
                if stdout_path is None:
                    result = subprocess.run(command, cwd=work, stdout=diagnostics, stderr=diagnostics)
                else:
                    with stdout_path.open("wb") as output:
                        result = subprocess.run(command, cwd=work, stdout=output, stderr=diagnostics)
            record.update(exit_code=result.returncode, finished_utc=utcnow())
            result.check_returncode()

        raw_path, denoised_path = work / f"{sample}.bin", work / f"{sample}-denoised.bin"
        execute("create", ["create", "-b", str(bam), "-o", sample])
        summary["raw_bin"] = raw = read_bin(raw_path)
        summary["denoised_bin"] = None
        if not raw["positions"]:
            status = "untestable_zero_raw_positions"
        else:
            execute("denoise", ["denoise", raw_path.name, "-o", sample + "-denoised"])
            summary["denoised_bin"] = denoised = read_bin(denoised_path)
            validate_denoised(raw, denoised)
            if not denoised["positions"]:
                status = "untestable_zero_denoised_positions"
            else:
                execute("lrdmodel", ["lrdmodel", "-t", str(threads), denoised_path.name], results)
                summary.update(parse_results(results, denoised_path.name))
                status = summary["status"]
        summary["status"] = status
        summary["lrdmodel_executed"] = any(item["step"] == "lrdmodel" for item in summary["commands"])
        if status.startswith("untestable_"):
            summary.update(best_supported_fixed_model=None, minimum_delta_models=[],
                           fixed_model_deltas=None, native_likelihoods=None)
            # An explicit skipped report is distinct from a successful model row.
            results.write_text(f"# NOT_RUN: {status}\n" + "\t".join(HEADER) + "\n")
        model = summary["best_supported_fixed_model"] or "NA"
        (work / f"{sample}_ploidy.tsv").write_text(f"{sample}\t{model}\n")
        summary.update(state="completed", finished_utc=utcnow(),
                       outputs={path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                                for path in sorted(work.iterdir()) if path.is_file()})
        write_json(work / f"{sample}_nquire.summary.json", summary)
        (work / "nquire.done").touch()
        os.replace(work, outdir)
        print(json.dumps({"sample": sample, "status": status, "best_supported_fixed_model":
                          summary["best_supported_fixed_model"], "output_directory": str(outdir)}), flush=True)
        return summary
    except BaseException as error:
        summary.update(state="failed_incomplete", error=str(error), traceback=traceback.format_exc(),
                       finished_utc=utcnow())
        write_json(work / "failure.json", summary)
        print(f"nQuire failed; diagnostic workspace retained at {work}", file=sys.stderr, flush=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "bam", "outdir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--revision", default="")
    args = parser.parse_args()
    run(args.executable, args.bam, args.outdir, args.sample, args.threads, args.revision)


if __name__ == "__main__":
    main()
