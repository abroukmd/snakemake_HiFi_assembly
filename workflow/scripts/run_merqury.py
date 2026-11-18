#!/usr/bin/env python3
"""Run one assembly in a fresh Merqury workspace and validate real outputs."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess
import tempfile


PINNED_PLOT_SHA256 = "3a8701dc6ac6eff1f9ea1e983e6fe9272b20a31e81e4abc0cc4cf893e28d17ae"
PLOT_ZERO_BLOCK = '''    dat_0=read.table(zero, header=FALSE)
    dat_0[,1]=as.factor(dat_0[,1])
    dat_0[,1]=factor(dat_0[,1], levels=rev(unique(dat_0[,1])), ordered=TRUE)'''
PLOT_ZERO_BLOCK_FIXED = '''    dat_0=read.table(zero, header=FALSE)
    # Omit zero-height bars only in memory; retain native histogram files.
    dat_0=dat_0[dat_0[,3] > 0,,drop=FALSE]
    if (nrow(dat_0) == 0) {
      dat_0=NULL
    } else {
      dat_0[,1]=as.factor(dat_0[,1])
      dat_0[,1]=factor(dat_0[,1], levels=rev(unique(dat_0[,1])), ordered=TRUE)
    }'''
PLOT_YMAX_BLOCK = '''  if (zero != "") {
    y_max=max(y_max, sum(dat_0[,3]*1.1))'''
PLOT_YMAX_BLOCK_FIXED = '''  if (!is.null(dat_0)) {
    y_max=max(y_max, sum(dat_0[,3]*1.1))'''


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def resolve_tools():
    """Use the Merqury tree selected by PATH, never an inherited MERQURY tree."""
    executable = shutil.which("merqury.sh")
    meryl = shutil.which("meryl")
    if not executable or not meryl:
        raise FileNotFoundError("merqury.sh and meryl must be available")
    executable, meryl = Path(executable).resolve(strict=True), Path(meryl).resolve(strict=True)
    root = executable.parent
    sources = [executable, root / "util/util.sh", root / "eval/spectra-cn.sh",
               root / "plot/plot_spectra_cn.R"]
    for source in sources:
        if not source.is_file():
            raise FileNotFoundError(f"Incomplete Merqury installation selected by PATH: {source}")
    return executable, meryl, {
        "merqury_root": str(root),
        "merqury_executable": str(executable),
        "meryl_executable": str(meryl),
        "source_sha256": {str(path): sha256(path) for path in [*sources, meryl]},
    }


def prepare_runtime(root, work):
    """Adapt only the pinned plot script in a private retained package overlay.

    Merqury emits assembly-only rows with zero counts, which can add a colour
    level absent from the main histogram and fail ggplot's manual palette.
    Native histograms, databases, metrics, and the installed package stay intact.
    """
    root, work = Path(root), Path(work)
    source = root / "plot/plot_spectra_cn.R"
    original_hash = sha256(source)
    if original_hash != PINNED_PLOT_SHA256:
        raise ValueError(f"Unsupported Merqury plot source SHA256 {original_hash}: {source}; "
                         "the private plotting adaptation requires the pinned Merqury 1.4.1 source")
    content = source.read_text()
    if content.count(PLOT_ZERO_BLOCK) != 1 or content.count(PLOT_YMAX_BLOCK) != 1:
        raise ValueError("Pinned Merqury plot source does not match the guarded adaptation")
    patched = content.replace(PLOT_ZERO_BLOCK, PLOT_ZERO_BLOCK_FIXED).replace(
        PLOT_YMAX_BLOCK, PLOT_YMAX_BLOCK_FIXED)
    runtime = work / "merqury-runtime"
    runtime.mkdir()
    for child in root.iterdir():
        target = runtime / child.name
        if child.name == "plot":
            target.mkdir()
            for plot_source in child.iterdir():
                plot_target = target / plot_source.name
                if plot_source.name == source.name:
                    plot_target.write_text(patched)
                else:
                    plot_target.symlink_to(plot_source.resolve(), target_is_directory=plot_source.is_dir())
        else:
            target.symlink_to(child.resolve(), target_is_directory=child.is_dir())
    patched_source = runtime / "plot" / source.name
    return runtime, {
        "purpose": "Omit zero-height assembly-only bars in plot memory; preserve native histograms and metrics",
        "original_path": str(source), "original_sha256": original_hash,
        "patched_path": str(patched_source.relative_to(work)),
        "patched_sha256": sha256(patched_source),
    }


def normalize_completeness(work):
    """Retain modern native output and expose the workflow's stable filename."""
    native, stable = work / "merqury.completeness.stats", work / "completeness.stats"
    if native.is_file():
        if stable.exists() and stable.read_bytes() != native.read_bytes():
            raise ValueError("Conflicting native Merqury completeness reports")
        shutil.copyfile(native, stable)


def validate_reports(work, label):
    """Reject missing, duplicate, mislabelled, or undefined single-assembly metrics."""
    qv = [line.split() for line in (work / "merqury.qv").read_text().splitlines() if line.strip()]
    completeness = [line.split() for line in (work / "completeness.stats").read_text().splitlines() if line.strip()]
    if len(qv) != 1 or len(qv[0]) != 5 or qv[0][0] != label:
        raise ValueError("Merqury QV must contain exactly the labelled assembly")
    if len(completeness) != 1 or len(completeness[0]) != 5 or completeness[0][0] != label:
        raise ValueError("Merqury completeness must contain exactly the labelled assembly")
    errors, total = int(qv[0][1]), int(qv[0][2])
    quality, error_rate = float(qv[0][3]), float(qv[0][4])
    recovered, solid = int(completeness[0][2]), int(completeness[0][3])
    percent = float(completeness[0][4])
    if total <= 0 or not 0 <= errors <= total or math.isnan(quality) or quality < 0:
        raise ValueError("Undefined or inconsistent Merqury QV")
    if not math.isfinite(error_rate) or not 0 <= error_rate <= 1:
        raise ValueError("Invalid Merqury error rate")
    if solid <= 0 or not 0 <= recovered <= solid or not math.isfinite(percent) or not 0 <= percent <= 100:
        raise ValueError("Undefined or inconsistent Merqury completeness")
    png = work / "merqury.spectra-asm.ln.png"
    with png.open("rb") as handle:
        if handle.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError("Merqury did not produce a valid spectra PNG")
    # Older Merqury versions can exit zero even when their R subprocess fails.
    for log in (work / "logs").glob("*.log"):
        if re.search(r"(?m)^(Error(?: in |:|$)|ERROR:|Unknown option |Execution halted|Traceback \(most recent call last\):)", log.read_text(errors="replace")):
            raise ValueError(f"Merqury reported a nested tool failure: {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--read-db", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--memory-gb", type=int, default=15)
    parser.add_argument("--k", type=int, default=21)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.label):
        parser.error("label must be a safe filename component")
    assembly, reads = args.assembly.resolve(strict=True), args.read_db.resolve(strict=True)
    assembly_sha256 = sha256(assembly)
    outdir = args.outdir.absolute()
    outdir.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    merqury, meryl, tools = resolve_tools()
    with tempfile.TemporaryDirectory(prefix=f".{args.label}.", dir=outdir.parent) as temporary:
        work = Path(temporary)
        runtime, adaptation = prepare_runtime(tools["merqury_root"], work)
        env["MERQURY"] = str(runtime)
        tools.update(runtime_root=str(runtime.relative_to(work)), plot_adaptation=adaptation)
        local_assembly = work / f"{args.label}.fasta"
        local_assembly.symlink_to(assembly)
        # Never let an assembly-derived .meryl name alias the read database.
        local_reads = work / "illumina.reads.meryl"
        local_reads.symlink_to(reads, target_is_directory=True)
        tool_bin = work / "tool-bin"
        tool_bin.mkdir()
        wrapper = tool_bin / "meryl"
        wrapper.write_text("#!/bin/sh\nexec " + shlex.quote(str(meryl)) +
                           f" threads={args.threads} memory={args.memory_gb} \"$@\"\n")
        wrapper.chmod(0o755)
        env["PATH"] = str(tool_bin) + os.pathsep + env["PATH"]
        print(json.dumps({"assembly": str(assembly), "read_db": str(reads),
                          "label": args.label, "k": args.k, "threads": args.threads}), flush=True)
        try:
            subprocess.run([meryl, f"threads={args.threads}", f"memory={args.memory_gb}",
                            f"k={args.k}", "count", "output", f"{args.label}.meryl",
                            str(local_assembly)], cwd=work, env=env, check=True)
            subprocess.run([str(merqury), str(local_reads), str(local_assembly), "merqury"],
                           cwd=work, env=env, check=True)
            normalize_completeness(work)
            validate_reports(work, args.label)
            if sha256(assembly) != assembly_sha256:
                raise ValueError("Assembly changed while Merqury was running")
        except Exception:
            # The tool wrapper can swallow inner errors. Preserve their diagnostics
            # in the rule log before TemporaryDirectory removes the scratch files.
            for nested_log in (work / "logs").glob("*.log"):
                print(f"\n--- {nested_log.name} ---\n" + nested_log.read_text(errors="replace"), flush=True)
            raise
        (work / "inputs.json").write_text(json.dumps({
            "schema_version": 2,
            "assembly": str(assembly), "assembly_sha256": assembly_sha256,
            "read_db": str(reads), "label": args.label,
            "k": args.k, "threads": args.threads,
            "tools": tools,
            "runner_sha256": sha256(__file__),
            "mode": "one assembly per Merqury call"}, indent=2) + "\n")
        # Only this rule owns outdir. Scratch is on the same filesystem for rename.
        if outdir.is_symlink():
            raise ValueError(f"Refusing to replace symlink output directory: {outdir}")
        if outdir.exists():
            shutil.rmtree(outdir)
        os.replace(work, outdir)


if __name__ == "__main__":
    main()
