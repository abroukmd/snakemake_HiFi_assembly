#!/usr/bin/env python3
"""Run pinned MGA from HiFi reads in a fresh workspace, retaining native stages."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def open_reads(path):
    with Path(path).open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="ascii", newline="") if compressed else Path(path).open(
        "r", encoding="ascii", newline="")


def clean_reads(source, destination):
    """Stream FASTA/FASTQ, trim only header metadata, and preserve bases/qualities."""
    reads, bases, trimmed = 0, 0, 0
    seen = set()
    content_hash = hashlib.sha256()

    def write_record(out, header, sequence, quality=None):
        nonlocal reads, bases, trimmed
        identifier = header[1:].split()[0] if header[1:].split() else ""
        if not identifier:
            raise ValueError("Read header has an empty identifier")
        if identifier in seen:
            raise ValueError(f"Duplicate read identifier after trimming: {identifier}")
        seen.add(identifier)
        if not sequence or any(c.isspace() for c in sequence):
            raise ValueError(f"Empty sequence or whitespace in bases for {identifier}")
        if quality is not None and (len(quality) != len(sequence) or
                                    any(not 33 <= ord(c) <= 126 for c in quality)):
            raise ValueError(f"Invalid FASTQ quality for {identifier}")
        text = f"{header[0]}{identifier}\n{sequence}\n"
        if quality is not None:
            text += f"+\n{quality}\n"
        out.write(text)
        content_hash.update(text.encode("ascii"))
        reads += 1
        bases += len(sequence)
        trimmed += header[1:] != identifier

    with open_reads(source) as handle, Path(destination).open("w", encoding="ascii", newline="\n") as out:
        first = handle.readline().rstrip("\r\n")
        if not first.startswith((">", "@")):
            raise ValueError("HiFi reads must be nonempty FASTA or FASTQ")
        format_name = "fasta" if first.startswith(">") else "fastq"
        if format_name == "fasta":
            header, sequence = first, []
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    write_record(out, header, "".join(sequence))
                    header, sequence = line, []
                elif line:
                    sequence.append(line)
            write_record(out, header, "".join(sequence))
        else:
            header = first
            while header:
                if not header.startswith("@"):
                    raise ValueError("Expected a FASTQ @ header")
                sequence = []
                while True:
                    line = handle.readline()
                    if not line:
                        raise ValueError("Truncated FASTQ sequence")
                    line = line.rstrip("\r\n")
                    if line.startswith("+"):
                        break
                    sequence.append(line)
                sequence = "".join(sequence)
                quality = []
                length = 0
                while length < len(sequence):
                    line = handle.readline()
                    if not line:
                        raise ValueError("Truncated FASTQ quality")
                    line = line.rstrip("\r\n")
                    quality.append(line)
                    length += len(line)
                write_record(out, header, sequence, "".join(quality))
                line = handle.readline()
                if line and not line.rstrip("\r\n"):
                    raise ValueError("Unexpected blank line between FASTQ records")
                header = line.rstrip("\r\n")
    return {"format": format_name, "reads": reads, "bases": bases,
            "headers_trimmed": trimmed, "prepared_reads_sha256": content_hash.hexdigest()}


def validate_fasta(path):
    count, length, current = 0, 0, 0
    identifiers = set()
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count and not current:
                    raise ValueError("MGA assembly contains an empty FASTA record")
                identifier = line[1:].split()[0] if line[1:].split() else ""
                if not identifier or identifier in identifiers:
                    raise ValueError("MGA assembly has an empty or duplicate FASTA ID")
                identifiers.add(identifier)
                count += 1
                current = 0
            else:
                if not count or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                    raise ValueError("MGA assembly has invalid nucleotide FASTA content")
                current += len(line)
                length += len(line)
    if not count or not current:
        raise ValueError("MGA polished assembly is empty")
    return {"records": count, "bases": length}


def prepare_output_links(native, scratch, prepared):
    """Make native graph links survive publication outside the scratch tree."""
    native, scratch, prepared = map(lambda path: Path(path).absolute(), (native, scratch, prepared))
    links = []
    for directory, dirs, files in os.walk(native, followlinks=False):
        links.extend(Path(directory) / name for name in dirs + files
                     if (Path(directory) / name).is_symlink())
    report = {"internal_links_relocated": 0, "prepared_reads_retained": False,
              "external_links": []}
    for link in sorted(links):
        target = Path(os.readlink(link))
        target = Path(os.path.abspath(target if target.is_absolute() else link.parent / target))
        if not link.exists():
            raise ValueError(f"MGA created a dangling or cyclic output link: {link} -> {target}")
        resolved = target.resolve()
        if target.is_relative_to(native):
            destination = target
            report["internal_links_relocated"] += 1
        elif resolved.is_relative_to(native):
            destination = resolved
            report["internal_links_relocated"] += 1
        elif target.is_relative_to(scratch) or resolved.is_relative_to(scratch):
            if resolved != prepared.resolve():
                raise ValueError(f"MGA output links to an unpreserved scratch artifact: {link} -> {target}")
            destination = native / "retained_inputs" / prepared.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                shutil.copyfile(prepared, destination)
            report["prepared_reads_retained"] = True
        else:
            report["external_links"].append({"link": str(link.relative_to(native)),
                                              "target": str(target)})
            continue
        link.unlink()
        link.symlink_to(os.path.relpath(destination, link.parent), target_is_directory=destination.is_dir())
    # Validate after all chains have been rewritten, before publishing success.
    for link in links:
        if not link.exists():
            raise ValueError(f"MGA graph link is invalid after relocation: {link}")
    return report


def preserve_failed_workspace(scratch, outparent, receipt):
    """Retain all failure evidence, including targets outside native/ itself.

    Keeping the complete workspace avoids losing prepared reads or temporary
    files referenced by a partial graph. Rewrite absolute internal links before
    moving the tree; pre-existing dangling links remain documented as such.
    """
    scratch = Path(scratch).absolute()
    links = []
    for directory, dirs, files in os.walk(scratch, followlinks=False):
        links.extend(Path(directory) / name for name in dirs + files
                     if (Path(directory) / name).is_symlink())
    report = []
    for link in sorted(links):
        original = os.readlink(link)
        target = Path(original)
        target = Path(os.path.abspath(target if target.is_absolute() else link.parent / target))
        existed = link.exists()
        destination = target if target.is_relative_to(scratch) else None
        if destination is None and existed:
            resolved = target.resolve()
            if resolved.is_relative_to(scratch):
                destination = resolved
        entry = {"link": str(link.relative_to(scratch)), "original_target": original,
                 "target_existed_before_preservation": existed}
        if destination is not None:
            relative = os.path.relpath(destination, link.parent)
            link.unlink()
            link.symlink_to(relative, target_is_directory=destination.is_dir())
            entry["retained_target"] = relative
        report.append(entry)
    failed = Path(tempfile.mkdtemp(prefix=".mga-failed-", dir=outparent))
    os.replace(scratch, failed)
    receipt["failed_workspace_links"] = report
    (failed / "mga_run.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return failed


def native_cpu_allocation(threads):
    """Choose only CPUs already permitted by the scheduler or parent process."""
    if not callable(getattr(os, "sched_getaffinity", None)) or not callable(
            getattr(os, "sched_setaffinity", None)):
        raise RuntimeError("MGA requires Linux CPU affinity support to bound nested tool CPU use")
    try:
        inherited = sorted(os.sched_getaffinity(0))
    except OSError as error:
        raise RuntimeError("Cannot inspect inherited CPU affinity for MGA") from error
    if threads < 1 or not inherited:
        raise ValueError("MGA needs a positive allocation and at least one inherited CPU")
    return {"method": "Linux sched_setaffinity in native child before exec; descendants inherit",
            "inherited_cpus": inherited, "selected_cpus": inherited[:threads]}


def exec_native_with_affinity(cpus, command):
    """Set affinity in a fresh child interpreter, then replace it with MGA.

    This avoids a Python preexec_fn in a forked process and leaves the calling
    workflow process's affinity unchanged.
    """
    try:
        selected = set(cpus)
        if not selected or not selected <= os.sched_getaffinity(0):
            raise RuntimeError("Requested native CPUs are outside the inherited allocation")
        os.sched_setaffinity(0, selected)
        if os.sched_getaffinity(0) != selected:
            raise RuntimeError("Native CPU affinity differs from the requested bounded allocation")
    except (AttributeError, OSError) as error:
        raise RuntimeError("Cannot enforce the native MGA CPU affinity limit") from error
    if not command:
        raise ValueError("Missing MGA native command")
    os.execv(command[0], command)


def run(executable, reads, outdir, threads, revision):
    executable, reads = Path(executable).resolve(strict=True), Path(reads).resolve(strict=True)
    outdir = Path(outdir).absolute()
    if threads < 1:
        raise ValueError("MGA thread allocation must be positive")
    cpu_allocation = native_cpu_allocation(threads)
    # Upstream constructs internal shell commands without quoting these paths.
    for path in [executable, outdir.parent]:
        if not re.fullmatch(r"[A-Za-z0-9_./-]+", str(path)):
            raise ValueError("MGA installation/output paths must contain only letters, digits, _, ., /, -")
    if outdir.is_symlink():
        raise ValueError(f"Refusing symlink output directory: {outdir}")
    if outdir.exists() and any(not path.is_dir() or path.is_symlink() for path in outdir.rglob("*")):
        raise FileExistsError(f"MGA output directory already contains files: {outdir}")
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".mga-run-", dir=outdir.parent) as temporary:
        scratch = Path(temporary)
        prepared = scratch / "reads.cleaned"
        read_stats = clean_reads(reads, prepared)
        prepared = prepared.rename(scratch / ("reads.cleaned." + read_stats["format"]))
        native = scratch / "native"
        command = [str(executable), "--reads", str(prepared), "--output", str(native),
                   "--threads", str(threads)]
        native_launcher = [sys.executable, str(Path(__file__).resolve()), "--native-cpu-affinity",
                           ",".join(map(str, cpu_allocation["selected_cpus"])), "--", *command]
        receipt = {"revision": revision, "executable": str(executable),
                   "reads": str(reads), "command": command, "threads": threads,
                   "cpu_allocation": cpu_allocation, "native_launcher": native_launcher,
                   "input_reads": read_stats,
                   "native_assembly": "5_polishing/assembly.fasta"}
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        print(json.dumps(receipt), flush=True)
        try:
            subprocess.run(native_launcher, env=env, cwd=scratch, check=True)
            receipt["graph_links"] = prepare_output_links(native, scratch, prepared)
            polished = native / "5_polishing/assembly.fasta"
            receipt["assembly"] = validate_fasta(polished)
            receipt["assembly_sha256"] = hashlib.sha256(polished.read_bytes()).hexdigest()
            shutil.copyfile(polished, native / "assembly.fasta")
            (native / "mga_run.json").write_text(json.dumps(receipt, indent=2) + "\n")
            (native / "mga.done").write_text(revision + "\n")
        except Exception:
            # Preserve the complete scratch layout: native partial graphs can
            # reference both their own records and prepared/temporary inputs.
            if native.exists():
                failed = preserve_failed_workspace(scratch, outdir.parent, receipt)
                print(f"Failed MGA workspace retained at {failed}", flush=True)
            raise
        if outdir.is_symlink():
            raise ValueError(f"Refusing symlink output directory: {outdir}")
        if outdir.exists():
            # Snakemake owns this directory; direct callers must also use a
            # fresh output location to avoid erasing earlier graph evidence.
            # Its output preparation can create empty nested parent folders.
            if any(not path.is_dir() or path.is_symlink() for path in outdir.rglob("*")):
                raise FileExistsError(f"MGA output directory already contains files: {outdir}")
            shutil.rmtree(outdir)
        os.replace(native, outdir)
    return receipt


def main():
    if sys.argv[1:2] == ["--native-cpu-affinity"]:
        parser = argparse.ArgumentParser(description="Internal bounded MGA child launcher")
        parser.add_argument("--native-cpu-affinity", required=True)
        parser.add_argument("command", nargs=argparse.REMAINDER)
        args = parser.parse_args()
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        exec_native_with_affinity([int(cpu) for cpu in args.native_cpu_affinity.split(",")], command)
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    run(args.executable, args.reads, args.outdir, args.threads, args.revision)


if __name__ == "__main__":
    main()
