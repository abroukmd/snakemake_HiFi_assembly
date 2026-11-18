import os
import re

# Source is immutable and installed separately for each configured revision.
MGA_REPO = "https://github.com/ZhangZhenmiao/MGA.git"
MGA_COMMIT = str(config.get("mga_commit", "f0f26da124e5e1362875dbf15c0de1bec2fd34f8")).lower()
if not re.fullmatch(r"[0-9a-f]{40}", MGA_COMMIT):
    raise ValueError("mga_commit must be a complete immutable 40-character Git SHA")
MGA_PREFIX = outpath("resources", "mga", MGA_COMMIT)
MGA_SRC = os.path.join(MGA_PREFIX, "src")
MGA_BIN = os.path.join(MGA_SRC, "bin", "MGA")
MGA_STAMP = os.path.join(MGA_PREFIX, ".installed")


rule install_mga:
    output:
        done=MGA_STAMP,
        binary=MGA_BIN,
        lja_cdb=os.path.join(MGA_SRC, "lib/LJA/bin/lja_cdb"),
        multidbg=os.path.join(MGA_SRC, "lib/LJA/bin/multiDBG"),
        jumbodbg=os.path.join(MGA_SRC, "lib/LJA/bin/jumboDBG"),
        align=os.path.join(MGA_SRC, "lib/LJA/bin/align_and_print"),
        compress=os.path.join(MGA_SRC, "lib/LJA/bin/compress"),
        polish=os.path.join(MGA_SRC, "lib/LJA/bin/run_polishing")
    params:
        repo=MGA_REPO,
        commit=MGA_COMMIT,
        prefix=MGA_PREFIX,
        src=MGA_SRC
    log:
        os.path.join(MGA_PREFIX, "install.log")
    conda:
        "../envs/mga_env.yaml"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.prefix:q}
        exec > {log:q} 2>&1
        if [[ ! -d {params.src:q}/.git ]]; then
            git init {params.src:q}
            git -C {params.src:q} remote add origin {params.repo:q}
        fi
        test "$(git -C {params.src:q} config --get remote.origin.url)" = {params.repo:q}
        git -C {params.src:q} fetch --depth 1 origin {params.commit:q}
        git -C {params.src:q} checkout --detach FETCH_HEAD
        test "$(git -C {params.src:q} rev-parse HEAD)" = {params.commit:q}
        # Upstream resolves bundled LJA executables and Python scripts relative
        # to bin/MGA, so preserve its documented in-source build layout.
        cmake -S {params.src:q} -B {params.src:q} -DCMAKE_BUILD_TYPE=Release
        cmake --build {params.src:q} --parallel {threads}
        for binary in {output.binary:q} {output.lja_cdb:q} {output.multidbg:q} {output.jumbodbg:q} \
                      {output.align:q} {output.compress:q} {output.polish:q}; do
            test -x "$binary"
        done
        {output.binary:q} --help
        python -c 'import Bio, numpy, pysam'
        minimap2 --version
        samtools --version
        printf '%s\n' {params.commit:q} > {output.done:q}
        """


rule mga:
    input:
        reads=lambda wc: sample_paths[wc.sample],
        installed=rules.install_mga.output.done,
        binary=rules.install_mga.output.binary,
        lja_cdb=rules.install_mga.output.lja_cdb,
        multidbg=rules.install_mga.output.multidbg,
        jumbodbg=rules.install_mga.output.jumbodbg,
        align=rules.install_mga.output.align,
        compress=rules.install_mga.output.compress,
        polish=rules.install_mga.output.polish,
        helper=workflow.source_path("../scripts/run_mga.py")
    output:
        results=directory(outpath("Assemblies/{sample}/MGA")),
        fasta=outpath("Assemblies/{sample}/MGA/assembly.fasta"),
        native=outpath("Assemblies/{sample}/MGA/5_polishing/assembly.fasta"),
        done=outpath("Assemblies/{sample}/MGA/mga.done")
    params:
        commit=MGA_COMMIT
    log:
        outpath("Assemblies/{sample}/MGA.log")
    conda:
        "../envs/mga_env.yaml"
    threads: 32
    resources:
        mem_mb=int(config.get("mga_mem_mb", 64000))
    shell:
        r"""
        set -euo pipefail
        python {input.helper:q} --executable {input.binary:q} --reads {input.reads:q} \
            --outdir {output.results:q} --threads {threads} --revision {params.commit:q} \
            > {log:q} 2>&1
        """
