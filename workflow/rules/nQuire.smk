import os

PROJECT_DIR = config["project_dir"]

# nQuire install location (project-local)
NQ_PREFIX = os.path.join(PROJECT_DIR, "resources", "nQuire")
NQ_SRC    = os.path.join(NQ_PREFIX, "src")
NQ_DONE   = os.path.join(NQ_PREFIX, ".installed")
NQ_BIN    = os.path.join(NQ_SRC, "nQuire")

# optional pinning
NQ_REPO   = config.get("nquire", {}).get("repo", "https://github.com/clwgg/nQuire")
NQ_COMMIT = config.get("nquire", {}).get("commit", "")


rule install_nquire:
    output:
        done=NQ_DONE,
        binary=NQ_BIN
    params:
        repo=NQ_REPO,
        commit=NQ_COMMIT,
        prefix=NQ_PREFIX,
        src=NQ_SRC,
        bin=NQ_BIN,
    conda:
        "../envs/genomescope2_env.yaml"
    log:
        os.path.join(NQ_PREFIX, "install.log")
    shell:
        r"""
        set -euo pipefail
        set -x

        mkdir -p "{params.prefix}"
        touch "{log}"
        exec > >(tee -a "{log}") 2>&1

        echo "[INFO] repo={params.repo}"
        echo "[INFO] commit={params.commit}"
        echo "[INFO] prefix={params.prefix}"
        echo "[INFO] src={params.src}"
        echo "[INFO] bin={params.bin}"

        echo "[INFO] which git:  $(command -v git  || echo MISSING)"
        echo "[INFO] which make: $(command -v make || echo MISSING)"
        echo "[INFO] which gcc:  $(command -v gcc  || echo MISSING)"

        # clone once (with submodules)
        if [[ ! -d "{params.src}/.git" ]]; then
            rm -rf "{params.src}" || true
            git clone --recursive "{params.repo}" "{params.src}"
        fi

        # optional pin
        if [[ -n "{params.commit}" ]]; then
            (cd "{params.src}" && git fetch --all --tags)
            (cd "{params.src}" && git checkout -f "{params.commit}")
        fi

        # ensure submodules present
        (cd "{params.src}" && git submodule update --init --recursive)

        # build
        (cd "{params.src}" && make submodules)
        (cd "{params.src}" && make)

        # verify binary exists
        test -x "{params.bin}"
        "{params.bin}" --help | head -n 5 || true

        touch "{output.done}"
        """


rule nquire_ploidy:
    input:
        # ensure nQuire is installed
        nq_done=NQ_DONE,
        binary=NQ_BIN,
        # ensure mapping exists
        mapping_done=outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done"),
        bam=get_nquire_bam,
        helper=workflow.source_path("../scripts/run_nquire.py"),
    output:
        results=outpath("Assemblies/{sample}/QC/nQuire/{sample}_nquire.results"),
        ploidy=outpath("Assemblies/{sample}/QC/nQuire/{sample}_ploidy.tsv"),
        summary=outpath("Assemblies/{sample}/QC/nQuire/{sample}_nquire.summary.json"),
        done=outpath("Assemblies/{sample}/QC/nQuire/nquire.done"),
    log:
        out=outpath("Assemblies/{sample}/logs/nquire.runner.log"),
    params:
        outdir=lambda wc: outpath(f"Assemblies/{wc.sample}/QC/nQuire"),
        revision=NQ_COMMIT,
    threads: 32
    conda:
        "../envs/genomescope2_env.yaml"
    shell:
        r"""
        set -euo pipefail
        python {input.helper:q} --executable {input.binary:q} --bam {input.bam:q} \
            --outdir {params.outdir:q} --sample {wildcards.sample:q} --threads {threads} \
            --revision {params.revision:q} > {log.out:q} 2>&1
        """
