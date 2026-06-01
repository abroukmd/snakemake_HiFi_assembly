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
        done=NQ_DONE
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
        # ensure mapping exists
        mapping_done=outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done"),
        bam=outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.bam"),
    output:
        results=outpath("Assemblies/{sample}/QC/nQuire/{sample}_nquire.results"),
        ploidy=outpath("Assemblies/{sample}/QC/nQuire/{sample}_ploidy.tsv"),
        done=outpath("Assemblies/{sample}/QC/nQuire/nquire.done"),
    log:
        out=outpath("Assemblies/{sample}/QC/nQuire/{sample}.nquire.log"),
    threads: 32
    conda:
        "../envs/genomescope2_env.yaml"
    shell:
        r"""
        set -euo pipefail
        set -x

        OUTDIR="$(dirname "{output.done}")"
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"

        "{NQ_BIN}" create -b "{input.bam}" -o "{wildcards.sample}" >> "{log.out}" 2>&1
        "{NQ_BIN}" denoise "{wildcards.sample}.bin" -o "{wildcards.sample}-denoised" >> "{log.out}" 2>&1
        "{NQ_BIN}" lrdmodel -t {threads} "{wildcards.sample}-denoised.bin" > "{output.results}" 2>> "{log.out}"

        awk -F'\t' 'NR==1 {{next}} {{min=$6; pl="2"; if($7<min){{min=$7; pl="3"}} if($8<min){{pl="4"}} print $1, pl}}' \
          "{output.results}" | sed 's/-denoised.bin//g' > "{output.ploidy}"

        touch "{output.done}"
        """
