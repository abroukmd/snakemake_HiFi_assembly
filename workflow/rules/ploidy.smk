import os

PLOIDY_RANGE = [1, 2, 3, 4, 5, 6]
KMER_K = 21

PROJECT_DIR = config["project_dir"]
GS2_PREFIX  = os.path.join(PROJECT_DIR, "resources", "genomescope2")
GS2_SRC     = os.path.join(GS2_PREFIX, "src")
GS2_RLIB    = os.path.join(GS2_PREFIX, "Rlib")
GS2_DONE    = os.path.join(GS2_PREFIX, ".installed")

rule genomescope2_ploidy:
    input:
        illumina_r1=lambda wc: sample_fq1[wc.sample],
        illumina_r2=lambda wc: sample_fq2[wc.sample],
        hifi=lambda wc: sample_paths[wc.sample],
        gs2_install=GS2_DONE,   # <<< forces install first
    output:
        histo=outpath("Assemblies/{sample}/QC/GenomeScope2/{sample}.histo"),
        summary=outpath("Assemblies/{sample}/QC/GenomeScope2/{sample}.genomescope2.summary.tsv"),
        done=outpath("Assemblies/{sample}/QC/GenomeScope2/ploidy.done"),
    params:
        k=KMER_K,
        ploidies=" ".join(map(str, PLOIDY_RANGE)),
        genomescope_R=os.path.join(GS2_SRC, "genomescope.R"),
        rlib=GS2_RLIB,
        parse_script=os.path.join(workflow.basedir, "scripts/parse_genomescope2.py"),
        kmc_k=21,
        kmc_ci=1,
        kmc_cs=10000,
        kmc_cx=10000,
        kmc_mem_gb=64,
        kmc_threads=10,
    conda:
        "../envs/genomescope2_env.yaml"
    log:
        kmc=outpath("Assemblies/{sample}/QC/GenomeScope2/{sample}.kmc.log"),
        gs_all=outpath("Assemblies/{sample}/QC/GenomeScope2/{sample}.genomescope2.all.log"),
        parse=outpath("Assemblies/{sample}/QC/GenomeScope2/{sample}.parse.log"),
    shell:
        r"""
        set -euo pipefail

        OUTDIR="$(dirname "{output.histo}")"
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"

        : > "{log.kmc}"
        : > "{log.gs_all}"
        : > "{log.parse}"

        echo "[INFO] Sample: {wildcards.sample}" >> "{log.kmc}"
        echo "[INFO] OUTDIR: $OUTDIR" >> "{log.kmc}"
        echo "[INFO] genomescope.R: {params.genomescope_R}" >> "{log.kmc}"
        echo "[INFO] project Rlib: {params.rlib}" >> "{log.kmc}"

        # Make sure R in this job sees the installed genomescope package
        export R_LIBS_USER="{params.rlib}"

        # Build FILES like your CLI
        : > FILES
        for f in "{input.illumina_r1}" "{input.illumina_r2}" "{input.hifi}"; do
            if [[ -n "$f" && -f "$f" ]]; then
                echo "$f" >> FILES
            else
                echo "[WARN] skipping missing: $f" >> "{log.kmc}"
            fi
        done
        if [[ ! -s FILES ]]; then
            echo "[ERROR] FILES empty" >> "{log.kmc}"
            exit 1
        fi

        # KMC DB + histogram
        rm -rf tmp reads.kmc_pre reads.kmc_suf 2>/dev/null || true
        mkdir -p tmp

        echo "[INFO] Running kmc..." >> "{log.kmc}"
        kmc -k{params.kmc_k} -t{params.kmc_threads} -m{params.kmc_mem_gb} -ci{params.kmc_ci} -cs{params.kmc_cs} \
            @FILES reads tmp/ >> "{log.kmc}" 2>&1

        HISTO_LOCAL="{wildcards.sample}.kmc.histo"
        echo "[INFO] Running kmc_tools histogram -> $HISTO_LOCAL" >> "{log.kmc}"
        kmc_tools transform reads histogram "$HISTO_LOCAL" -cx{params.kmc_cx} >> "{log.kmc}" 2>&1

        if [[ ! -s "$HISTO_LOCAL" ]]; then
            echo "[ERROR] Histogram not created: $HISTO_LOCAL" >> "{log.kmc}"
            ls -lah >> "{log.kmc}" 2>&1
            exit 1
        fi

        # GenomeScope2 loop
        ok=0
        for p in {params.ploidies}; do
            odir="{wildcards.sample}_p${{p}}"
            mkdir -p "$odir"
            echo "[INFO] GenomeScope: p=${{p}} k={params.k} histo=$HISTO_LOCAL outdir=$odir" >> "{log.gs_all}"
            "{params.genomescope_R}" -k {params.k} -i "$HISTO_LOCAL" -o "$odir" -p "${{p}}" >> "{log.gs_all}" 2>&1 || true
            if ls "$odir"/* 1>/dev/null 2>&1; then
                ok=1
            fi
        done

        if [[ "$ok" -ne 1 ]]; then
            echo "[ERROR] GenomeScope2 produced no outputs" >> "{log.gs_all}"
            exit 1
        fi

        # Parse summary
        python "{params.parse_script}" "{output.summary}" "{log.gs_all}" >> "{log.parse}" 2>&1
        if [[ ! -s "{output.summary}" ]]; then
            echo "[ERROR] Summary TSV empty" >> "{log.parse}"
            exit 1
        fi

        # Publish histogram at end
        mv -f "$HISTO_LOCAL" "{output.histo}"
        touch "{output.done}"
        """
