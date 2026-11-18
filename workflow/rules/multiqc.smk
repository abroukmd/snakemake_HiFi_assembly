rule multiqc:
    input:
        expand(outpath("Assemblies/{sample}/QC/nanoplot/nanoplot.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done"), sample=HIFIASM_SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done"), sample=LJA_SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.done"), sample=HIFIASM_SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.done"), sample=LJA_SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/QUAST/quast.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/merqury/merqury.done"), sample=SAMPLES)
    output:
        html = outpath("Assemblies/multiqc/multiqc_report.html"),
        done = outpath("Assemblies/multiqc/multiqc.done")
    log:
        stdout = outpath("Assemblies/multiqc/multiqc.log")
    conda:
        "../envs/multiqc_env.yaml"
    threads: 1
    params:
        qc_dirs = [outpath(f"Assemblies/{sample}/QC") for sample in SAMPLES],
        outdir = outpath("Assemblies/multiqc")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir:q}
        # BUSCO exposes run_* directories as symlink aliases of auto_lineage.
        # Prune both spellings; retain generic/specific summaries at result roots.
        multiqc {params.qc_dirs:q} --dirs --dirs-depth 3 \
            --ignore '*/auto_lineage/*' \
            --ignore '*/BUSCO-*/*/run_*' \
            --ignore '*/BUSCO-*/*/run_*/*' \
            --force -o {params.outdir:q} &> {log.stdout:q}
        test -s {output.html:q}
        touch {output.done}
        """
