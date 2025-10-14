rule multiqc:
    input:
        expand(outpath("Assemblies/{sample}/QC/nanoplot/nanoplot.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/QUAST/quast.done"), sample=SAMPLES),
        expand(outpath("Assemblies/{sample}/QC/merqury/merqury.done"), sample=SAMPLES)
    output:
        html = outpath("Assemblies/multiqc/multiqc_report.html"),
        done = outpath("Assemblies/multiqc/multiqc.done")
    log:
        stdout = outpath("Assemblies/multiqc/multiqc.log")
    conda:
        "../envs/multiqc_env.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p Assemblies/multiqc
        multiqc Assemblies/ -o Assemblies/multiqc &> {log.stdout}
        touch {output.done}
        """
