rule multiqc:
    input:
        expand("Assemblies/{sample}/QC/nanoplot/nanoplot.done", sample=SAMPLES),
        expand("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done", sample=SAMPLES),
        expand("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done", sample=SAMPLES),
        expand("Assemblies/{sample}/QC/QUAST/quast.done", sample=SAMPLES),
        expand("Assemblies/{sample}/QC/merqury/merqury.done", sample=SAMPLES)
    output:
        html = "Assemblies/multiqc/multiqc_report.html",
        done = "Assemblies/multiqc/multiqc.done"
    log:
        stdout = "Assemblies/multiqc/multiqc.log"
    conda:
        "../envs/multiqc_env.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p Assemblies/multiqc
        multiqc Assemblies/ -o Assemblies/multiqc &> {log.stdout}
        touch {output.done}
        """
