
rule barrnap:
    input:
        asm_hifiasm = "Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta",
        asm_lja     = "Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta",
        hifiasm_done = "Assemblies/{sample}/HIFIASM/hifiasm.done",
        lja_done     = "Assemblies/{sample}/LJA/lja.done"
    output:
        gff3_hifiasm = "Assemblies/{sample}/FILTERED/{sample}-barrnap-hifiasm.gff3",
        gff3_lja     = "Assemblies/{sample}/FILTERED/{sample}-barrnap-lja.gff3",
        done         = "Assemblies/{sample}/FILTERED/barrnap.done"
    log:
        barrnap_hifiasm_log  =   "Assemblies/{sample}/FILTERED/barrnap/barrnap_hifiasm.log",
        barrnap_lja_log  =   "Assemblies/{sample}/FILTERED/barrnap/barrnap_lja.log"
    conda:
        "../envs/barrnap_env.yaml"
    shell:
        """
        set -euo pipefail

        barrnap --kingdom euk --threads 8 < {input.asm_hifiasm} > {output.gff3_hifiasm} 2>&1 | tee {log.barrnap_hifiasm_log}
        barrnap --kingdom euk --threads 8 < {input.asm_lja} > {output.gff3_lja} 2>&1 | tee {log.barrnap_lja_log}

        touch {output.done}
        """

rule telosearch:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        fasta = "Assemblies/{sample}/TELOSEARCH/{sample}.hifireads.fasta",
        outdir = directory("Assemblies/{sample}/TELOSEARCH")
    params:
        threads = 4
    conda:
        "../envs/telosearch_env.yaml"
    log:
        stdout = "Assemblies/{sample}/TELOSEARCH/telosearchlr.log",
        stderr = "Assemblies/{sample}/TELOSEARCH/telosearchlr.err"
    shell:
        """
        mkdir -p {output.outdir}
        zcat {input.fastq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fasta}
        cd {output.outdir}
        TeloSearchLR.py -f {wildcards.sample}.hifireads.fasta -c {params.threads} -k 4 -K 20 -m 1 -M 100 -n 6000 &> telosearchlr.log
        """
