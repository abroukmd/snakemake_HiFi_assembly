import os

# Define accessors for dynamic paths
def get_fastq(wc):
    return sample_paths[wc.sample]

def get_fq1(wc):
    return sample_fq1[wc.sample]

def get_fq2(wc):
    return sample_fq2[wc.sample]

def get_asm_hifiasm(wc):
    return outpath(f"Assemblies/{wc.sample}/FILTERED/{wc.sample}-HIFI-hifiasm.fasta")

def get_asm_lja(wc):
    return outpath(f"Assemblies/{wc.sample}/FILTERED/{wc.sample}-HIFI-lja.fasta")

rule busco:
    input:
        asm_hifiasm = get_asm_hifiasm,
        asm_lja     = get_asm_lja
    output:
        busco_hifiasm_dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/{sample}-HIFI-hifiasm")),
        busco_hifiasm_done = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done"),
        busco_lja_dir      = directory(outpath("Assemblies/{sample}/QC/BUSCO-lja/{sample}-HIFI-lja")),
        busco_lja_done     = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done")
    log:
        busco_hifiasm_log  = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_hifiasm.log"),
        busco_lja_log  = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.log")
    params:
        cores = 32,
        hifiasm_outpath = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm"),
        lja_outpath     = outpath("Assemblies/{sample}/QC/BUSCO-lja")
    conda:
        "../envs/busco_env.yaml"
    shell:
        """
        set -euo pipefail

        mkdir -p {params.hifiasm_outpath}
        mkdir -p {params.lja_outpath}

        busco -i {input.asm_hifiasm} -c {params.cores} -m geno -f --auto-lineage-euk \
            --out_path {params.hifiasm_outpath} -o {wildcards.sample}-HIFI-hifiasm &> {log.busco_hifiasm_log}

        touch {output.busco_hifiasm_done}

        busco -i {input.asm_lja} -c {params.cores} -m geno -f --auto-lineage-euk \
            --out_path {params.lja_outpath} -o {wildcards.sample}-HIFI-lja &> {log.busco_lja_log}

        touch {output.busco_lja_done}
        """

rule quast:
    input:
        asm_hifiasm = get_asm_hifiasm,
        asm_lja     = get_asm_lja
    output:
        quast_done = outpath("Assemblies/{sample}/QC/QUAST/quast.done")
    log:
        stdout = outpath("Assemblies/{sample}/QC/QUAST/quast.log")
    params:
        cores = 32
    conda:
        "../envs/quast_env.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.quast_done})
        quast -o $(dirname {output.quast_done}) -t 32 {input.asm_hifiasm} {input.asm_lja} &> {log.stdout}
        touch {output.quast_done}
        """
