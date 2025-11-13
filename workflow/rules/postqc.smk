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

def get_asm_purged(wc):
    return outpath(f"Assemblies/{wc.sample}/PURGE_DUPS/purged.fa")


rule busco:
    input:
        asm_hifiasm = get_asm_hifiasm,
        asm_purged  = get_asm_purged
    output:
        busco_hifiasm_dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/{sample}-HIFI-hifiasm")),
        busco_hifiasm_done = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done"),
        busco_purged_dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-purged/{sample}-HIFI-purged")),
        busco_purged_done = outpath("Assemblies/{sample}/QC/BUSCO-purged/busco_purged.done")
    log:
        busco_hifiasm_log  = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_hifiasm.log"),
        busco_purged_log  = outpath("Assemblies/{sample}/QC/BUSCO-purged/busco_purged.log")
    params:
        cores = 32,
        hifiasm_outpath = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm"),
        purged_outpath  = outpath("Assemblies/{sample}/QC/BUSCO-purged")
    conda:
        "../envs/busco_env.yaml"
    shell:
        """
        set -euo pipefail

        mkdir -p {params.hifiasm_outpath}
        busco -i {input.asm_hifiasm} -c {params.cores} -m geno -f --auto-lineage-euk \
            --out_path {params.hifiasm_outpath} -o {wildcards.sample}-HIFI-hifiasm &> {log.busco_hifiasm_log}
        touch {output.busco_hifiasm_done}

        mkdir -p {params.purged_outpath}
        busco -i {input.asm_purged} -c {params.cores} -m geno -f --auto-lineage-euk \
            --out_path {params.purged_outpath} -o {wildcards.sample}-HIFI-purged &> {log.busco_purged_log}
        touch {output.busco_purged_done}

        """


rule quast:
    input:
        asm_hifiasm = get_asm_hifiasm,
        asm_purged  = get_asm_purged
    output:
        quast_done = outpath("Assemblies/{sample}/QC/QUAST/quast.done")
    log:
        stdout = outpath("Assemblies/{sample}/QC/QUAST/quast.log")
    params:
        cores = 32
    conda:
        "../envs/quast_env.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.quast_done})
        mkdir -p "$OUTDIR"

        echo "[INFO] Running QUAST on:"
        echo "  - {input.asm_hifiasm}"
        echo "  - {input.asm_purged}"

        quast \
            -o "$OUTDIR" \
            -t {params.cores} \
            {input.asm_hifiasm} \
            {input.asm_purged} \
            &> {log.stdout}

        touch {output.quast_done}
        """
