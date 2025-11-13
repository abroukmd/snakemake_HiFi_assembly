import os

def normalize_sample(wc):
    return wc.sample.split("/")[-1]

def get_asm_hifiasm_abs(wc):
    sample = normalize_sample(wc)
    return outpath(f"Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta")

def get_asm_purged_abs(wc):
    sample = normalize_sample(wc)
    return outpath(f"Assemblies/{sample}/PURGE_DUPS/purged.fa")

def get_fq1_abs(wc):
    return sample_fq1[normalize_sample(wc)]

def get_fq2_abs(wc):
    return sample_fq2[normalize_sample(wc)]


rule merqury:
    input:
        fq1     = get_fq1_abs,
        fq2     = get_fq2_abs,
        hifiasm = get_asm_hifiasm_abs,
        purged  = get_asm_purged_abs
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/merqury.done")
    log:
        stdout = outpath("Assemblies/{sample}/QC/merqury/merqury.log")
    conda:
        "../envs/merqury_env.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"

        SAMPLE="{wildcards.sample}"

        # Absolute paths
        HIFIASM=$(realpath {input.hifiasm})
        PURGED=$(realpath {input.purged})

        cd "$OUTDIR"

        echo "[INFO] Running Meryl k-mer counting..." >> merqury.log
        meryl k=21 count {input.fq1} {input.fq2} output ${{SAMPLE}}.meryl &>> merqury.log

        echo "[INFO] Running Merqury on: $HIFIASM and $PURGED" >> merqury.log
        merqury.sh ${{SAMPLE}}.meryl $HIFIASM $PURGED ${{SAMPLE}}_merqury_out &>> merqury.log

        touch merqury.done
        """
