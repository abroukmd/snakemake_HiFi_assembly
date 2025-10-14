import os

def normalize_sample(wc):
    return wc.sample.split("/")[-1]

def get_asm_hifiasm_abs(wc):
    sample = normalize_sample(wc)
    return os.path.join(wc.sample, "FILTERED", f"{sample}-HIFI-hifiasm.fasta")

def get_asm_lja_abs(wc):
    sample = normalize_sample(wc)
    return os.path.join(wc.sample, "FILTERED", f"{sample}-HIFI-lja.fasta")

def get_fq1_abs(wc):
    sample = normalize_sample(wc)
    return sample_fq1[sample]

def get_fq2_abs(wc):
    sample = normalize_sample(wc)
    return sample_fq2[sample]


rule merqury:
    input:
        fq1 = get_fq1_abs,
        fq2 = get_fq2_abs,
        hifiasm = get_asm_hifiasm_abs,
        lja = get_asm_lja_abs
    output:
        done = "{sample}/QC/merqury/merqury.done"
    log:
        stdout = "{sample}/QC/merqury/merqury.log"
    conda:
        "../envs/merqury_env.yaml"
    shell:
        """
    set -euo pipefail

    mkdir -p {wildcards.sample}/QC/merqury/

    SAMPLE_NAME=$(basename {wildcards.sample})

    # Get absolute paths before cd
    HIFIASM=$(realpath {input.hifiasm})
    LJA=$(realpath {input.lja})

    cd {wildcards.sample}/QC/merqury/

    # Build k-mer DB (relative output)
    meryl k=21 count {input.fq1} {input.fq2} output ${{SAMPLE_NAME}}.meryl &>> merqury.log

    # Run Merqury with absolute paths
    merqury.sh ${{SAMPLE_NAME}}.meryl $HIFIASM $LJA ${{SAMPLE_NAME}}_merqury_out &>> merqury.log

    # Mark completion
    touch merqury.done
        """
