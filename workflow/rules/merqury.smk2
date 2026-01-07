include: "../accessors.smk"

def get_fq1_abs(wc):
    return sample_fq1[wc.sample]

def get_fq2_abs(wc):
    return sample_fq2[wc.sample]

rule merqury:
    input:
        fq1     = get_fq1_abs,
        fq2     = get_fq2_abs,
        hifiasm = get_asm_hifiasm_abs,
        lja     = get_asm_lja_abs,
        purged_hifiasm = get_purged_hifiasm_or_original,
        purged_lja     = get_purged_lja_or_original
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

        cd "$OUTDIR"

        echo "[INFO] Meryl k-mer counting..." >> {log.stdout}
        meryl k=21 count {input.fq1} {input.fq2} \
            output {wildcards.sample}.meryl \
            &>> {log.stdout}

        echo "[INFO] Running Merqury..." >> {log.stdout}
        merqury.sh \
            {wildcards.sample}.meryl \
            {input.hifiasm} \
            {input.lja} \
            {input.purged_hifiasm} \
            {input.purged_lja} \
            {wildcards.sample}_merqury_out \
            &>> {log.stdout}

        touch {output.done}
        """
