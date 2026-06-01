include: "../accessors.smk"

###############################################################################
# Helper accessors
###############################################################################

def get_fq1(wc):
    return sample_fq1[wc.sample]


def get_fq2(wc):
    return sample_fq2[wc.sample]


def asm_hifiasm_raw(wc):
    return outpath(f"Assemblies/{wc.sample}/HIFIASM/{wc.sample}.fasta")


def asm_lja_raw(wc):
    return outpath(f"Assemblies/{wc.sample}/LJA/assembly.fasta")


###############################################################################
# Build k-mer DB
###############################################################################

rule merqury_meryl:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2
    output:
        meryl = directory(outpath("Assemblies/{sample}/QC/merqury/{sample}.meryl"))
    log:
        out = outpath("Assemblies/{sample}/QC/merqury/meryl.log")
    conda:
        "../envs/merqury_env.yaml"
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.meryl})
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"
        meryl k=21 count {input.fq1} {input.fq2} \
            output {wildcards.sample}.meryl &>> {log.out}
        """


###############################################################################
# Merqury runs
###############################################################################

rule merqury_raw:
    input:
        meryl   = rules.merqury_meryl.output.meryl,
        hifiasm = asm_hifiasm_raw,
        lja     = asm_lja_raw
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/raw/raw.done")
    log:
        out = outpath("Assemblies/{sample}/QC/merqury/raw/raw.log")
    conda:
        "../envs/merqury_env.yaml"
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"
        merqury.sh {input.meryl} \
                   $(realpath {input.hifiasm}) \
                   $(realpath {input.lja}) raw_merqury &>> {log.out}
        touch {output.done}
        """


rule merqury_final:
    input:
        meryl   = rules.merqury_meryl.output.meryl,
        hifiasm = get_final_hifiasm_abs,
        lja     = get_final_lja_abs
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/final/final.done")
    log:
        out = outpath("Assemblies/{sample}/QC/merqury/final/final.log")
    conda:
        "../envs/merqury_env.yaml"
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"
        merqury.sh {input.meryl} \
                   $(realpath {input.hifiasm}) \
                   $(realpath {input.lja}) final_merqury &>> {log.out}
        touch {output.done}
        """


###############################################################################
# Final marker
###############################################################################
rule merqury:
    input:
        rules.merqury_raw.output.done,
        rules.merqury_final.output.done
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/merqury.done")
    shell:
        "touch {output.done}"
