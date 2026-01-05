###############################################################
# PLOIDY PLOT — USE GLOBAL ACCESSORS
###############################################################
include: "../accessors.smk"

def get_purged_for_asm(wc):
    """Return purged assembly based on wildcard asm."""
    if wc.asm == "hifiasm":
        return get_purged_hifiasm_or_original(wc)
    elif wc.asm == "lja":
        return get_purged_lja_or_original(wc)
    else:
        raise ValueError(f"Unknown asm type: {wc.asm}")


###############################################################
# Build FAI
###############################################################
rule build_fai:
    input:
        fa = get_purged_for_asm
    output:
        fai = outpath("Assemblies/{sample}/Ploidy_plots/{asm}.fai")
    conda:
        "../envs/mapping_qc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fai})
        cp {input.fa} $(dirname {output.fai})/{wildcards.asm}.fa
        samtools faidx $(dirname {output.fai})/{wildcards.asm}.fa
        mv $(dirname {output.fai})/{wildcards.asm}.fa.fai {output.fai}
        """


###############################################################
# Compute depth from BAM → depth file under Ploidy_plots
###############################################################
rule compute_depth:
    input:
        bam = lambda wc: outpath(f"Assemblies/{wc.sample}/QC/Mapping/purged-{wc.asm}.bam"),
        fai = rules.build_fai.output.fai
    output:
        depth = outpath("Assemblies/{sample}/Ploidy_plots/{asm}.depth")
    conda:
        "../envs/mapping_qc.yaml"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.depth})
        samtools depth -@ {threads} {input.bam} > {output.depth}
        """


###############################################################
# Ploidy Plot PDF
###############################################################
rule ploidy_plot:
    input:
        fai   = rules.build_fai.output.fai,
        depth = rules.compute_depth.output.depth
    output:
        pdf = outpath("Assemblies/{sample}/Ploidy_plots/{asm}_ploidy_plot.pdf")
    conda:
        "../envs/ploidy_env.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.pdf})

        Rscript {workflow.basedir}/scripts/ploidy_plot.R \
            "{wildcards.sample}_{wildcards.asm}" \
            "{input.fai}" \
            "{input.depth}" \
            "{output.pdf}"
        """
    
    
