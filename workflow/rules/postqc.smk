#######################################################################
# POSTQC (BUSCO, QUAST, MAPPING, MULTIQC)
#
# This module now uses ONLY the unified accessors from accessors.smk:
#
#   get_asm_hifiasm_abs(sample)
#   get_asm_lja_abs(sample)
#   get_purged_hifiasm(wc)
#   get_purged_lja(wc)
#
#######################################################################

#######################################################################
# BUSCO (raw + purged assemblies)
#######################################################################
include: "../accessors.smk"

rule busco_hifiasm:
    input:
        lambda wc: get_asm_hifiasm_abs(wc)
    output:
        dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/{sample}-hifiasm")),
        done = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.log")
    conda: "../envs/busco_env.yaml"
    threads: 32
    shell:
        """
        mkdir -p $(dirname {output.done})
        busco -i {input} -o {wildcards.sample}-hifiasm \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.dir}) \
              &> {log.out}
        touch {output.done}
        """


rule busco_lja:
    input:
        lambda wc: get_asm_lja_abs(wc)
    output:
        dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-lja/{sample}-lja")),
        done = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.log")
    conda: "../envs/busco_env.yaml"
    threads: 32
    shell:
        """
        mkdir -p $(dirname {output.done})
        busco -i {input} -o {wildcards.sample}-lja \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.dir}) \
              &> {log.out}
        touch {output.done}
        """


rule busco_purged_hifiasm:
    input:
        get_purged_hifiasm_or_original
    output:
        dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/{sample}-purged-hifiasm")),
        done = outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.log")
    conda: "../envs/busco_env.yaml"
    threads: 32
    shell:
        """
        mkdir -p $(dirname {output.done})
        busco -i {input} -o {wildcards.sample}-purged-hifiasm \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.dir}) \
              &> {log.out}
        touch {output.done}
        """


rule busco_purged_lja:
    input:
        get_purged_lja_or_original
    output:
        dir  = directory(outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/{sample}-purged-lja")),
        done = outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.log")
    conda: "../envs/busco_env.yaml"
    threads: 32
    shell:
        """
        mkdir -p $(dirname {output.done})
        busco -i {input} -o {wildcards.sample}-purged-lja \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.dir}) \
              &> {log.out}
        touch {output.done}
        """


#######################################################################
# QUAST — run all 4 assemblies in one command
#######################################################################
rule quast:
    input:
        hifiasm        = lambda wc: get_asm_hifiasm_abs(wc),
        lja            = lambda wc: get_asm_lja_abs(wc),
        purged_hifiasm = get_purged_hifiasm_or_original,
        purged_lja     = get_purged_lja_or_original
    output:
        done = outpath("Assemblies/{sample}/QC/QUAST/quast.done")
    log:
        out = outpath("Assemblies/{sample}/QC/QUAST/quast.log")
    conda:
        "../envs/quast_env.yaml"
    threads: 32
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"

        echo "[INFO] Running QUAST on raw + purged assemblies" >> {log.out}

        quast -t {threads} -o "$OUTDIR" \
            {input.hifiasm} \
            {input.lja} \
            {input.purged_hifiasm} \
            {input.purged_lja} \
            &>> {log.out}

        touch {output.done}
        """


#######################################################################
# MAPPING QC — Illumina reads mapped to purged assemblies
#######################################################################
rule mapping_qc:
    input:
        fq1  = get_fq1_abs,
        fq2  = get_fq2_abs,
        phif = get_purged_hifiasm_or_original,
        plja = get_purged_lja_or_original
    output:
        bam_hifiasm  = outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.bam"),
        bam_lja      = outpath("Assemblies/{sample}/QC/Mapping/purged-lja.bam"),
        stat_hifiasm = outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.flagstat.txt"),
        stat_lja     = outpath("Assemblies/{sample}/QC/Mapping/purged-lja.flagstat.txt"),
        done         = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done")
    log:
        out = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.log"),
        err = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.err")
    conda:
        "../envs/mapping_qc.yaml"
    threads: 16
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"

        mapit() {{
            ASM=$1
            PREFIX=$2

            bwa index "$ASM" >> {log.out} 2>> {log.err} || true

            bwa mem -t {threads} "$ASM" {input.fq1} {input.fq2} \
                | samtools sort -@ {threads} -o "$OUTDIR/$PREFIX.bam"

            samtools index "$OUTDIR/$PREFIX.bam"
            samtools flagstat "$OUTDIR/$PREFIX.bam" > "$OUTDIR/$PREFIX.flagstat.txt"
        }}

        mapit "{input.phif}" "purged-hifiasm"
        mapit "{input.plja}" "purged-lja"

        touch {output.done}
        """


#######################################################################
# MULTIQC
#######################################################################
rule multiqc_sample:
    input:
        quast      = outpath("Assemblies/{sample}/QC/QUAST/quast.done"),
        busco_hif  = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done"),
        busco_lja  = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done"),
        busco_phif = outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.done"),
        busco_plja = outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.done"),
        mapping    = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done"),
        merqury    = outpath("Assemblies/{sample}/QC/merqury/merqury.done")
    output:
        report = outpath("Assemblies/{sample}/QC/multiqc/multiqc_report.html"),
        done   = outpath("Assemblies/{sample}/QC/multiqc/multiqc.done")
    log:
        out = outpath("Assemblies/{sample}/QC/multiqc/multiqc.log"),
        err = outpath("Assemblies/{sample}/QC/multiqc/multiqc.err")
    conda:
        "../envs/multiqc_env.yaml"
    threads: 12
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.report})
        mkdir -p "$OUTDIR"

        SAMPLE_QC_DIR=$(realpath {outpath("Assemblies/{sample}/QC")})

        echo "[INFO] Running MultiQC on QC folder" >> {log.out}

        multiqc "$SAMPLE_QC_DIR" \
            --outdir "$OUTDIR" \
            --filename multiqc_report.html \
            >> {log.out} 2>> {log.err}

        touch {output.done}
        """
