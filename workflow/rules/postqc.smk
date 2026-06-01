#######################################################################
# POSTQC (BUSCO, QUAST, MAPPING)
# Uses final mito-updated assemblies for downstream QC.
#######################################################################
include: "../accessors.smk"

rule busco_hifiasm:
    input:
        assembly = get_asm_hifiasm_abs
    output:
        outdir = directory(outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/{sample}-hifiasm")),
        done   = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-hifiasm/busco_hifiasm.log")
    conda:
        "../envs/busco_env.yaml"
    threads: 32
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.done})
        busco -i {input.assembly} -o {wildcards.sample}-hifiasm \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.outdir}) \
              &> {log.out}
        touch {output.done}
        """

rule busco_lja:
    input:
        assembly = get_asm_lja_abs
    output:
        outdir = directory(outpath("Assemblies/{sample}/QC/BUSCO-lja/{sample}-lja")),
        done   = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-lja/busco_lja.log")
    conda:
        "../envs/busco_env.yaml"
    threads: 32
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.done})
        busco -i {input.assembly} -o {wildcards.sample}-lja \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.outdir}) \
              &> {log.out}
        touch {output.done}
        """

rule busco_purged_hifiasm:
    input:
        assembly = get_final_hifiasm_abs
    output:
        outdir = directory(outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/{sample}-purged-hifiasm")),
        done   = outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-purged-hifiasm/busco_purged_hifiasm.log")
    conda:
        "../envs/busco_env.yaml"
    threads: 32
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.done})
        busco -i {input.assembly} -o {wildcards.sample}-purged-hifiasm \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.outdir}) \
              &> {log.out}
        touch {output.done}
        """

rule busco_purged_lja:
    input:
        assembly = get_final_lja_abs
    output:
        outdir = directory(outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/{sample}-purged-lja")),
        done   = outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.done")
    log:
        out = outpath("Assemblies/{sample}/QC/BUSCO-purged-lja/busco_purged_lja.log")
    conda:
        "../envs/busco_env.yaml"
    threads: 32
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.done})
        busco -i {input.assembly} -o {wildcards.sample}-purged-lja \
              -m geno -c {threads} --auto-lineage-euk \
              --out_path $(dirname {output.outdir}) \
              &> {log.out}
        touch {output.done}
        """

rule quast:
    input:
        hifiasm   = get_asm_hifiasm_abs,
        lja       = get_asm_lja_abs,
        final_hif = get_final_hifiasm_abs,
        final_lja = get_final_lja_abs
    output:
        done = outpath("Assemblies/{sample}/QC/QUAST/quast.done")
    log:
        out = outpath("Assemblies/{sample}/QC/QUAST/quast.log")
    conda:
        "../envs/quast_env.yaml"
    threads: 16
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.done})
        mkdir -p "$OUTDIR"

        quast.py -t {threads} -o "$OUTDIR" \
            {input.hifiasm} \
            {input.lja} \
            {input.final_hif} \
            {input.final_lja} \
            &> {log.out}

        touch {output.done}
        """

rule mapping_qc:
    input:
        fq1 = get_fq1_abs,
        fq2 = get_fq2_abs,
        hif = get_final_hifiasm_abs,
        lja = get_final_lja_abs
    output:
        bam_hif  = outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.bam"),
        bai_hif  = outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.bam.bai"),
        flag_hif = outpath("Assemblies/{sample}/QC/Mapping/purged-hifiasm.flagstat.txt"),
        bam_lja  = outpath("Assemblies/{sample}/QC/Mapping/purged-lja.bam"),
        bai_lja  = outpath("Assemblies/{sample}/QC/Mapping/purged-lja.bam.bai"),
        flag_lja = outpath("Assemblies/{sample}/QC/Mapping/purged-lja.flagstat.txt"),
        done     = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done")
    log:
        out = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.out"),
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

        echo "[INFO] Mapping reads to mito-updated final assemblies" >> {log.out}
        mapit "{input.hif}" "purged-hifiasm"
        mapit "{input.lja}" "purged-lja"

        touch {output.done}
        """
