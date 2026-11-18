#######################################################################
# POSTQC (BUSCO, QUAST, MAPPING)
# Uses final mito-updated assemblies for downstream QC.
#######################################################################
include: "../accessors.smk"

import shlex

def busco_options():
    options = []
    if config.get("busco_download_path"):
        options += ["--download_path", str(config["busco_download_path"])]
    if config.get("busco_offline", False):
        options.append("--offline")
    if config.get("busco_lineage"):
        options += ["-l", str(config["busco_lineage"])]
    else:
        options.append("--auto-lineage-euk")
    return " ".join(shlex.quote(value) for value in options)


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
    params:
        options=busco_options()
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done:q})"
        busco -i {input.assembly:q} -o {wildcards.sample}-hifiasm \
              -m geno -c {threads} {params.options} \
              --out_path "$(dirname {output.outdir:q})" \
              &> {log.out:q}
        touch {output.done:q}
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
    params:
        options=busco_options()
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done:q})"
        busco -i {input.assembly:q} -o {wildcards.sample}-lja \
              -m geno -c {threads} {params.options} \
              --out_path "$(dirname {output.outdir:q})" \
              &> {log.out:q}
        touch {output.done:q}
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
    params:
        options=busco_options()
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done:q})"
        busco -i {input.assembly:q} -o {wildcards.sample}-purged-hifiasm \
              -m geno -c {threads} {params.options} \
              --out_path "$(dirname {output.outdir:q})" \
              &> {log.out:q}
        touch {output.done:q}
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
    params:
        options=busco_options()
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done:q})"
        busco -i {input.assembly:q} -o {wildcards.sample}-purged-lja \
              -m geno -c {threads} {params.options} \
              --out_path "$(dirname {output.outdir:q})" \
              &> {log.out:q}
        touch {output.done:q}
        """

def selected_mapping_assembly(wc):
    if wc.asm not in enabled_assemblers(wc):
        raise ValueError(f"Assembler {wc.asm} is disabled for {wc.sample}")
    return get_final_hifiasm_abs(wc) if wc.asm == "hifiasm" else get_final_lja_abs(wc)


def selected_mapping_reports(wc):
    return [outpath(f"Assemblies/{wc.sample}/QC/Mapping/purged-{asm}.flagstat.txt")
            for asm in enabled_assemblers(wc)]


rule quast:
    input:
        raw = get_enabled_raw_assemblies,
        final = get_enabled_final_assemblies
    output:
        report = outpath("Assemblies/{sample}/QC/QUAST/report.tsv"),
        done = outpath("Assemblies/{sample}/QC/QUAST/quast.done")
    params:
        labels=lambda wc: ",".join(
            [f"raw_{asm}" for asm in enabled_assemblers(wc)] +
            [f"final_{asm}" for asm in enabled_assemblers(wc)])
    log:
        out = outpath("Assemblies/{sample}/QC/QUAST/quast.log")
    conda:
        "../envs/quast_env.yaml"
    threads: 16
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.done:q})
        mkdir -p "$OUTDIR"
        quast.py -t {threads} -o "$OUTDIR" --labels {params.labels:q} \
            {input.raw:q} {input.final:q} > {log.out:q} 2>&1
        test -s {output.report:q}
        touch {output.done:q}
        """


rule mapping_assembly:
    input:
        fq1 = get_fq1_abs,
        fq2 = get_fq2_abs,
        assembly = selected_mapping_assembly
    output:
        bam = outpath("Assemblies/{sample}/QC/Mapping/purged-{asm}.bam"),
        bai = outpath("Assemblies/{sample}/QC/Mapping/purged-{asm}.bam.bai"),
        flag = outpath("Assemblies/{sample}/QC/Mapping/purged-{asm}.flagstat.txt")
    wildcard_constraints:
        asm="hifiasm|lja"
    log:
        out = outpath("Assemblies/{sample}/QC/Mapping/{asm}.out"),
        err = outpath("Assemblies/{sample}/QC/Mapping/{asm}.err")
    conda:
        "../envs/mapping_qc.yaml"
    threads: 16
    resources:
        mem_mb=8000
    params:
        mapping_threads=lambda wc, threads: max(1, threads - 2),
        sort_threads=lambda wc, threads: 1 if threads > 1 else 0
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.bam:q})
        mkdir -p "$OUTDIR"
        WORK=$(mktemp -d "$OUTDIR/.{wildcards.asm}.mapping.XXXXXX")
        trap 'rm -rf "$WORK"' EXIT
        # Index a local link: source FASTAs and adjacent production indexes stay untouched.
        ln -s "$(realpath {input.assembly:q})" "$WORK/assembly.fasta"
        bwa index "$WORK/assembly.fasta" > {log.out:q} 2> {log.err:q}
        bwa mem -t {params.mapping_threads} "$WORK/assembly.fasta" {input.fq1:q} {input.fq2:q} \
            2>> {log.err:q} \
            | samtools sort -@ {params.sort_threads} -m 1G -T "$WORK/sort" -o {output.bam:q}
        samtools quickcheck -v {output.bam:q}
        samtools index {output.bam:q}
        samtools flagstat {output.bam:q} > {output.flag:q}
        test -s {output.flag:q}
        """


rule mapping_qc:
    input:
        selected_mapping_reports
    output:
        done = outpath("Assemblies/{sample}/QC/Mapping/mapping_qc.done")
    shell:
        "touch {output.done:q}"
