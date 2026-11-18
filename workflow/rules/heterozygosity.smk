# Optional standardized diploid-model estimate, conditional on the final assembly.
import json


def heterozygosity_final_assembly(wc):
    if wc.asm not in enabled_assemblers(wc):
        raise ValueError(f"Assembler {wc.asm} is disabled for {wc.sample}")
    return get_final_hifiasm_abs(wc) if wc.asm == "hifiasm" else get_final_lja_abs(wc)


def heterozygosity_extra_ids(wc):
    by_assembler = config.get("heterozygosity_exclude_contigs", {}).get(wc.sample, {})
    if not isinstance(by_assembler, dict):
        raise ValueError("heterozygosity_exclude_contigs must map sample -> assembler -> list of exact IDs")
    identifiers = by_assembler.get(wc.asm, [])
    if not isinstance(identifiers, list) or any(not isinstance(value, str) for value in identifiers):
        raise ValueError("Each assembler's heterozygosity exclusion must be a list of exact contig IDs")
    return json.dumps(identifiers)


rule heterozygosity:
    input:
        assembly=heterozygosity_final_assembly,
        bam=outpath("Assemblies/{sample}/QC/Mapping/purged-{asm}.bam"),
        bai=outpath("Assemblies/{sample}/QC/Mapping/purged-{asm}.bam.bai"),
        mito=outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        helper=workflow.source_path("../scripts/estimate_heterozygosity.py")
    output:
        summary=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/summary.json"),
        table=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/summary.tsv"),
        callable=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/callable.bed"),
        nuclear=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/nuclear_acgt.bed"),
        excluded=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/excluded_contigs.txt"),
        variants=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/heterozygous_snps.vcf.gz"),
        variant_index=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/heterozygous_snps.vcf.gz.tbi"),
        sites=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/heterozygous_snps.tsv.gz"),
        balance=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/allele_balance.tsv"),
        receipt=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/receipt.json"),
        done=outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/heterozygosity.done")
    wildcard_constraints:
        asm="hifiasm|lja"
    log:
        outpath("Assemblies/{sample}/QC/heterozygosity/{asm}/run.log")
    conda:
        "../envs/heterozygosity_env.yaml"
    threads: 4
    resources:
        mem_mb=8000
    params:
        outdir=lambda wc: outpath(f"Assemblies/{wc.sample}/QC/heterozygosity/{wc.asm}"),
        extra_ids=heterozygosity_extra_ids,
        min_callable_fraction=float(config.get("heterozygosity_min_callable_fraction", 0.5))
    shell:
        r"""
        set -euo pipefail
        python {input.helper:q} --reference {input.assembly:q} --bam {input.bam:q} \
            --bam-index {input.bai:q} --mito-source {input.mito:q} \
            --sample {wildcards.sample:q} --assembler {wildcards.asm:q} \
            --outdir {params.outdir:q} --extra-excluded-ids-json {params.extra_ids:q} \
            --min-callable-fraction {params.min_callable_fraction} --threads {threads} \
            > {log:q} 2>&1
        """
