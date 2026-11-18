include: "../accessors.smk"

rule integrate_mito_hifiasm:
    input:
        assembly = get_purged_hifiasm_or_original,
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        helper = workflow.source_path("../scripts/integrate_mito.py")
    output:
        final = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/{sample}.purged.final.fasta"),
        paf = outpath("Assemblies/{sample}/OATK/mito_integration/hifiasm/mito_vs_assembly.paf"),
        removed = outpath("Assemblies/{sample}/OATK/mito_integration/hifiasm/contigs_to_remove.txt"),
        report = outpath("Assemblies/{sample}/OATK/mito_integration/hifiasm/integration.tsv")
    params:
        min_identity = config.get("mito_min_identity", 0.95),
        min_qcov = config.get("mito_min_query_coverage", 0.80),
        min_tcov = config.get("mito_min_target_coverage", 0.95)
    log:
        outpath("Assemblies/{sample}/OATK/mito_integration/hifiasm/integration.log")
    threads: 4
    conda:
        "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.paf:q})" "$(dirname {output.final:q})"
        minimap2 -x asm5 -t {threads} {input.assembly:q} {input.mito:q} > {output.paf:q} 2> {log:q}
        python {input.helper:q} --assembly {input.assembly:q} --mito {input.mito:q} \
            --paf {output.paf:q} --sample {wildcards.sample:q} --output {output.final:q} \
            --removed {output.removed:q} --report {output.report:q} \
            --min-identity {params.min_identity} --min-query-coverage {params.min_qcov} \
            --min-target-coverage {params.min_tcov} >> {log:q} 2>&1
        test -s {output.final:q}
        """


rule integrate_mito_lja:
    input:
        assembly = get_purged_lja_or_original,
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        helper = workflow.source_path("../scripts/integrate_mito.py")
    output:
        final = outpath("{sample}_FINAL_ASSEMBLIES/lja/{sample}.purged.final.fasta"),
        paf = outpath("Assemblies/{sample}/OATK/mito_integration/lja/mito_vs_assembly.paf"),
        removed = outpath("Assemblies/{sample}/OATK/mito_integration/lja/contigs_to_remove.txt"),
        report = outpath("Assemblies/{sample}/OATK/mito_integration/lja/integration.tsv")
    params:
        min_identity = config.get("mito_min_identity", 0.95),
        min_qcov = config.get("mito_min_query_coverage", 0.80),
        min_tcov = config.get("mito_min_target_coverage", 0.95)
    log:
        outpath("Assemblies/{sample}/OATK/mito_integration/lja/integration.log")
    threads: 4
    conda:
        "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.paf:q})" "$(dirname {output.final:q})"
        minimap2 -x asm5 -t {threads} {input.assembly:q} {input.mito:q} > {output.paf:q} 2> {log:q}
        python {input.helper:q} --assembly {input.assembly:q} --mito {input.mito:q} \
            --paf {output.paf:q} --sample {wildcards.sample:q} --output {output.final:q} \
            --removed {output.removed:q} --report {output.report:q} \
            --min-identity {params.min_identity} --min-query-coverage {params.min_qcov} \
            --min-target-coverage {params.min_tcov} >> {log:q} 2>&1
        test -s {output.final:q}
        """
