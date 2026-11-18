include: "../accessors.smk"


def merqury_assembly(wc):
    if wc.asm not in enabled_assemblers(wc):
        raise ValueError(f"Assembler {wc.asm} is disabled for {wc.sample}")
    if wc.stage == "raw":
        return get_asm_hifiasm_abs(wc) if wc.asm == "hifiasm" else get_asm_lja_abs(wc)
    return get_final_hifiasm_abs(wc) if wc.asm == "hifiasm" else get_final_lja_abs(wc)


def merqury_reports(wc, stage):
    return [outpath(f"Assemblies/{wc.sample}/QC/merqury/{stage}/{asm}/merqury.qv")
            for asm in enabled_assemblers(wc)]


rule merqury_meryl:
    input:
        fq1 = get_fq1_abs,
        fq2 = get_fq2_abs
    output:
        meryl = directory(outpath("Assemblies/{sample}/QC/merqury/{sample}.reads.meryl"))
    log:
        out = outpath("Assemblies/{sample}/QC/merqury/meryl.log")
    conda:
        "../envs/merqury_env.yaml"
    threads: 8
    resources:
        mem_mb=16000
    params:
        k=int(config.get("merqury_k", 21))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.meryl:q})
        meryl threads={threads} memory=15 k={params.k} count {input.fq1:q} {input.fq2:q} \
            output {output.meryl:q} > {log.out:q} 2>&1
        test -s {output.meryl:q}/merylIndex
        """


rule merqury_assembly:
    input:
        meryl = rules.merqury_meryl.output.meryl,
        assembly = merqury_assembly,
        runner = workflow.source_path("../scripts/run_merqury.py")
    output:
        results = directory(outpath("Assemblies/{sample}/QC/merqury/{stage}/{asm}")),
        qv = outpath("Assemblies/{sample}/QC/merqury/{stage}/{asm}/merqury.qv"),
        completeness = outpath("Assemblies/{sample}/QC/merqury/{stage}/{asm}/completeness.stats"),
        plot = outpath("Assemblies/{sample}/QC/merqury/{stage}/{asm}/merqury.spectra-asm.ln.png")
    wildcard_constraints:
        asm="hifiasm|lja",
        stage="raw|final"
    log:
        out = outpath("Assemblies/{sample}/QC/merqury/{stage}/{asm}.log")
    conda:
        "../envs/merqury_env.yaml"
    threads: 8
    resources:
        mem_mb=16000
    params:
        label=lambda wc: f"{wc.sample}.{wc.stage}.{wc.asm}",
        k=int(config.get("merqury_k", 21))
    shell:
        r"""
        set -euo pipefail
        python {input.runner:q} --read-db {input.meryl:q} --assembly {input.assembly:q} \
            --outdir {output.results:q} --label {params.label:q} --k {params.k} \
            --threads {threads} --memory-gb 15 > {log.out:q} 2>&1
        """


rule merqury_raw:
    input:
        lambda wc: merqury_reports(wc, "raw")
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/raw/raw.done")
    shell:
        "touch {output.done:q}"


rule merqury_final:
    input:
        lambda wc: merqury_reports(wc, "final")
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/final/final.done")
    shell:
        "touch {output.done:q}"


rule merqury:
    input:
        rules.merqury_raw.output.done,
        rules.merqury_final.output.done
    output:
        done = outpath("Assemblies/{sample}/QC/merqury/merqury.done")
    shell:
        "touch {output.done:q}"
