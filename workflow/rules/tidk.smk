# Optional final-assembly telomere-repeat screen. Enabled in the main Snakefile.
import os


def tidk_final_assembly(wc):
    if wc.asm not in enabled_assemblers(wc):
        raise ValueError(f"Assembler {wc.asm} is disabled for {wc.sample}")
    return get_final_hifiasm_abs(wc) if wc.asm == "hifiasm" else get_final_lja_abs(wc)


def tidk_memory_mb(wc, input):
    # The helper retains FASTA strings and per-motif interval tables in memory.
    assembly_mb = os.path.getsize(input.assembly) / 1024**2
    return max(2048, int(12 * assembly_mb) + 512)


rule tidk:
    input:
        assembly = tidk_final_assembly,
        runner = workflow.source_path("../scripts/tidk_telomeres.py")
    output:
        results = directory(outpath("Assemblies/{sample}/QC/tidk/{asm}")),
        summary = outpath("Assemblies/{sample}/QC/tidk/{asm}/assembly_summary.tsv"),
        candidates = outpath("Assemblies/{sample}/QC/tidk/{asm}/candidate_motifs.tsv"),
        sequences = outpath("Assemblies/{sample}/QC/tidk/{asm}/telomere_counts_by_sequence.tsv"),
        arrays = outpath("Assemblies/{sample}/QC/tidk/{asm}/repeat_arrays.tsv"),
        config = outpath("Assemblies/{sample}/QC/tidk/{asm}/run_config.json"),
        receipt = outpath("Assemblies/{sample}/QC/tidk/{asm}/complete.json"),
        done = outpath("Assemblies/{sample}/QC/tidk/{asm}/tidk.done")
    log:
        out = outpath("Assemblies/{sample}/QC/tidk/logs/{asm}.stdout.log"),
        err = outpath("Assemblies/{sample}/QC/tidk/logs/{asm}.stderr.log")
    params:
        executable = str(config.get("tidk_executable", "tidk")),
        minimum = int(config.get("tidk_minimum", 5)),
        maximum = int(config.get("tidk_maximum", 30)),
        top_n = int(config.get("tidk_top_n", 10)),
        end_size = int(config.get("tidk_end_size", 10000)),
        window = int(config.get("tidk_window", 10000)),
        copies = int(config.get("tidk_copies", 5)),
        threshold = int(config.get("tidk_threshold", 20)),
        distance = float(config.get("tidk_distance", 0.05)),
        min_sequence_length = int(config.get("tidk_min_sequence_length", 0))
    wildcard_constraints:
        asm = "hifiasm|lja"
    threads: 1
    resources:
        mem_mb = tidk_memory_mb
    conda:
        "../envs/tidk_env.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log.out:q})"
        export RAYON_NUM_THREADS={threads}
        python {input.runner:q} {input.assembly:q} -o {output.results:q} \
            --tidk {params.executable:q} -m {params.minimum} -x {params.maximum} \
            -n {params.top_n} -e {params.end_size} -w {params.window} \
            -c {params.copies} -t {params.threshold} -d {params.distance} \
            -l {params.min_sequence_length} > {log.out:q} 2> {log.err:q}
        test -s {output.summary:q}
        test -s {output.receipt:q}
        touch {output.done:q}
        """
