rule barrnap:
    input:
        asm_hifiasm = get_purged_hifiasm_or_original,
        hifiasm_done = outpath("Assemblies/{sample}/HIFIASM/hifiasm.done")
    output:
        gff3_hifiasm = outpath("Assemblies/{sample}/FILTERED/{sample}-barrnap-hifiasm.gff3"),
        done         = outpath("Assemblies/{sample}/FILTERED/barrnap.done")
    log:
        barrnap_hifiasm_log  =   outpath("Assemblies/{sample}/FILTERED/barrnap/barrnap_hifiasm.log")
    conda:
        "../envs/barrnap_env.yaml"
    threads: 8
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.gff3_hifiasm:q})" "$(dirname {log.barrnap_hifiasm_log:q})"
        barrnap --kingdom euk --threads {threads} < {input.asm_hifiasm:q} \
            > {output.gff3_hifiasm:q} 2> {log.barrnap_hifiasm_log:q}
        grep -q '^##gff-version 3' {output.gff3_hifiasm:q}
        awk '!/^#/ && NF && NF != 9 {{ exit 1 }}' FS='\t' {output.gff3_hifiasm:q}

        touch {output.done}
        """

rule telosearch:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        fasta = outpath("Assemblies/{sample}/TELOSEARCH/{sample}.hifireads.fasta"),
        ranked = outpath("Assemblies/{sample}/TELOSEARCH/TeloSearchLR.k4.K20.t1000.{sample}.hifireads.rankedRepeatsTable.txt"),
        done = outpath("Assemblies/{sample}/TELOSEARCH/telosearch.done"),
        outdir = directory(outpath("Assemblies/{sample}/TELOSEARCH"))
    threads: 4
    conda:
        "../envs/telosearch_env.yaml"
    log:
        stdout = outpath("Assemblies/{sample}/logs/telosearchlr.log"),
        stderr = outpath("Assemblies/{sample}/logs/telosearchlr.err")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.outdir:q} "$(dirname {log.stdout:q})"
        if [[ {input.fastq:q} == *.gz ]]; then
            gzip -cd -- {input.fastq:q}
        else
            cat -- {input.fastq:q}
        fi | sed -n '1~4s/^@/>/p;2~4p' > {output.fasta:q}
        grep -q '^>' {output.fasta:q}
        cd {output.outdir:q}
        TeloSearchLR.py -f {wildcards.sample}.hifireads.fasta -c {threads} \
            -k 4 -K 20 -t 1000 -m 1 -M 100 -n 6000 \
            > {log.stdout:q} 2> {log.stderr:q}
        test -s {output.ranked:q}
        grep -Fq 'Script finished.' {log.stdout:q}
        touch {output.done:q}
        """
