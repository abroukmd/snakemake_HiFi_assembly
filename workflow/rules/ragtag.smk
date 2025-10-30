rule ragtag_scaffold:
    """
    Run RagTag scaffolding for both hifiasm and lja filtered assemblies.
    """
    input:
        hifiasm = lambda wc: outpath(f"Assemblies/{wc.sample}/FILTERED/{wc.sample}-HIFI-hifiasm.fasta"),
        lja = lambda wc: outpath(f"Assemblies/{wc.sample}/FILTERED/{wc.sample}-HIFI-lja.fasta"),
        reference = lambda wc: reference_from_tsv(wc.sample)
    output:
        hifiasm_scaf = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.fasta"),
        lja_scaf = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_lja.fasta"),
        hifiasm_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.done"),
        lja_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_lja.done")
    threads: 8
    conda:
        "../envs/ragtag.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.hifiasm_scaf})

        echo ">>> Running RagTag for {wildcards.sample} <<<"

        # Scaffold Hifiasm assembly
        ragtag.py scaffold {input.reference} {input.hifiasm} \
            -o $(dirname {output.hifiasm_scaf})/hifiasm_tmp
        cp $(dirname {output.hifiasm_scaf})/hifiasm_tmp/ragtag.scaffold.fasta \
           {output.hifiasm_scaf}
        touch {output.hifiasm_done}

        # Scaffold LJA assembly
        ragtag.py scaffold {input.reference} {input.lja} \
            -o $(dirname {output.lja_scaf})/lja_tmp
        cp $(dirname {output.lja_scaf})/lja_tmp/ragtag.scaffold.fasta \
           {output.lja_scaf}
        touch {output.lja_done}

        echo ">>> RagTag completed for {wildcards.sample} <<<"
        """
