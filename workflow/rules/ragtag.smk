import pandas as pd

# Load sample sheet
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample")

# Samples requiring RagTag
RAGTAG_SAMPLES = samples_df.query("ragtag == 'Y'").index.tolist()


# -----------------------------------------------------
# RagTag on purge_dups "purged.fa"
# -----------------------------------------------------
rule ragtag_purged:
    input:
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.fa")
    output:
        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
        fai   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta.fai"),
        done  = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    threads: 32
    conda: "../envs/ragtag.yaml"
    wildcard_constraints:
        sample="|".join(RAGTAG_SAMPLES)
    log:
        out = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.log"),
        err = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.err")
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})

        ragtag.py scaffold \
            -t {threads} \
            -o $(dirname {output.fasta}) \
            {input.ref} {input.purge} \
            > {log.out} 2> {log.err}

        # Index output
        samtools faidx {output.fasta}

        touch {output.done}
        """
