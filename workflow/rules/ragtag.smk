rule ragtag_purged:
    input:
        ref   = lambda wildcards: SAMPLES_DF.loc[wildcards.sample, "Ref"],
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.fa")
    output:
        outdir = directory(outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag"))
    threads: 32
    conda: "../envs/ragtag.yaml"
    log:
        out = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.log"),
        err = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.err")
    shell:
        r"""
        mkdir -p {output.outdir}
        ragtag.py scaffold -t {threads} -o {output.outdir} {input.ref} {input.purge} \
            > {log.out} 2> {log.err}
        """