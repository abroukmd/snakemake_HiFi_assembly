# Optional RagTag comparisons; included only when explicitly enabled.
import os

rule chromeister_ragtag_purged:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.hifiasm.reforder.fasta"),
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    output:
        mat = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.mat"),
        png = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.mat.filt.png"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag_score.txt"),
        done = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.purged_ragtag.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    wildcard_constraints:
        sample = "|".join(RAGTAG_SAMPLES) or "(?!)"
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """


rule chromeister_ragtag_lja:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.lja.reforder.fasta"),
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}_ragtag_lja.done")
    output:
        mat = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_lja_ragtag.mat"),
        png = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_lja_ragtag.mat.filt.png"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_lja_ragtag_score.txt"),
        done = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.purged_lja_ragtag.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_lja_ragtag.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_lja_ragtag.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    wildcard_constraints:
        sample = "|".join(RAGTAG_SAMPLES) or "(?!)"
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """
