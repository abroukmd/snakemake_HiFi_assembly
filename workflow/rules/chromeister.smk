# Primary assembly comparisons. The helper isolates compressed references per job.
import os

rule chromeister_hifiasm:
    input:
        query = get_asm_hifiasm_abs,
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/HIFIASM/hifiasm.done")
    output:
        mat = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.mat"),
        png = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.mat.filt.png"),
        score = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm_score.txt"),
        done = outpath("Assemblies/{sample}/CHROMEISTER/.hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """


rule chromeister_lja:
    input:
        query = get_asm_lja_abs,
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/LJA/lja.done")
    output:
        mat = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-lja.mat"),
        png = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-lja.mat.filt.png"),
        score = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-lja_score.txt"),
        done = outpath("Assemblies/{sample}/CHROMEISTER/.lja.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-lja.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """
