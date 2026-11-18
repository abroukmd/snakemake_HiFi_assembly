# Mitochondria-integrated final assembly comparisons and optional haplotig plots.
import os

rule chromeister_final_hifiasm:
    input:
        query = get_final_hifiasm_abs,
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh")
    output:
        mat = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.mat"),
        png = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.mat.filt.png"),
        score = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm_score.txt"),
        done = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/.final_hifiasm.done")
    log:
        out = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.log"),
        err = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """


rule chromeister_final_lja:
    input:
        query = get_final_lja_abs,
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh")
    output:
        mat = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.mat"),
        png = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.mat.filt.png"),
        score = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja_score.txt"),
        done = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/.final_lja.done")
    log:
        out = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.log"),
        err = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 0
        """


rule chromeister_final_hap_hifiasm:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/{sample}.hap.final.fasta"),
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    output:
        mat = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.mat"),
        png = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.mat.filt.png"),
        score = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm_score.txt"),
        done = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/.final_hap_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.log"),
        err = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 1
        """


rule chromeister_final_hap_lja:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/{sample}.hap.final.fasta"),
        ref = lambda wc: sample_refs[wc.sample],
        runner = os.path.join(workflow.basedir, "scripts", "run_chromeister.sh"),
        assembly_done = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}_ragtag_lja.done")
    output:
        mat = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.mat"),
        png = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.mat.filt.png"),
        score = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja_score.txt"),
        done = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/.final_hap_lja.done")
    log:
        out = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.log"),
        err = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        bash {input.runner:q} {input.query:q} {input.ref:q} \
            {output.mat:q} {output.png:q} {output.score:q} {output.done:q} \
            {log.out:q} {log.err:q} 1
        """
