# rules/oatk.smk

project_dir = config["project_dir"]

rule oatk:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        gfa = outpath("Assemblies/{sample}/OATK/{sample}.mito.gfa"),
        ctg = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        done = outpath("Assemblies/{sample}/OATK/.done")
    params:
        outprefix = outpath("Assemblies/{sample}/OATK/{sample}"),
        db = config["oatk_db"]
    conda:
        "../envs/oatk_env.yaml"
    log:
        outpath("Assemblies/{sample}/OATK/oatk.log")
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        oatk -o {params.outprefix} -t 32 -m {params.db} {input.fastq} > {log} 2>&1
        touch {output.done}
        """

#rule bandage:
#    input:
#        gfa = outpath("Assemblies/{sample}/OATK/{sample}.mito.gfa")
#    output:
#        image = outpath("Assemblies/{sample}/OATK/{sample}.bandage.png")
#    conda:
#        "../envs/bandage_env.yaml"
#    shell:
#        "Bandage image {input.gfa} {output.image}"
