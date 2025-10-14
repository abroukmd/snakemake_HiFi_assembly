# rules/oatk.smk

project_dir = config["project_dir"]

rule oatk:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        gfa = "Assemblies/{sample}/OATK/{sample}.mito.gfa",
        ctg = "Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta",
        done = "Assemblies/{sample}/OATK/.done"
    params:
        outprefix = "Assemblies/{sample}/OATK/{sample}",
        db = config["oatk_db"]
    conda:
        "../envs/oatk_env.yaml"
    log:
        "Assemblies/{sample}/OATK/oatk.log"
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        oatk -o {params.outprefix} -t 32 -m {params.db} {input.fastq} > {log} 2>&1
        touch {output.done}
        """

#rule bandage:
#    input:
#        gfa = "Assemblies/{sample}/OATK/{sample}.mito.gfa"
#    output:
#        image = "Assemblies/{sample}/OATK/{sample}.bandage.png"
#    conda:
#        "../envs/bandage_env.yaml"
#    shell:
#        "Bandage image {input.gfa} {output.image}"
