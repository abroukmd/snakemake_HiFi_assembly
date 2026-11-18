# rules/oatk.smk

project_dir = config["project_dir"]

rule oatk:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample],
        db = config["oatk_db"]
    output:
        gfa = outpath("Assemblies/{sample}/OATK/{sample}.mito.gfa"),
        ctg = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        done = outpath("Assemblies/{sample}/OATK/.done")
    params:
        outprefix = outpath("Assemblies/{sample}/OATK/{sample}")
    threads: 32
    conda:
        "../envs/oatk_env.yaml"
    log:
        outpath("Assemblies/{sample}/OATK/oatk.log")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {params.outprefix:q})"
        oatk -o {params.outprefix:q} -t {threads} -m {input.db:q} {input.fastq:q} > {log:q} 2>&1
        test -s {output.gfa:q}
        grep -q '^>' {output.ctg:q}
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
