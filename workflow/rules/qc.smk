rule nanoplot:
    input:
        lambda wildcards: sample_paths[wildcards.sample]
    log:
        stdout = outpath("Assemblies/{sample}/QC/nanoplot/nanoplot.log"),
        stderr = outpath("Assemblies/{sample}/QC/nanoplot/nanoplot.err")
    output:
        nanoplot_dir = directory(outpath("Assemblies/{sample}/QC/nanoplot"))
    conda:
        "../envs/nanoplot_env.yaml"
    params:
        cores=12,
        mem="10gb",
        time="01:00:00"
    shell:
        "NanoPlot -o {output} -t {params.cores} --fastq {input} --maxlength 40000 --plots dot"

#rule jellyfish_histo:
#    input:
#        fastq = lambda wildcards: sample_paths[wildcards.sample]
#    log:
#        stdout = outpath("Assemblies/{sample}/QC/genomescope/genomescope.log"),
#        stderr = outpath("Assemblies/{sample}/QC/genomescope/genomescope.err")
#    output:
#        histo = outpath("Assemblies/{sample}/QC/genomescope/{sample}_21mer.histo")
#    params:
#        k = 21,
#        outdir = outpath("Assemblies/{sample}/QC/genomescope")
#    conda:
#        "../envs/genomescope_env.yaml"
#    shell:
#        """
#        mkdir -p {params.outdir}
#        zcat {input.fastq} | jellyfish count -C -m {params.k} -s 100M -t 10 -o {params.outdir}/{wildcards.sample}_{params.k}mer.jf /dev/fd/0
#        jellyfish histo -t 10 {params.outdir}/{wildcards.sample}_{params.k}mer.jf > {output.histo}
#        """

#rule genomescope:
#    input:
#        histo = outpath("Assemblies/{sample}/QC/genomescope/{sample}_21mer.histo")
#    output:
#        model = outpath("Assemblies/{sample}/QC/genomescope/ploidy_{ploidy}/model.txt")
#    params:
#        k = 21,
#        outdir = outpath("Assemblies/{sample}/QC/genomescope/ploidy_{ploidy}")
#    conda:
#        "../envs/genomescope_env.yaml"
#    log:
#        outpath("Assemblies/{sample}/QC/genomescope/genomescope-ploidy_{ploidy}.log")
#    shell:
#        """
#        mkdir -p {params.outdir}
#        Rscript DB/genomescope2/genomescope.R \
#            -i {input.histo} \
#            -o {params.outdir} \
#            -p {wildcards.ploidy} \
#            -k {params.k} > {log} 2>&1
#        """
