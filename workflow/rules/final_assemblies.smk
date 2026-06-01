include: "../accessors.smk"

rule integrate_mito_hifiasm:
    input:
        assembly = get_purged_hifiasm_or_original,
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta")
    output:
        final = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/{sample}.purged.final.fasta")
    params:
        workdir = outpath("Assemblies/{sample}/OATK/mito_integration/hifiasm"),
        min_identity = 0.95,
        min_qcov = 0.80
    threads: 4
    conda:
        "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}
        mkdir -p $(dirname {output.final})

        minimap2 -x asm5 -t {threads} {input.assembly} {input.mito} > {params.workdir}/mito_vs_assembly.paf

        awk -v min_ident={params.min_identity} -v min_qcov={params.min_qcov} 'BEGIN{{OFS="\t"}} {{
            qname=$1; qlen=$2; qstart=$3; qend=$4;
            matches=$10; alnlen=$11;
            if (alnlen > 0 && qlen > 0) {{
                identity=matches/alnlen;
                qcov=(qend-qstart)/qlen;
                if (identity >= min_ident && qcov >= min_qcov) print qname;
            }}
        }}' {params.workdir}/mito_vs_assembly.paf | sort -u > {params.workdir}/contigs_to_remove.txt

        awk 'BEGIN {{
                 while ((getline line < "{params.workdir}/contigs_to_remove.txt") > 0) {{
                     remove[line]=1
                 }}
                 close("{params.workdir}/contigs_to_remove.txt")
                 keep=1
             }}
             /^>/ {{
                 header=$0
                 name=substr($1,2)
                 sub(/[[:space:]].*$/, "", name)
                 keep = !(name in remove)
             }}
             keep {{ print }}' {input.assembly} > {params.workdir}/filtered.fa

        awk 'BEGIN{{printed=0}}
             /^>/ {{
                 if (!printed) {{
                     print ">{wildcards.sample}_mito"
                     printed=1
                     next
                 }}
             }}
             {{ print }}' {input.mito} > {params.workdir}/mito_renamed.fa

        cat {params.workdir}/filtered.fa {params.workdir}/mito_renamed.fa > {output.final}
        """


rule integrate_mito_lja:
    input:
        assembly = get_purged_lja_or_original,
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta")
    output:
        final = outpath("{sample}_FINAL_ASSEMBLIES/lja/{sample}.purged.final.fasta")
    params:
        workdir = outpath("Assemblies/{sample}/OATK/mito_integration/lja"),
        min_identity = 0.95,
        min_qcov = 0.80
    threads: 4
    conda:
        "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.workdir}
        mkdir -p $(dirname {output.final})

        minimap2 -x asm5 -t {threads} {input.assembly} {input.mito} > {params.workdir}/mito_vs_assembly.paf

        awk -v min_ident={params.min_identity} -v min_qcov={params.min_qcov} 'BEGIN{{OFS="\t"}} {{
            qname=$1; qlen=$2; qstart=$3; qend=$4;
            matches=$10; alnlen=$11;
            if (alnlen > 0 && qlen > 0) {{
                identity=matches/alnlen;
                qcov=(qend-qstart)/qlen;
                if (identity >= min_ident && qcov >= min_qcov) print qname;
            }}
        }}' {params.workdir}/mito_vs_assembly.paf | sort -u > {params.workdir}/contigs_to_remove.txt

        awk 'BEGIN {{
                 while ((getline line < "{params.workdir}/contigs_to_remove.txt") > 0) {{
                     remove[line]=1
                 }}
                 close("{params.workdir}/contigs_to_remove.txt")
                 keep=1
             }}
             /^>/ {{
                 header=$0
                 name=substr($1,2)
                 sub(/[[:space:]].*$/, "", name)
                 keep = !(name in remove)
             }}
             keep {{ print }}' {input.assembly} > {params.workdir}/filtered.fa

        awk 'BEGIN{{printed=0}}
             /^>/ {{
                 if (!printed) {{
                     print ">{wildcards.sample}_mito"
                     printed=1
                     next
                 }}
             }}
             {{ print }}' {input.mito} > {params.workdir}/mito_renamed.fa

        cat {params.workdir}/filtered.fa {params.workdir}/mito_renamed.fa > {output.final}
        """
