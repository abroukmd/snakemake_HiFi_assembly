rule purge_dups:
    input:
        hifiasm = outpath("Assemblies/{sample}/HIFIASM/{sample}.fasta")
    output:
        split = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.fasta"),
        paf   = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.self.paf.gz"),
        bed   = outpath("Assemblies/{sample}/PURGE_DUPS/dups_{sample}.bed"),
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.hifiasm.fa"),
        hap   = outpath("Assemblies/{sample}/PURGE_DUPS/hap.hifiasm.fa"),
    log:
        out = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}-purge_dups.log"),
        err = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}-purge_dups.err")
    conda:
        "../envs/purge_dups.yaml"
    threads: 24
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.purge})"

        cat <<EOF > cutoffs
low=5
mid=35
high=90
EOF

        split_fa {input.hifiasm} > {output.split} 2>> {log.err}

        minimap2 -xasm5 -DP -t {threads} {output.split} {output.split} \
            | gzip -c > {output.paf} 2>> {log.err}

        purge_dups -2 -T cutoffs -c PB.base.cov {output.paf} > {output.bed} 2>> {log.err}

        get_seqs {output.bed} {input.hifiasm} 2>> {log.err}

        mv purged.fa {output.purge}
        mv hap.fa    {output.hap}

        samtools faidx {output.purge}
        samtools faidx {output.hap}
        """


rule purge_dups_lja:
    input:
        lja = outpath("Assemblies/{sample}/LJA/assembly.fasta")
    output:
        split = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.fasta"),
        paf   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.self.paf.gz"),
        bed   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/dups_{sample}.bed"),
        purge = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.lja.fa"),
        hap   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/hap.lja.fa"),
        done  = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purge_lja.done")
    log:
        out = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}-purge_dups_lja.log"),
        err = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}-purge_dups_lja.err")
    conda:
        "../envs/purge_dups.yaml"
    threads: 24
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.purge})"

        cat <<EOF > cutoffs
low=5
mid=35
high=90
EOF

        split_fa {input.lja} > {output.split} 2>> {log.err}

        minimap2 -xasm5 -DP -t {threads} {output.split} {output.split} \
            | gzip -c > {output.paf} 2>> {log.err}

        purge_dups -2 -T cutoffs -c PB.base.cov {output.paf} > {output.bed} 2>> {log.err}

        if [[ ! -s {output.bed} ]]; then
            cp {input.lja} {output.purge}
            echo -e ">hap_empty\nN" > {output.hap}
        else
            get_seqs {output.bed} {input.lja} 2>> {log.err}

            if [[ ! -s purged.fa ]]; then
                cp {input.lja} {output.purge}
                echo -e ">hap_empty\nN" > {output.hap}
            else
                mv purged.fa {output.purge}
                mv hap.fa    {output.hap}
            fi
        fi

        samtools faidx {output.purge}
        samtools faidx {output.hap}

        touch {output.done}
        """
