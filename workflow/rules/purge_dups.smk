rule purge_dups:
    input:
        hifiasm = outpath("Assemblies/{sample}/HIFIASM/{sample}.fasta")
    output:
        split = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.fasta"),
        paf   = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.self.paf.gz"),
        bed   = outpath("Assemblies/{sample}/PURGE_DUPS/dups_{sample}.bed"),
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.fa"),
        hap   = outpath("Assemblies/{sample}/PURGE_DUPS/hap.fa"),
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

        # 1. Create cutoffs file for purge_dups
        cat <<EOF > cutoffs
low=5
mid=35
high=90
EOF

        echo ">>> Splitting contigs" >> {log.out}
        split_fa {input.hifiasm} > {output.split} 2>> {log.err}

        echo ">>> Self-alignment" >> {log.out}
        minimap2 -xasm5 -DP -t {threads} {output.split} {output.split} \
            | gzip -c > {output.paf} 2>> {log.err}

        echo ">>> Running purge_dups" >> {log.out}
        purge_dups -2 -T cutoffs -c PB.base.cov {output.paf} > {output.bed} 2>> {log.err}

        echo ">>> Extracting purged/haplotig FASTAs" >> {log.out}
        get_seqs {output.bed} {input.hifiasm} 2>> {log.err}

        # get_seqs writes purged.fa and hap.fa in CWD
        mv purged.fa {output.purge}
        mv hap.fa {output.hap}

        samtools faidx {output.purge}
        samtools faidx {output.hap}
        """
