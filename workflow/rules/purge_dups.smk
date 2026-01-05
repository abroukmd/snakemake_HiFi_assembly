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

#############################################
# purge_dups for LJA assemblies
#############################################
#rule clean_lja_headers:
#    input:
#        lja = outpath("Assemblies/{sample}/LJA/assembly.fasta")
#    output:
#        clean = outpath("Assemblies/{sample}/LJA/assembly.clean.fasta")
#    shell:
#        r'''
#        awk '
#            /^>/ {{ printf(">contig_%d\n", ++c); next }}
#            {{ print }}
#        ' {input.lja} > {output.clean}
#        '''

rule purge_dups_lja:
    input:
        lja = outpath("Assemblies/{sample}/LJA/assembly.fasta")
    output:
        split = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.fasta"),
        paf   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.self.paf.gz"),
        bed   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/dups_{sample}.bed"),
        purge = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.fa"),
        hap   = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/hap.fa"),
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

        # cutoff file
        cat <<EOF > cutoffs
low=5
mid=35
high=90
EOF

        echo ">>> Splitting LJA contigs" >> {log.out}
        split_fa {input.lja} > {output.split} 2>> {log.err}

        echo ">>> Self alignment (LJA)" >> {log.out}
        minimap2 -xasm5 -DP -t {threads} {output.split} {output.split} \
            | gzip -c > {output.paf} 2>> {log.err}

        echo ">>> Running purge_dups (LJA)" >> {log.out}
        purge_dups -2 -T cutoffs -c PB.base.cov {output.paf} > {output.bed} 2>> {log.err}

echo ">>> Extracting purged sequences (LJA)" >> {log.out}

if [[ ! -s {output.bed} ]]; then
    echo "[INFO] No duplications detected — copying original assembly." >> {log.out}

    # purged = full assembly
    cp {output.split} {output.purge}

    # hap = minimal valid fasta (avoid samtools faidx crash)
    echo -e ">hap_empty\nN" > {output.hap}

    else
    get_seqs {output.bed} {input.lja} 2>> {log.err}

    # If get_seqs failed to produce purged.fa → fall back
    if [[ ! -f purged.fa || ! -s purged.fa ]]; then
        echo "[WARN] get_seqs produced no purged.fa — falling back to full assembly." >> {log.out}

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
