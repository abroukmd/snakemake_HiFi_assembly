rule filter_hifiasm:
    input:
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        hifiasm = outpath("Assemblies/{sample}/HIFIASM/{sample}.fasta")
    output:
        filtered = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta"),
        paf = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta.paf"),
        list = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta.list")
    log:
        out = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.log"),
        err = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.err")
    conda:
        "../envs/hifiasm_env.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.filtered})

        minimap2 -x asm5 -t 24 {input.mito} {input.hifiasm} > {output.paf} 2>> {log.out}

        awk '{{P=($11/$2)*100; if (P >= 50) print P"\\t" $0}}' {output.paf} | cut -f 2 | sort | uniq > {output.list}

        if [ -s {output.list} ]; then
            filterbyname.sh in={input.hifiasm} out={output.filtered} names={output.list} overwrite=true exclude 2>> {log.out}
        else
            cp {input.hifiasm} {output.filtered}
        fi

        tmpout=$(mktemp --suffix=.fasta)
        sortbyname.sh in={output.filtered} out=$tmpout overwrite=true length descending 2>> {log.out}

        # Concatenate sorted and mito into a new temp final output
        final_tmp=$(mktemp --suffix=.fasta)
        cat $tmpout {input.mito} > $final_tmp

        # Move to final output
        mv $final_tmp {output.filtered}
        rm $tmpout

        if ! grep -q '^>' {output.filtered}; then
            echo "ERROR: Output does not contain valid FASTA headers (>)" >&2
            exit 1
        fi
        """

rule filter_lja:
    input:
        mito = outpath("Assemblies/{sample}/OATK/{sample}.mito.ctg.fasta"),
        asm = outpath("Assemblies/{sample}/LJA/assembly.fasta")
    output:
        filtered = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta"),
        paf = temp(outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta.paf")),
        list = temp(outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta.list"))
    log:
        stdout = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.log"),
        stderr = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.err")
    conda:
        "../envs/hifiasm_env.yaml"
    threads: 24
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.filtered})

        minimap2 -x asm5 -t {threads} {input.mito} {input.asm} > {output.paf} 2>> {log.stdout}

        awk '{{P=($11/$2)*100; if (P >= 50) print P"\\t" $0}}' {output.paf} | cut -f 2 | sort | uniq > {output.list}

        if [ -s {output.list} ]; then
            filterbyname.sh in={input.asm} out={output.filtered} names={output.list} overwrite=true exclude 2>> {log.stdout}
        else
            cp {input.asm} {output.filtered}
        fi

        tmpout=$(mktemp --suffix=.fasta)
        sortbyname.sh in={output.filtered} out=$tmpout overwrite=true length descending 2>> {log.stdout}

        # Concatenate sorted and mito into a new temp final output
        final_tmp=$(mktemp --suffix=.fasta)
        cat $tmpout {input.mito} > $final_tmp

        # Move to final output
        mv $final_tmp {output.filtered}
        rm $tmpout

        if ! grep -q '^>' {output.filtered}; then
            echo "ERROR: Output does not contain valid FASTA headers (>)" >&2
            exit 1
        fi
        """
