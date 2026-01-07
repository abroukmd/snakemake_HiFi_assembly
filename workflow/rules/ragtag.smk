import pandas as pd

# Load sample sheet
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample")

# Samples requiring RagTag
RAGTAG_SAMPLES = samples_df.query("ragtag == 'Y'").index.tolist()

# -----------------------------------------------------
# RagTag on purge_dups "purged.fa"
# -----------------------------------------------------
rule ragtag_purged:
    input:
        ref   = lambda wc: samples[wc.sample]["Ref"],
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.hifiasm.fa")
    output:
        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
        agp   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.agp"),
        fai   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta.fai"),
        done  = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    threads: 32
    conda: "../envs/ragtag.yaml"
    log:
        out = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.log"),
        err = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.err")
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.fasta})
        mkdir -p "$OUTDIR"

        ragtag.py scaffold \
            -t {threads} \
            -o "$OUTDIR" \
            {input.ref} {input.purge} \
            >> {log.out} 2>> {log.err}

        # RagTag writes ragtag.scaffold.fasta and ragtag.scaffold.agp in OUTDIR
        samtools faidx {output.fasta}

        touch {output.done}
        """

#rule ragtag_purged:
#    input:
#        ref   = lambda wc: samples[wc.sample]["Ref"],
#        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.fa")
#    output:
#        output:
#        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
#        agp   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.agp"),
#        fai   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta.fai"),
#        done  = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
        #fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
        #fai   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta.fai"),
        #done  = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
#    threads: 32
#    conda: "../envs/ragtag.yaml"
#    log:
#        out = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.log"),
#        err = outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.err")
#    shell:
#        r"""
#        set -euo pipefail
#
#        OUTDIR=$(dirname {output.fasta})
#        mkdir -p "$OUTDIR"
#
#        ragtag.py scaffold \
#            -t {threads} \
#            -o "$OUTDIR" \
#            {input.ref} {input.purge} \
#            >> {log.out} 2>> {log.err}
#
#        # RagTag writes ragtag.scaffold.fasta in OUTDIR
#        samtools faidx {output.fasta}
#
#        touch {output.done}
#        """


# -----------------------------------------------------
# sort RagTag scaffolds by reference order
# -----------------------------------------------------
rule ragtag_sort_by_ref:
    input:
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        agp   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.agp"),
        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta")
    output:
        sorted = outpath(
            "Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/"
            "ragtag.scaffold.hifiasm.reforder.fasta"
        )
    threads: 2
    conda: "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname "{output.sorted}")
        mkdir -p "$OUTDIR"

        REFORDER="$OUTDIR/ref.order.txt"
        SCAFFORDER="$OUTDIR/scaffold.order.txt"

        if [[ "{input.ref}" == *.gz ]]; then
            zgrep '^>' "{input.ref}" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        else
            grep '^>' "{input.ref}" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        fi

        awk '
            NR==FNR { ord[$1]=NR; next }
            /^#/ { next }
            $1 != prev {
                name=$1
                sub(/_RagTag$/, "", name)
                if (name in ord)
                    print ord[name] "\t" $1
                prev=$1
            }
        ' "$REFORDER" "{input.agp}" \
        | sort -k1,1n | cut -f2 > "$SCAFFORDER"

        samtools faidx "{input.fasta}"
        samtools faidx "{input.fasta}" $(cat "$SCAFFORDER") > "{output.sorted}"
        """



################ LJA #######################
#rule ragtag_purged_lja:
#    input:
#        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
#        purge = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.fa")
#    output:
#        fasta = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.fasta"),
#        agp   = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.agp"),
#        done  = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}_ragtag_lja.done")
#    conda: "../envs/ragtag.yaml"
#    threads: 32
#    log:
#        out = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.log"),
#        err = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.err")
#    shell:
#        r"""
#        set -euo pipefail
#        OUTDIR=$(dirname {output.fasta})
#        mkdir -p "$OUTDIR"
#
#        ragtag.py scaffold \
#            -t {threads} \
#            -o "$OUTDIR" \
#            {input.ref} {input.purge} \
#            >> {log.out} 2>> {log.err}
#
#        # samtools index
#        samtools faidx {output.fasta}
#
#        touch {output.done}
#        """

rule ragtag_purged_lja:
    input:
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        purge = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.lja.fa")
    output:
        fasta = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.fasta"),
        agp   = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.agp"),
        done  = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}_ragtag_lja.done")
    threads: 32
    conda: "../envs/ragtag.yaml"
    log:
        out = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.log"),
        err = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.err")
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.fasta})
        mkdir -p "$OUTDIR"

        # IMPORTANT: remove any stale RagTag outputs so we don't reuse old .paf/.agp
        rm -f "$OUTDIR"/ragtag.scaffold.*

        echo "[INFO] Running RagTag LJA for {wildcards.sample}" >> {log.out}

        ragtag.py scaffold \
            -t {threads} \
            -o "$OUTDIR" \
            {input.ref} {input.purge} \
            >> {log.out} 2>> {log.err}

        # Sanity check: RagTag should have produced the scaffold FASTA
        if [[ ! -s {output.fasta} ]]; then
            echo "[ERROR] RagTag did not produce a non-empty scaffold FASTA for {wildcards.sample}" >&2
            exit 1
        fi

        # Index scaffold FASTA
        samtools faidx {output.fasta}

        touch {output.done}
        """


rule ragtag_sort_by_ref_lja:
    input:
        ref   = lambda wc: sample_refs[wc.sample],
        agp   = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.agp"),
        fasta = outpath("Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/ragtag.scaffold.fasta")
    output:
        sorted = outpath(
            "Assemblies/{sample}/RAGTAG_LJA/purged-lja_ragtag/"
            "ragtag.scaffold.lja.reforder.fasta"
        )
    conda: "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname "{output.sorted}")
        mkdir -p "$OUTDIR"

        REFORDER="$OUTDIR/ref.order.txt"
        SCAFFORDER="$OUTDIR/scaffold.order.txt"

        if [[ "{input.ref}" == *.gz ]]; then
            zgrep '^>' "{input.ref}" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        else
            grep '^>' "{input.ref}" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        fi

        awk '
            NR==FNR { ord[$1]=NR; next }
            /^#/ { next }
            $1 != prev {
                name=$1
                sub(/_RagTag$/, "", name)
                if (name in ord)
                    print ord[name] "\t" $1
                prev=$1
            }
        ' "$REFORDER" "{input.agp}" \
        | sort -k1,1n | cut -f2 > "$SCAFFORDER"

        samtools faidx "{input.fasta}"
        samtools faidx "{input.fasta}" $(cat "$SCAFFORDER") > "{output.sorted}"
        """
