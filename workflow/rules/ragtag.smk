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
        purge = outpath("Assemblies/{sample}/PURGE_DUPS/purged.fa")
    output:
        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
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

        # RagTag writes ragtag.scaffold.fasta in OUTDIR
        samtools faidx {output.fasta}

        touch {output.done}
        """


# -----------------------------------------------------
# sort RagTag scaffolds by reference order
# -----------------------------------------------------
rule ragtag_sort_by_ref:
    input:
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        agp   = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.agp"),
        fasta = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta")
    output:
        sorted = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.reforder.fasta")
    threads: 2
    conda: "../envs/ragtag.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname "{output.sorted}")
        mkdir -p "$OUTDIR"

        REF="{input.ref}"
        AGP="{input.agp}"
        FASTA="{input.fasta}"

        REFORDER="$OUTDIR/ref.order.txt"
        SCAFFORDER="$OUTDIR/scaffold.order.txt"

        # ----------------------------------------------------
        # 1) Extract reference order from FASTA headers
        #    (keep order, take first token after '>')
        # ----------------------------------------------------
        if [[ "$REF" == *.gz ]]; then
            zgrep '^>' "$REF" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        else
            grep  '^>' "$REF" | cut -d' ' -f1 | sed 's/^>//' > "$REFORDER"
        fi

        # ----------------------------------------------------
        # 2) Build scaffold order:
        #    - AGP col1 = <ref>_RagTag
        #    - strip _RagTag, match to REFORDER
        #    - output col1 (index) and full scaffold name
        # ----------------------------------------------------
        awk '
            NR==FNR {{ ord[$1]=NR; next }}
            /^#/  {{ next }}
            $1 != prev {{
                name=$1
                sub(/_RagTag$/, "", name)
                if (name in ord)
                    print ord[name] "\t" $1
                prev=$1
            }}
        ' "$REFORDER" "$AGP" \
        | sort -k1,1n \
        | cut -f2 \
        > "$SCAFFORDER"

        if [[ ! -s "$SCAFFORDER" ]]; then
            echo "[ERROR] scaffold.order.txt is empty for sample {wildcards.sample}" >&2
            echo "[HINT] Check that reference IDs in $REF match AGP col1 (without _RagTag)" >&2
            exit 1
        fi

        # ----------------------------------------------------
        # 3) Reorder RagTag scaffold FASTA
        # ----------------------------------------------------
        samtools faidx "$FASTA"

        samtools faidx "$FASTA" $(cat "$SCAFFORDER") \
            > "{output.sorted}"
        """
