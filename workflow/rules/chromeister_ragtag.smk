import pandas as pd

# Load sample sheet
samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample")

# Samples that require RagTag
RAGTAG_SAMPLES = samples_df.query("ragtag == 'Y'").index.tolist()


##############################################
#   Chromeister on Purged RagTag Assembly
##############################################
rule chromeister_ragtag_purged:
    input:
        asm = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.reforder.fasta"),
        ref = lambda wc: samples_df.loc[wc.sample, "Ref"],
        ragtag_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.mat"),
        png   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.mat.filt.png"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.purged_ragtag.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    wildcard_constraints:
        sample="|".join(RAGTAG_SAMPLES)
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        # Resolve scaffold path
        ASM_ABS=$(realpath {input.asm}) || {{ echo "Missing assembly: {input.asm}" >&2; exit 1; }}

        # Prepare reference
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            BASENAME=$(basename "$REF_FILE" .gz)
            gunzip -c "$REF_FILE" > "$OUTDIR/$BASENAME"
            REF_FILE="$OUTDIR/$BASENAME"
        fi

        # Run Chromeister
        CHROMEISTER -query "$REF_FILE" \
                    -db "$ASM_ABS" \
                    -out "{output.mat}" \
                    -dimension 2000 \
                    > "{log.out}" 2>> "{log.err}"

        # Compute score
        Rscript $(which compute_score.R) "{output.mat}" 2000 \
            > "{output.score}" 2>> "{log.err}"

        # PNG rendering by Chromeister
        if [[ -f "{output.mat}.png" ]]; then
            mv "{output.mat}.png" "{output.png}"
        fi

        # Remove decompressed reference (temporary)
        if [[ "$REF_FILE" == "$OUTDIR/"* && "$REF_FILE" != *.fa && "$REF_FILE" != *.fna && "$REF_FILE" != *.fasta ]]; then
            rm -f "$REF_FILE"
        fi

        touch "{output.done}"
        """
