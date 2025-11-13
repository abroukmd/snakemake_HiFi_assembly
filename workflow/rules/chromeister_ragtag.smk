import pandas as pd

samples = pd.read_csv(config["sample_sheet"], sep="\t", index_col="sample").to_dict(orient="index")


##############################################
#   Chromeister on Purged RagTag Assembly
##############################################
rule chromeister_ragtag_purged:
    input:
        asm = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/ragtag.scaffold.fasta"),
        ref = lambda wc: samples[wc.sample]["Ref"],
        ragtag_done = outpath("Assemblies/{sample}/RAGTAG/purged-hifiasm_ragtag/.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.mat"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.purged_ragtag.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-purged_ragtag.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail

        # Create output directory
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        # Resolve absolute path of scaffold
        ASM_ABS=$(realpath {input.asm}) || {{ echo "Missing assembly {input.asm}" >&2; exit 1; }}

        # Handle gzipped reference
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            BASENAME=$(basename "$REF_FILE" .gz)
            gunzip -c "$REF_FILE" > "$OUTDIR/$BASENAME"
            REF_FILE="$OUTDIR/$BASENAME"
        fi

        # Run Chromeister
        CHROMEISTER -query "$REF_FILE" -db "$ASM_ABS" \
            -out "{output.mat}" -dimension 2000 \
            > "{log.out}" 2>> "{log.err}"

        # Compute score
        Rscript $(which compute_score.R) "{output.mat}" 2000 \
            > "{output.score}" 2>> "{log.err}"

        # Remove decompressed temporary reference if it was created
        if [[ "$REF_FILE" == "$OUTDIR/"* && "$REF_FILE" != *.fa && "$REF_FILE" != *.fna && "$REF_FILE" != *.fasta ]]; then
            rm -f "$REF_FILE"
        fi

        touch "{output.done}"
        """







#import pandas as pd
#
## Load the same sample sheet and use the same reference mapping logic
#samples = pd.read_csv(config["sample_sheet"], sep="\t", index_col="sample").to_dict(orient="index")
#
#rule chromeister_ragtag_hifiasm:
#    input:
#        query = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.fasta"),
#        ref = lambda wc: samples[wc.sample]["Ref"],
#        ragtag_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.done")
#    output:
#        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-hifiasm.mat"),
#        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}_score.txt"),
#        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.hifiasm.done")
#    log:
#        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}.log"),
#        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}.err")
#    conda:
#        "../envs/chromeister_env.yaml"
#    shell:
#        r"""
#        set -euo pipefail
#
#        OUTDIR=$(dirname {output.mat})
#        mkdir -p "$OUTDIR"
#
#        QUERY_ABS=$(realpath {input.query}) || {{ echo "Missing {input.query}" >&2; exit 1; }}
#        REF_FILE="{input.ref}"
#
#        # Unzip reference if compressed
#        if [[ "$REF_FILE" == *.gz ]]; then
#            BASENAME=$(basename "$REF_FILE" .gz)
#            gunzip -c "$REF_FILE" > "$OUTDIR/$BASENAME"
#            REF_FILE="$OUTDIR/$BASENAME"
#        fi
#
#        CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2>> "{log.err}"
#        Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"
#
#        # Cleanup temporary unzipped ref
#        if [[ "$REF_FILE" == "$OUTDIR/"* && "$REF_FILE" != *.fa && "$REF_FILE" != *.fasta ]]; then
#            rm -f "$REF_FILE"
#        fi
#
#        touch "{output.done}"
#        """

