import pandas as pd

# Load the same sample sheet and use the same reference mapping logic
samples = pd.read_csv(config["sample_sheet"], sep="\t", index_col="sample").to_dict(orient="index")

rule chromeister_ragtag_hifiasm:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.fasta"),
        ref = lambda wc: samples[wc.sample]["Ref"],
        ragtag_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_hifiasm.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-hifiasm.mat"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}.err")
    conda:
        "../envs/chromeister_env.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        QUERY_ABS=$(realpath {input.query}) || {{ echo "Missing {input.query}" >&2; exit 1; }}
        REF_FILE="{input.ref}"

        # Unzip reference if compressed
        if [[ "$REF_FILE" == *.gz ]]; then
            BASENAME=$(basename "$REF_FILE" .gz)
            gunzip -c "$REF_FILE" > "$OUTDIR/$BASENAME"
            REF_FILE="$OUTDIR/$BASENAME"
        fi

        CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2>> "{log.err}"
        Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"

        # Cleanup temporary unzipped ref
        if [[ "$REF_FILE" == "$OUTDIR/"* && "$REF_FILE" != *.fa && "$REF_FILE" != *.fasta ]]; then
            rm -f "$REF_FILE"
        fi

        touch "{output.done}"
        """

rule chromeister_ragtag_lja:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_lja.fasta"),
        ref = lambda wc: samples[wc.sample]["Ref"],
        ragtag_done = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_lja.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-lja.mat"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}_lja_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.lja.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}_lja.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}_lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        QUERY_ABS=$(realpath {input.query}) || {{ echo "Missing {input.query}" >&2; exit 1; }}
        REF_FILE="{input.ref}"

        if [[ "$REF_FILE" == *.gz ]]; then
            BASENAME=$(basename "$REF_FILE" .gz)
            gunzip -c "$REF_FILE" > "$OUTDIR/$BASENAME"
            REF_FILE="$OUTDIR/$BASENAME"
        fi

        CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2>> "{log.err}"
        Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"

        if [[ "$REF_FILE" == "$OUTDIR/"* && "$REF_FILE" != *.fa && "$REF_FILE" != *.fasta ]]; then
            rm -f "$REF_FILE"
        fi

        touch "{output.done}"
        """

rule chromeister_ragtag_compare:
    input:
        hifiasm_mat = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-hifiasm.mat"),
        lja_mat = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/{sample}-lja.mat"),
        hifiasm_done = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.hifiasm.done"),
        lja_done = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG/.lja.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG_COMPARE/{sample}-compare.mat"),
        score = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG_COMPARE/{sample}_compare_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG_COMPARE/.compare.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG_COMPARE/{sample}_compare.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER_RAGTAG_COMPARE/{sample}_compare.err")
    conda:
        "../envs/chromeister_env.yaml"
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        QUERY_ABS=$(realpath {input.hifiasm_mat}) || {{ echo "Missing {input.hifiasm_mat}" >&2; exit 1; }}
        REF_ABS=$(realpath {input.lja_mat}) || {{ echo "Missing {input.lja_mat}" >&2; exit 1; }}

        cd "$OUTDIR"

        CHROMEISTER -query "$REF_ABS" -db "$QUERY_ABS" -out "$(basename {output.mat})" -dimension 2000 > "$(basename {log.out})" 2>> "$(basename {log.err})"
        Rscript $(which compute_score.R) "$(basename {output.mat})" 2000 > "$(basename {output.score})" 2>> "$(basename {log.err})"

        touch "$(basename {output.done})"
        """
