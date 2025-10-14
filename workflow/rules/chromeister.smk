
import pandas as pd

samples = pd.read_csv(config["sample_sheet"], sep="\t", index_col="sample").to_dict(orient="index")

rule chromeister_hifiasm:
    input:
        query = "Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta",
        ref = lambda wc: samples[wc.sample]["Ref"],
        hifiasm_done = "Assemblies/{sample}/HIFIASM/hifiasm.done"
    output:
        mat   = "Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.mat",
        score = "Assemblies/{sample}/CHROMEISTER/{sample}_score.txt",
        done  = "Assemblies/{sample}/CHROMEISTER/.hifiasm.done"
    log:
        out = "Assemblies/{sample}/CHROMEISTER/{sample}.log",
        err = "Assemblies/{sample}/CHROMEISTER/{sample}.err"
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

rule chromeister_lja:
    input:
        query = "Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta",
        ref = lambda wc: samples[wc.sample]["Ref"],
        lja_done = "Assemblies/{sample}/LJA/lja.done"
    output:
        mat   = "Assemblies/{sample}/CHROMEISTER/{sample}-lja.mat",
        score = "Assemblies/{sample}/CHROMEISTER/{sample}_lja_score.txt",
        done  = "Assemblies/{sample}/CHROMEISTER/.lja.done"
    log:
        out = "Assemblies/{sample}/CHROMEISTER/{sample}_lja.log",
        err = "Assemblies/{sample}/CHROMEISTER/{sample}_lja.err"
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

rule chromeister_compare:
    input:
        hifiasm = "Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta",
        lja = "Assemblies/{sample}/FILTERED/{sample}-HIFI-lja.fasta",
        hifiasm_done = "Assemblies/{sample}/HIFIASM/hifiasm.done",
        lja_done = "Assemblies/{sample}/LJA/lja.done"
    output:
        mat   = "Assemblies/{sample}/CHROMEISTER_COMPARE/{sample}-compare.mat",
        score = "Assemblies/{sample}/CHROMEISTER_COMPARE/{sample}_compare_score.txt",
        done  = "Assemblies/{sample}/CHROMEISTER_COMPARE/.compare.done"
    log:
        out = "Assemblies/{sample}/CHROMEISTER_COMPARE/{sample}_compare.log",
        err = "Assemblies/{sample}/CHROMEISTER_COMPARE/{sample}_compare.err"
    conda:
        "../envs/chromeister_env.yaml"
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"

        QUERY_ABS=$(realpath {input.hifiasm}) || {{ echo "Missing {input.hifiasm}" >&2; exit 1; }}
        REF_ABS=$(realpath {input.lja}) || {{ echo "Missing {input.lja}" >&2; exit 1; }}

        cd "$OUTDIR"

        CHROMEISTER -query "$REF_ABS" -db "$QUERY_ABS" -out "$(basename {output.mat})" -dimension 2000 > "$(basename {log.out})" 2>> "$(basename {log.err})"
        Rscript $(which compute_score.R) "$(basename {output.mat})" 2000 > "$(basename {output.score})" 2>> "$(basename {log.err})"

        touch "$(basename {output.done})"
        """
