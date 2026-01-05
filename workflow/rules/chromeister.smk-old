import pandas as pd

samples = pd.read_csv(config["sample_sheet"], sep="\t", index_col="sample").to_dict(orient="index")

rule chromeister_hifiasm:
    input:
        query = outpath("Assemblies/{sample}/FILTERED/{sample}-HIFI-hifiasm.fasta"),
        ref = lambda wc: samples[wc.sample]["Ref"],
        hifiasm_done = outpath("Assemblies/{sample}/HIFIASM/hifiasm.done")
    output:
        mat   = outpath("Assemblies/{sample}/CHROMEISTER/{sample}-hifiasm.mat"),
        score = outpath("Assemblies/{sample}/CHROMEISTER/{sample}_score.txt"),
        done  = outpath("Assemblies/{sample}/CHROMEISTER/.hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/CHROMEISTER/{sample}.log"),
        err = outpath("Assemblies/{sample}/CHROMEISTER/{sample}.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
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

