include: "../accessors.smk"
import pandas as pd

samples_df = pd.read_csv(config["sample_sheet"], sep="\t").set_index("sample")

rule chromeister_final_hifiasm:
    input:
        query = get_final_hifiasm_abs,
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        done  = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/{sample}.purged.final.fasta")
    output:
        mat   = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.mat"),
        png   = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.mat.png"),
        score = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm_score.txt"),
        done  = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/.final_hifiasm.done")
    log:
        out = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.log"),
        err = outpath("{sample}_FINAL_ASSEMBLIES/hifiasm/CHROMEISTER/{sample}-final_hifiasm.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"
        QUERY_ABS=$(realpath {input.query})
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            REF_LOCAL="$OUTDIR/$(basename "$REF_FILE" .gz)"
            gunzip -c "$REF_FILE" > "$REF_LOCAL"
            REF_FILE="$REF_LOCAL"
        fi
        CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2> "{log.err}"
        Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"
        if [[ -f "{output.mat}.png" ]]; then
            mv "{output.mat}.png" "{output.png}"
        else
            touch "{output.png}"
        fi
        touch "{output.done}"
        """

rule chromeister_final_lja:
    input:
        query = get_final_lja_abs,
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        done  = outpath("{sample}_FINAL_ASSEMBLIES/lja/{sample}.purged.final.fasta")
    output:
        mat   = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.mat"),
        png   = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.mat.png"),
        score = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja_score.txt"),
        done  = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/.final_lja.done")
    log:
        out = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.log"),
        err = outpath("{sample}_FINAL_ASSEMBLIES/lja/CHROMEISTER/{sample}-final_lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"
        QUERY_ABS=$(realpath {input.query})
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            REF_LOCAL="$OUTDIR/$(basename "$REF_FILE" .gz)"
            gunzip -c "$REF_FILE" > "$REF_LOCAL"
            REF_FILE="$REF_LOCAL"
        fi
        CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2> "{log.err}"
        Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"
        if [[ -f "{output.mat}.png" ]]; then
            mv "{output.mat}.png" "{output.png}"
        else
            touch "{output.png}"
        fi
        touch "{output.done}"
        """

rule chromeister_final_hap_hifiasm:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/{sample}.hap.final.fasta"),
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        done  = outpath("Assemblies/{sample}/RAGTAG/{sample}_ragtag_purged.done")
    output:
        mat   = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.mat"),
        png   = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.mat.png"),
        score = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm_score.txt"),
        done  = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/.final_hap_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.log"),
        err = outpath("Assemblies/{sample}/RAGTAG/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_hifiasm.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"
        QUERY_ABS=$(realpath {input.query})
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            REF_LOCAL="$OUTDIR/$(basename "$REF_FILE" .gz)"
            gunzip -c "$REF_FILE" > "$REF_LOCAL"
            REF_FILE="$REF_LOCAL"
        fi
        if grep -q '^>hap_empty' "$QUERY_ABS"; then
            : > "{output.mat}"
            : > "{output.png}"
            echo "SKIPPED hap_empty" > "{output.score}"
        else
            CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2> "{log.err}"
            Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"
            if [[ -f "{output.mat}.png" ]]; then
                mv "{output.mat}.png" "{output.png}"
            else
                touch "{output.png}"
            fi
        fi
        touch "{output.done}"
        """

rule chromeister_final_hap_lja:
    input:
        query = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/{sample}.hap.final.fasta"),
        ref   = lambda wc: samples_df.loc[wc.sample, "Ref"],
        done  = outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}_ragtag_lja.done")
    output:
        mat   = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.mat"),
        png   = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.mat.png"),
        score = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja_score.txt"),
        done  = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/.final_hap_lja.done")
    log:
        out = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.log"),
        err = outpath("Assemblies/{sample}/RAGTAG_LJA/FINAL_{sample}/CHROMEISTER/{sample}-final_hap_lja.err")
    conda:
        "../envs/chromeister_env.yaml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        OUTDIR=$(dirname {output.mat})
        mkdir -p "$OUTDIR"
        QUERY_ABS=$(realpath {input.query})
        REF_FILE="{input.ref}"
        if [[ "$REF_FILE" == *.gz ]]; then
            REF_LOCAL="$OUTDIR/$(basename "$REF_FILE" .gz)"
            gunzip -c "$REF_FILE" > "$REF_LOCAL"
            REF_FILE="$REF_LOCAL"
        fi
        if grep -q '^>hap_empty' "$QUERY_ABS"; then
            : > "{output.mat}"
            : > "{output.png}"
            echo "SKIPPED hap_empty" > "{output.score}"
        else
            CHROMEISTER -query "$REF_FILE" -db "$QUERY_ABS" -out "{output.mat}" -dimension 2000 > "{log.out}" 2> "{log.err}"
            Rscript $(which compute_score.R) "{output.mat}" 2000 > "{output.score}" 2>> "{log.err}"
            if [[ -f "{output.mat}.png" ]]; then
                mv "{output.mat}.png" "{output.png}"
            else
                touch "{output.png}"
            fi
        fi
        touch "{output.done}"
        """
