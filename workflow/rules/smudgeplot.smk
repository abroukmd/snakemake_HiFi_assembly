import os

def normalize_sample(wc):
    return wc.sample.split("/")[-1]

def get_fq1_abs(wc):
    sample = normalize_sample(wc)
    return sample_fq1[sample]

def get_fq2_abs(wc):
    sample = normalize_sample(wc)
    return sample_fq2[sample]

rule smudgeplot:
    input:
        fq1 = get_fq1_abs,
        fq2 = get_fq2_abs
    log:
        stdout = outpath("Assemblies/{sample}/QC/smudgeplot/smudgeplot.log"),
        stderr = outpath("Assemblies/{sample}/QC/smudgeplot/smudgeplot.err")
    output:
        done = outpath("Assemblies/{sample}/QC/smudgeplot/smudgeplot.done"),
        outdir = directory(outpath("Assemblies/{sample}/QC/smudgeplot"))
    conda:
        "../envs/smudgeplot_env.yaml"
    params:
        cores=51,
        mem="100gb",
        time="01:00:00"
    shell:
        """
    set -euo pipefail

    mkdir -p $(dirname {output.done})
    SAMPLE_NAME=$(basename {wildcards.sample})

    for k in 31 51; do
        #echo "smudeplot - Processing k=$k...";
        mkdir -p {output[1]}/smudge_k$k;
        FastK -v -t{params.cores} -k$k -M64 -T{params.cores} {input.fq1} {input.fq2} -N{output[1]}/smudge_k$k/FastK_Table 2>>{log.stdout};
        smudgeplot.py hetmers -L 5 -t {params.cores} -o {output[1]}/smudge_k$k/kmerpairs --verbose {output[1]}/smudge_k$k/FastK_Table 2>>{log.stdout};
        if [ -f {output[1]}/smudge_k$k/kmerpairs_text.smu ]; then
            smudgeplot.py all -o {output[1]}/smudge_k$k/$SAMPLE_NAME {output[1]}/smudge_k$k/kmerpairs_text.smu 2>>{log.stdout};
        else
            echo "WARNING: kmerpairs_text.smu not found for k=$k" >&2;
        fi;
        #echo "smudeplot - Finished k=$k.";
    done

    touch {output.done}

        """
