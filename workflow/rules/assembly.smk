#############################################
# HiFi assembly rules
#############################################

# ---------- HIFIASM ----------
rule hifiasm:
    input:
        fastq = lambda wc: sample_paths[wc.sample]
    output:
        fasta = outpath("Assemblies/{sample}/HIFIASM/{sample}.fasta"),
        done  = outpath("Assemblies/{sample}/HIFIASM/hifiasm.done")
    conda:
        "../envs/hifiasm_env.yaml"
    threads: 32
    log:
        outpath("Assemblies/{sample}/HIFIASM/hifiasm.log")
    shell:
        r"""
        set -euo pipefail

        mkdir -p Assemblies/{wildcards.sample}/HIFIASM

        echo "[INFO] Running Hifiasm for {wildcards.sample}" > {log}
        hifiasm -o Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample} \
                -t {threads} {input.fastq} >> {log} 2>&1

        # Convert GFA to FASTA
        awk '/^S/{{print ">"$2"\n"$3}}' \
            Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample}.bp.p_ctg.gfa \
            > {output.fasta}

        touch {output.done}
        """


# ---------- LJA ----------
rule lja:
    input:
        fastq = lambda wc: sample_paths[wc.sample]
    output:
        fasta = outpath("Assemblies/{sample}/LJA/assembly.fasta"),
        done  = outpath("Assemblies/{sample}/LJA/lja.done")
    conda:
        "../envs/lja_env.yaml"
    threads: 32
    params:
        out_dir = lambda wc: outpath(f"Assemblies/{wc.sample}/LJA")
    log:
        outpath("Assemblies/{sample}/LJA/lja.log")
    shell:
        r"""
        set -euo pipefail

        echo "[INFO] Running LJA for {wildcards.sample}" > {log}
        mkdir -p {params.out_dir}

        # Run LJA with full absolute output path
        lja --reads {input.fastq} \
            -o {params.out_dir} \
            --threads {threads} >> {log} 2>&1

        echo "[INFO] Listing contents of {params.out_dir}" >> {log}
        ls -lh {params.out_dir} >> {log}

        # Check expected output
        if [[ ! -f {params.out_dir}/assembly.fasta ]]; then
            echo "[ERROR] LJA did not produce assembly.fasta" >&2
            echo "[ERROR] LJA did not produce assembly.fasta" >> {log}
            exit 1
        fi

        # Mark completion
        touch {output.done}
        """
