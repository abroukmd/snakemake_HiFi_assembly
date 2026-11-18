#############################################
# HiFi and LJA assembly rules
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

        echo "[INFO] Running Hifiasm for {wildcards.sample}" > {log:q}
        hifiasm -o Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample} \
                 -t {threads} {input.fastq:q} >> {log:q} 2>&1

        # Convert GFA to FASTA
        awk '/^S/{{print ">"$2"\n"$3}}' \
            Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample}.bp.p_ctg.gfa \
            > {output.fasta:q}

        test -s {output.fasta:q}
        grep -q '^>' {output.fasta:q}
        touch {output.done:q}
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

        echo "[INFO] Running LJA for {wildcards.sample}" > {log:q}
        mkdir -p {params.out_dir:q}

        # Run LJA with full absolute output path
        lja --reads {input.fastq:q} \
            -o {params.out_dir:q} \
            --threads {threads} >> {log:q} 2>&1

        echo "[INFO] Listing LJA output files" >> {log:q}
        ls -lh {params.out_dir:q} >> {log:q}

        # Check expected output
        if [[ ! -s {output.fasta:q} ]] || ! grep -q '^>' {output.fasta:q}; then
            echo "[ERROR] LJA did not produce a nonempty assembly FASTA" >&2
            echo "[ERROR] LJA did not produce a nonempty assembly FASTA" >> {log:q}
            exit 1
        fi

        # Mark completion
        touch {output.done:q}
        """

