
rule hifiasm:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        fasta = "Assemblies/{sample}/HIFIASM/{sample}.fasta",
        done = "Assemblies/{sample}/HIFIASM/hifiasm.done"
    conda:
        "../envs/hifiasm_env.yaml"
    params:
        cores = 32
    log:
        "Assemblies/{sample}/HIFIASM/hifiasm.log"
    shell:
        """
        set -euo pipefail

        mkdir -p Assemblies/{wildcards.sample}/HIFIASM
        hifiasm -o Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample} -t {params.cores} {input.fastq} > {log} 2>&1

        # Convert GFA to FASTA
        awk '/^S/{{print ">"$2"\\n"$3}}' \
            Assemblies/{wildcards.sample}/HIFIASM/{wildcards.sample}.bp.p_ctg.gfa > {output.fasta}

        touch {output.done}
        """


rule lja:
    input:
        fastq = lambda wildcards: sample_paths[wildcards.sample]
    output:
        fasta = "Assemblies/{sample}/LJA/assembly.fasta",
        done = "Assemblies/{sample}/LJA/lja.done"
    conda:
        "../envs/lja_env.yaml"
    threads: 32
    log:
        "Assemblies/{sample}/LJA/lja.log"
    shell:
        """
        set -euo pipefail

        lja --reads {input.fastq} -o Assemblies/{wildcards.sample}/LJA --threads {threads} > {log} 2>&1

        # Debug: list output files
        ls -lh Assemblies/{wildcards.sample}/LJA >> {log}

        # Check and move the expected output, fail gracefully if not found
        if [[ ! -f Assemblies/{wildcards.sample}/LJA/assembly.fasta ]]; then
            echo "[ERROR] LJA did not produce assembly.fasta" >&2
            exit 1
        fi

#        mv Assemblies/{wildcards.sample}/LJA/assembly.fasta {output.fasta}
        touch {output.done}
        """
