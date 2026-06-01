PURGE_LOW = int(config.get("purge_dups_low", 5))
PURGE_MID = int(config.get("purge_dups_mid", 35))
PURGE_HIGH = int(config.get("purge_dups_high", 90))


def _purge_len_sorted_path(sample, asm):
    if asm == "hifiasm":
        return outpath(f"Assemblies/{sample}/PURGE_DUPS/purged-hifiasm.sorted.fa")
    if asm == "lja":
        return outpath(f"Assemblies/{sample}/PURGE_DUPS_LJA/purged-lja.sorted.fa")
    raise ValueError(f"Unsupported asm: {asm}")


rule purge_dups:
    input:
        asm = outpath("Assemblies/{sample}/HIFIASM/{sample}.fasta")
    output:
        split      = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.fasta"),
        paf        = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.split.self.paf.gz"),
        bed        = outpath("Assemblies/{sample}/PURGE_DUPS/dups_{sample}.bed"),
        clean_bed  = outpath("Assemblies/{sample}/PURGE_DUPS/dups_{sample}.clean.bed"),
        oneline    = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}.asm.clean.oneline.fa"),
        purge      = outpath("Assemblies/{sample}/PURGE_DUPS/purged.hifiasm.fa"),
        hap        = outpath("Assemblies/{sample}/PURGE_DUPS/hap.hifiasm.fa"),
        purge_fai  = outpath("Assemblies/{sample}/PURGE_DUPS/purged.hifiasm.fa.fai"),
        hap_fai    = outpath("Assemblies/{sample}/PURGE_DUPS/hap.hifiasm.fa.fai"),
        sorted     = outpath("Assemblies/{sample}/PURGE_DUPS/purged-hifiasm.sorted.fa"),
        sorted_fai = outpath("Assemblies/{sample}/PURGE_DUPS/purged-hifiasm.sorted.fa.fai"),
        done       = outpath("Assemblies/{sample}/PURGE_DUPS/purge_hifiasm.done")
    log:
        out = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}-purge_dups.log"),
        err = outpath("Assemblies/{sample}/PURGE_DUPS/{sample}-purge_dups.err")
    conda:
        "../envs/purge_dups.yaml"
    threads: 24
    params:
        low  = PURGE_LOW,
        mid  = PURGE_MID,
        high = PURGE_HIGH
    shell:
        r"""
        set -euo pipefail

        OUTDIR="$(dirname "{output.purge}")"
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"

        fasta_to_oneline() {{
            local in="$1"
            local out="$2"
            awk '
                BEGIN {{ seq=""; hdr="" }}
                /^>/ {{
                    if (hdr != "") {{ print hdr; print seq }}
                    hdr=$0
                    seq=""
                    next
                }}
                {{
                    gsub(/\r$/, "", $0)
                    if ($0 != "") seq = seq $0
                }}
                END {{
                    if (hdr != "") {{ print hdr; print seq }}
                }}
            ' "$in" > "$out"
        }}

        sort_fasta_by_len_desc() {{
            local in="$1"
            local out="$2"
            local tmp_oneline
            tmp_oneline=$(mktemp -p . {wildcards.sample}.purge_dups.oneline.XXXXXX.fa)
            fasta_to_oneline "$in" "$tmp_oneline"
            awk '
                NR % 2 == 1 {{ hdr=$0; next }}
                NR % 2 == 0 {{ print hdr "\t" length($0) "\t" $0 }}
            ' "$tmp_oneline" \
            | sort -k2,2nr \
            | awk '{{ print $1; print $3 }}' > "$out"
            rm -f "$tmp_oneline"
        }}

        cat > cutoffs <<EOF
low={params.low}
mid={params.mid}
high={params.high}
EOF

        echo "[INFO] split_fa on {input.asm}" > "{log.out}"
        split_fa "{input.asm}" > "{output.split}" 2>> "{log.err}"

        nseq=$(grep -c '^>' "{output.split}" || true)
        if [[ "$nseq" -eq 0 ]]; then
            echo "[ERROR] split_fa produced an empty FASTA" >> "{log.out}"
            exit 1
        fi

        echo "[INFO] minimap2 self-alignment" >> "{log.out}"
        minimap2 -xasm5 -DP -t {threads} "{output.split}" "{output.split}" \
            | gzip -c > "{output.paf}" 2>> "{log.err}"

        echo "[INFO] purge_dups -2 -T cutoffs (no PB.base.cov)" >> "{log.out}"
        purge_dups -2 -T cutoffs "{output.paf}" > "{output.bed}" 2>> "{log.err}"

        awk 'BEGIN{{OFS="\t"}} $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3>$2{{print}}' \
            "{output.bed}" > "{output.clean_bed}"

        clean_bed_n=$(wc -l < "{output.clean_bed}" | tr -d ' ')
        echo "[INFO] clean BED lines: $clean_bed_n" >> "{log.out}"

        fasta_to_oneline "{input.asm}" "{output.oneline}"
        samtools faidx "{output.oneline}" >> "{log.out}" 2>> "{log.err}" || true

        if [[ "$clean_bed_n" -eq 0 ]]; then
            echo "[INFO] BED empty after cleaning; copying original assembly" >> "{log.out}"
            cp -f "{input.asm}" "{output.purge}"
            echo -e ">hap_empty\nN" > "{output.hap}"
        else
            echo "[INFO] Running get_seqs -e -c" >> "{log.out}"
            rm -f purged.fa hap.fa
            set +e
            get_seqs -e -c "{output.clean_bed}" "{output.oneline}" >> "{log.out}" 2>> "{log.err}"
            rc=$?
            set -e

            if [[ $rc -ne 0 ]]; then
                echo "[WARN] get_seqs failed (exit $rc); falling back to original assembly" >> "{log.out}"
                cp -f "{input.asm}" "{output.purge}"
                echo -e ">hap_empty\nN" > "{output.hap}"
            elif [[ ! -f purged.fa || ! -s purged.fa ]]; then
                echo "[WARN] purged.fa missing/empty after get_seqs; falling back" >> "{log.out}"
                cp -f "{input.asm}" "{output.purge}"
                echo -e ">hap_empty\nN" > "{output.hap}"
            else
                mv -f purged.fa "{output.purge}"
                if [[ -f hap.fa && -s hap.fa ]]; then
                    mv -f hap.fa "{output.hap}"
                else
                    echo -e ">hap_empty\nN" > "{output.hap}"
                fi
            fi
        fi

        samtools faidx "{output.purge}" >> "{log.out}" 2>> "{log.err}"
        samtools faidx "{output.hap}"   >> "{log.out}" 2>> "{log.err}"

        echo "[INFO] Creating length-sorted purged assembly" >> "{log.out}"
        sort_fasta_by_len_desc "{output.purge}" "{output.sorted}"
        samtools faidx "{output.sorted}" >> "{log.out}" 2>> "{log.err}"

        touch "{output.done}"
        """


rule purge_dups_lja:
    input:
        asm = outpath("Assemblies/{sample}/LJA/assembly.fasta")
    output:
        split      = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.fasta"),
        paf        = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.split.self.paf.gz"),
        bed        = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/dups_{sample}.bed"),
        clean_bed  = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/dups_{sample}.clean.bed"),
        oneline    = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}.asm.clean.oneline.fa"),
        purge      = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.lja.fa"),
        hap        = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/hap.lja.fa"),
        purge_fai  = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged.lja.fa.fai"),
        hap_fai    = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/hap.lja.fa.fai"),
        sorted     = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged-lja.sorted.fa"),
        sorted_fai = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purged-lja.sorted.fa.fai"),
        done       = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/purge_lja.done")
    log:
        out = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}-purge_dups_lja.log"),
        err = outpath("Assemblies/{sample}/PURGE_DUPS_LJA/{sample}-purge_dups_lja.err")
    conda:
        "../envs/purge_dups.yaml"
    threads: 24
    params:
        low  = PURGE_LOW,
        mid  = PURGE_MID,
        high = PURGE_HIGH
    shell:
        r"""
        set -euo pipefail

        OUTDIR="$(dirname "{output.purge}")"
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"

        fasta_to_oneline() {{
            local in="$1"
            local out="$2"
            awk '
                BEGIN {{ seq=""; hdr="" }}
                /^>/ {{
                    if (hdr != "") {{ print hdr; print seq }}
                    hdr=$0
                    seq=""
                    next
                }}
                {{
                    gsub(/\r$/, "", $0)
                    if ($0 != "") seq = seq $0
                }}
                END {{
                    if (hdr != "") {{ print hdr; print seq }}
                }}
            ' "$in" > "$out"
        }}

        sort_fasta_by_len_desc() {{
            local in="$1"
            local out="$2"
            local tmp_oneline
            tmp_oneline=$(mktemp -p . {wildcards.sample}.purge_dups_lja.oneline.XXXXXX.fa)
            fasta_to_oneline "$in" "$tmp_oneline"
            awk '
                NR % 2 == 1 {{ hdr=$0; next }}
                NR % 2 == 0 {{ print hdr "\t" length($0) "\t" $0 }}
            ' "$tmp_oneline" \
            | sort -k2,2nr \
            | awk '{{ print $1; print $3 }}' > "$out"
            rm -f "$tmp_oneline"
        }}

        cat > cutoffs <<EOF
low={params.low}
mid={params.mid}
high={params.high}
EOF

        echo "[INFO] split_fa on {input.asm}" > "{log.out}"
        split_fa "{input.asm}" > "{output.split}" 2>> "{log.err}"

        nseq=$(grep -c '^>' "{output.split}" || true)
        if [[ "$nseq" -eq 0 ]]; then
            echo "[ERROR] split_fa produced an empty FASTA" >> "{log.out}"
            exit 1
        fi

        echo "[INFO] minimap2 self-alignment" >> "{log.out}"
        minimap2 -xasm5 -DP -t {threads} "{output.split}" "{output.split}" \
            | gzip -c > "{output.paf}" 2>> "{log.err}"

        echo "[INFO] purge_dups -2 -T cutoffs (no PB.base.cov)" >> "{log.out}"
        purge_dups -2 -T cutoffs "{output.paf}" > "{output.bed}" 2>> "{log.err}"

        awk 'BEGIN{{OFS="\t"}} $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $3>$2{{print}}' \
            "{output.bed}" > "{output.clean_bed}"

        clean_bed_n=$(wc -l < "{output.clean_bed}" | tr -d ' ')
        echo "[INFO] clean BED lines: $clean_bed_n" >> "{log.out}"

        fasta_to_oneline "{input.asm}" "{output.oneline}"
        samtools faidx "{output.oneline}" >> "{log.out}" 2>> "{log.err}" || true

        if [[ "$clean_bed_n" -eq 0 ]]; then
            echo "[INFO] BED empty after cleaning; copying original assembly" >> "{log.out}"
            cp -f "{input.asm}" "{output.purge}"
            echo -e ">hap_empty\nN" > "{output.hap}"
        else
            echo "[INFO] Running get_seqs -e -c" >> "{log.out}"
            rm -f purged.fa hap.fa
            set +e
            get_seqs -e -c "{output.clean_bed}" "{output.oneline}" >> "{log.out}" 2>> "{log.err}"
            rc=$?
            set -e

            if [[ $rc -ne 0 ]]; then
                echo "[WARN] get_seqs failed (exit $rc); falling back to original assembly" >> "{log.out}"
                cp -f "{input.asm}" "{output.purge}"
                echo -e ">hap_empty\nN" > "{output.hap}"
            elif [[ ! -f purged.fa || ! -s purged.fa ]]; then
                echo "[WARN] purged.fa missing/empty after get_seqs; falling back" >> "{log.out}"
                cp -f "{input.asm}" "{output.purge}"
                echo -e ">hap_empty\nN" > "{output.hap}"
            else
                mv -f purged.fa "{output.purge}"
                if [[ -f hap.fa && -s hap.fa ]]; then
                    mv -f hap.fa "{output.hap}"
                else
                    echo -e ">hap_empty\nN" > "{output.hap}"
                fi
            fi
        fi

        samtools faidx "{output.purge}" >> "{log.out}" 2>> "{log.err}"
        samtools faidx "{output.hap}"   >> "{log.out}" 2>> "{log.err}"

        echo "[INFO] Creating length-sorted purged assembly" >> "{log.out}"
        sort_fasta_by_len_desc "{output.purge}" "{output.sorted}"
        samtools faidx "{output.sorted}" >> "{log.out}" 2>> "{log.err}"

        touch "{output.done}"
        """
