#############################################
# rules/ragtag.smk
# Purged-only RagTag
#############################################

def ragtag_outputs(sample_wildcard, asm_label, rootdir):
    ordered_ext = "hifiasm" if asm_label == "hifiasm" else "lja"
    suffix = "purged" if asm_label == "hifiasm" else "lja"
    return {
        "fasta": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.fasta"),
        "agp": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.agp"),
        "fai": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.fasta.fai"),
        "ordered": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.{ordered_ext}.reforder.fasta"),
        "ordered_fai": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.{ordered_ext}.reforder.fasta.fai"),
        "compat": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.reforder.fasta"),
        "compat_fai": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/purged-{asm_label}_ragtag/ragtag.scaffold.reforder.fasta.fai"),
        "final_purged": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.purged.final.fasta"),
        "final_purged_fai": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.purged.final.fasta.fai"),
        "final_purged_map": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.purged.name_map.tsv"),
        "final_hap": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.hap.final.fasta"),
        "final_hap_fai": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.hap.final.fasta.fai"),
        "final_hap_map": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/FINAL_{sample_wildcard}/{sample_wildcard}.hap.name_map.tsv"),
        "done": outpath(f"Assemblies/{sample_wildcard}/{rootdir}/{sample_wildcard}_ragtag_{suffix}.done"),
    }

rule ragtag_purged:
    input:
        ref=lambda wc: sample_refs[wc.sample],
        purge=lambda wc: outpath(f"Assemblies/{wc.sample}/PURGE_DUPS/purged.hifiasm.fa"),
        hap=lambda wc: outpath(f"Assemblies/{wc.sample}/PURGE_DUPS/hap.hifiasm.fa"),
    output:
        **ragtag_outputs("{sample}", "hifiasm", "RAGTAG")
    threads: 32
    params:
        asm_label="hifiasm",
        threshold=32
    log:
        out=outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.log"),
        err=outpath("Assemblies/{sample}/RAGTAG/{sample}-ragtag-purged.err")
    conda:
        "../envs/ragtag.yaml"
    shell:
        r'''
        set -euo pipefail

        PURGED_DIR="$(dirname "{output.fasta}")"
        FINAL_DIR="$(dirname "{output.final_purged}")"
        mkdir -p "$PURGED_DIR" "$FINAL_DIR"

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
            local tmp_tab
            tmp_oneline=$(mktemp -p "$PURGED_DIR" hifiasm.oneline.XXXXXX.fa)
            tmp_tab=$(mktemp -p "$PURGED_DIR" hifiasm.tab.XXXXXX.tsv)
            fasta_to_oneline "$in" "$tmp_oneline"
            awk 'BEGIN{{OFS="\t"}}
                NR % 2 == 1 {{ hdr=$0; next }}
                NR % 2 == 0 {{ print hdr, length($0), $0 }}
            ' "$tmp_oneline" > "$tmp_tab"
            sort -t $'\t' -k2,2nr "$tmp_tab" | awk -F '\t' '{{ print $1; print $3 }}' > "$out"
            rm -f "$tmp_oneline" "$tmp_tab"
        }}

        reorder_scaffolds_by_ref_keep_unmatched() {{
            local ref="$1"
            local agp="$2"
            local fasta="$3"
            local out_sorted="$4"
            local dir
            dir="$(dirname "$out_sorted")"
            mkdir -p "$dir"

            local reforder="$dir/ref.order.txt"
            local matched_order="$dir/scaffold.order.matched.txt"
            local fasta_all="$dir/scaffold.order.all.txt"
            local unmatched_set="$dir/scaffold.order.unmatched.set.txt"
            local final_order="$dir/scaffold.order.final.txt"
            local present_ids="$dir/scaffold.present.txt"
            local filtered_order="$dir/scaffold.order.filtered.txt"

            if [[ "$ref" == *.gz ]]; then
                zgrep '^>' "$ref" | cut -d' ' -f1 | sed 's/^>//' > "$reforder"
            else
                grep '^>' "$ref" | cut -d' ' -f1 | sed 's/^>//' > "$reforder"
            fi

            awk '
                NR==FNR {{ ord[$1]=NR; next }}
                /^#/ {{ next }}
                $1 != prev {{
                    name=$1
                    sub(/_RagTag$/, "", name)
                    if (name in ord) print ord[name] "\t" $1
                    prev=$1
                }}
            ' "$reforder" "$agp" | sort -k1,1n | cut -f2 > "$matched_order"

            grep '^>' "$fasta" | sed 's/^>//' | cut -d' ' -f1 > "$fasta_all"
            comm -23 <(sort "$fasta_all") <(sort "$matched_order") > "$unmatched_set"
            cp -f "$matched_order" "$final_order"
            awk 'NR==FNR{{u[$1]=1; next}} ($1 in u){{print $1}}' "$unmatched_set" "$fasta_all" >> "$final_order"

            if [[ ! -s "$final_order" ]]; then
                cp -f "$fasta" "$out_sorted"
                return 0
            fi

            samtools faidx "$fasta" >> "{log.out}" 2>> "{log.err}"
            cut -f1 "${{fasta}}.fai" > "$present_ids"
            grep -Fxf "$present_ids" "$final_order" > "$filtered_order" || true

            if [[ ! -s "$filtered_order" ]]; then
                echo "[WARN] No valid reordered sequence IDs found; keeping original fasta order" >> "{log.out}"
                cp -f "$fasta" "$out_sorted"
                return 0
            fi

            : > "$out_sorted"
            while IFS= read -r seqid; do
                [[ -n "$seqid" ]] || continue
                if ! samtools faidx "$fasta" "$seqid" >> "$out_sorted" 2>> "{log.err}"; then
                    echo "[WARN] Could not extract sequence '$seqid'; skipping" >> "{log.out}"
                fi
            done < "$filtered_order"

            if [[ ! -s "$out_sorted" ]]; then
                echo "[WARN] Reordered fasta ended up empty; keeping original fasta order" >> "{log.out}"
                cp -f "$fasta" "$out_sorted"
            fi
        }}

        rename_fasta_to_sample_scaffs() {{
            local in="$1"
            local out="$2"
            local map="$3"
            local prefix="$4"
            local tmp_oneline
            tmp_oneline=$(mktemp -p "$FINAL_DIR" hifiasm.rename.XXXXXX.fa)
            fasta_to_oneline "$in" "$tmp_oneline"
            : > "$out"
            : > "$map"
            awk -v P="$prefix" -v OUT="$out" -v MAP="$map" '
                function clean(h, x) {{
                    x=h
                    sub(/^>/, "", x)
                    sub(/[ \t].*$/, "", x)
                    return x
                }}
                NR % 2 == 1 {{ old=clean($0); next }}
                NR % 2 == 0 {{
                    n++
                    new=P "_scaff" n
                    print old "\t" new >> MAP
                    print ">" new >> OUT
                    print $0 >> OUT
                }}
            ' "$tmp_oneline"
            rm -f "$tmp_oneline"
        }}

        require_nonempty() {{
            local f="$1"
            local label="$2"
            if [[ ! -s "$f" ]]; then
                echo "[ERROR] Missing or empty: $label ($f)" >> "{log.err}"
                exit 1
            fi
        }}

        REF_NSEQ=$(grep -c '^>' <(if [[ "{input.ref}" == *.gz ]]; then gunzip -c "{input.ref}"; else cat "{input.ref}"; fi) || true)
        echo "[INFO] Reference sequence count: $REF_NSEQ" > "{log.out}"

        echo "[INFO] RagTag scaffolding of purged hifiasm" >> "{log.out}"
        ragtag.py scaffold -t {threads} -o "$PURGED_DIR" "{input.ref}" "{input.purge}" >> "{log.out}" 2>> "{log.err}"

        require_nonempty "{output.fasta}" "purged ragtag scaffold fasta"
        samtools faidx "{output.fasta}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.fai}" "purged ragtag scaffold fai"

        if [[ "$REF_NSEQ" -gt {params.threshold} ]]; then
            echo "[INFO] Reference has > {params.threshold} sequences; ordering purged scaffolds by length" >> "{log.out}"
            sort_fasta_by_len_desc "{output.fasta}" "{output.ordered}"
        else
            echo "[INFO] Ordering purged scaffolds by reference order; unmatched scaffolds kept at end" >> "{log.out}"
            reorder_scaffolds_by_ref_keep_unmatched "{input.ref}" "{output.agp}" "{output.fasta}" "{output.ordered}"
        fi

        require_nonempty "{output.ordered}" "purged ordered fasta"
        cp -f "{output.ordered}" "{output.compat}"
        samtools faidx "{output.ordered}" >> "{log.out}" 2>> "{log.err}"
        samtools faidx "{output.compat}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.ordered_fai}" "purged ordered fasta fai"
        require_nonempty "{output.compat_fai}" "purged compat fasta fai"

        cp -f "{input.hap}" "{output.final_hap}"
        samtools faidx "{output.final_hap}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.final_hap}" "final hap fasta"
        require_nonempty "{output.final_hap_fai}" "final hap fasta fai"
        awk '/^>/{{gsub(/^>/,"",$1); print $1 "\t" "{wildcards.sample}_hap"}}' "{input.hap}" > "{output.final_hap_map}"

        rename_fasta_to_sample_scaffs "{output.ordered}" "{output.final_purged}" "{output.final_purged_map}" "{wildcards.sample}"
        require_nonempty "{output.final_purged}" "final purged fasta"
        samtools faidx "{output.final_purged}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.final_purged_fai}" "final purged fasta fai"

        touch "{output.done}"
        '''

rule ragtag_purged_lja:
    input:
        ref=lambda wc: sample_refs[wc.sample],
        purge=lambda wc: outpath(f"Assemblies/{wc.sample}/PURGE_DUPS_LJA/purged.lja.fa"),
        hap=lambda wc: outpath(f"Assemblies/{wc.sample}/PURGE_DUPS_LJA/hap.lja.fa"),
    output:
        **ragtag_outputs("{sample}", "lja", "RAGTAG_LJA")
    threads: 32
    params:
        asm_label="lja",
        threshold=32
    log:
        out=outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.log"),
        err=outpath("Assemblies/{sample}/RAGTAG_LJA/{sample}-ragtag-lja.err")
    conda:
        "../envs/ragtag.yaml"
    shell:
        r'''
        set -euo pipefail

        PURGED_DIR="$(dirname "{output.fasta}")"
        FINAL_DIR="$(dirname "{output.final_purged}")"
        mkdir -p "$PURGED_DIR" "$FINAL_DIR"

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
            local tmp_tab
            tmp_oneline=$(mktemp -p "$PURGED_DIR" lja.oneline.XXXXXX.fa)
            tmp_tab=$(mktemp -p "$PURGED_DIR" lja.tab.XXXXXX.tsv)
            fasta_to_oneline "$in" "$tmp_oneline"
            awk 'BEGIN{{OFS="\t"}}
                NR % 2 == 1 {{ hdr=$0; next }}
                NR % 2 == 0 {{ print hdr, length($0), $0 }}
            ' "$tmp_oneline" > "$tmp_tab"
            sort -t $'\t' -k2,2nr "$tmp_tab" | awk -F '\t' '{{ print $1; print $3 }}' > "$out"
            rm -f "$tmp_oneline" "$tmp_tab"
        }}

        reorder_scaffolds_by_ref_keep_unmatched() {{
            local ref="$1"
            local agp="$2"
            local fasta="$3"
            local out_sorted="$4"
            local dir
            dir="$(dirname "$out_sorted")"
            mkdir -p "$dir"

            local reforder="$dir/ref.order.txt"
            local matched_order="$dir/scaffold.order.matched.txt"
            local fasta_all="$dir/scaffold.order.all.txt"
            local unmatched_set="$dir/scaffold.order.unmatched.set.txt"
            local final_order="$dir/scaffold.order.final.txt"
            local present_ids="$dir/scaffold.present.txt"
            local filtered_order="$dir/scaffold.order.filtered.txt"

            if [[ "$ref" == *.gz ]]; then
                zgrep '^>' "$ref" | cut -d' ' -f1 | sed 's/^>//' > "$reforder"
            else
                grep '^>' "$ref" | cut -d' ' -f1 | sed 's/^>//' > "$reforder"
            fi

            awk '
                NR==FNR {{ ord[$1]=NR; next }}
                /^#/ {{ next }}
                $1 != prev {{
                    name=$1
                    sub(/_RagTag$/, "", name)
                    if (name in ord) print ord[name] "\t" $1
                    prev=$1
                }}
            ' "$reforder" "$agp" | sort -k1,1n | cut -f2 > "$matched_order"

            grep '^>' "$fasta" | sed 's/^>//' | cut -d' ' -f1 > "$fasta_all"
            comm -23 <(sort "$fasta_all") <(sort "$matched_order") > "$unmatched_set"
            cp -f "$matched_order" "$final_order"
            awk 'NR==FNR{{u[$1]=1; next}} ($1 in u){{print $1}}' "$unmatched_set" "$fasta_all" >> "$final_order"

            if [[ ! -s "$final_order" ]]; then
                cp -f "$fasta" "$out_sorted"
                return 0
            fi

            samtools faidx "$fasta" >> "{log.out}" 2>> "{log.err}"
            cut -f1 "${{fasta}}.fai" > "$present_ids"
            grep -Fxf "$present_ids" "$final_order" > "$filtered_order" || true

            if [[ ! -s "$filtered_order" ]]; then
                echo "[WARN] No valid reordered sequence IDs found; keeping original fasta order" >> "{log.out}"
                cp -f "$fasta" "$out_sorted"
                return 0
            fi

            : > "$out_sorted"
            while IFS= read -r seqid; do
                [[ -n "$seqid" ]] || continue
                if ! samtools faidx "$fasta" "$seqid" >> "$out_sorted" 2>> "{log.err}"; then
                    echo "[WARN] Could not extract sequence '$seqid'; skipping" >> "{log.out}"
                fi
            done < "$filtered_order"

            if [[ ! -s "$out_sorted" ]]; then
                echo "[WARN] Reordered fasta ended up empty; keeping original fasta order" >> "{log.out}"
                cp -f "$fasta" "$out_sorted"
            fi
        }}

        rename_fasta_to_sample_scaffs() {{
            local in="$1"
            local out="$2"
            local map="$3"
            local prefix="$4"
            local tmp_oneline
            tmp_oneline=$(mktemp -p "$FINAL_DIR" lja.rename.XXXXXX.fa)
            fasta_to_oneline "$in" "$tmp_oneline"
            : > "$out"
            : > "$map"
            awk -v P="$prefix" -v OUT="$out" -v MAP="$map" '
                function clean(h, x) {{
                    x=h
                    sub(/^>/, "", x)
                    sub(/[ \t].*$/, "", x)
                    return x
                }}
                NR % 2 == 1 {{ old=clean($0); next }}
                NR % 2 == 0 {{
                    n++
                    new=P "_scaff" n
                    print old "\t" new >> MAP
                    print ">" new >> OUT
                    print $0 >> OUT
                }}
            ' "$tmp_oneline"
            rm -f "$tmp_oneline"
        }}

        require_nonempty() {{
            local f="$1"
            local label="$2"
            if [[ ! -s "$f" ]]; then
                echo "[ERROR] Missing or empty: $label ($f)" >> "{log.err}"
                exit 1
            fi
        }}

        REF_NSEQ=$(grep -c '^>' <(if [[ "{input.ref}" == *.gz ]]; then gunzip -c "{input.ref}"; else cat "{input.ref}"; fi) || true)
        echo "[INFO] Reference sequence count: $REF_NSEQ" > "{log.out}"

        echo "[INFO] RagTag scaffolding of purged lja" >> "{log.out}"
        ragtag.py scaffold -t {threads} -o "$PURGED_DIR" "{input.ref}" "{input.purge}" >> "{log.out}" 2>> "{log.err}"

        require_nonempty "{output.fasta}" "purged ragtag scaffold fasta"
        samtools faidx "{output.fasta}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.fai}" "purged ragtag scaffold fai"

        if [[ "$REF_NSEQ" -gt {params.threshold} ]]; then
            echo "[INFO] Reference has > {params.threshold} sequences; ordering purged scaffolds by length" >> "{log.out}"
            sort_fasta_by_len_desc "{output.fasta}" "{output.ordered}"
        else
            echo "[INFO] Ordering purged scaffolds by reference order; unmatched scaffolds kept at end" >> "{log.out}"
            reorder_scaffolds_by_ref_keep_unmatched "{input.ref}" "{output.agp}" "{output.fasta}" "{output.ordered}"
        fi

        require_nonempty "{output.ordered}" "purged ordered fasta"
        cp -f "{output.ordered}" "{output.compat}"
        samtools faidx "{output.ordered}" >> "{log.out}" 2>> "{log.err}"
        samtools faidx "{output.compat}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.ordered_fai}" "purged ordered fasta fai"
        require_nonempty "{output.compat_fai}" "purged compat fasta fai"

        cp -f "{input.hap}" "{output.final_hap}"
        samtools faidx "{output.final_hap}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.final_hap}" "final hap fasta"
        require_nonempty "{output.final_hap_fai}" "final hap fasta fai"
        awk '/^>/{{gsub(/^>/,"",$1); print $1 "\t" "{wildcards.sample}_hap"}}' "{input.hap}" > "{output.final_hap_map}"

        rename_fasta_to_sample_scaffs "{output.ordered}" "{output.final_purged}" "{output.final_purged_map}" "{wildcards.sample}"
        require_nonempty "{output.final_purged}" "final purged fasta"
        samtools faidx "{output.final_purged}" >> "{log.out}" 2>> "{log.err}"
        require_nonempty "{output.final_purged_fai}" "final purged fasta fai"

        touch "{output.done}"
        '''