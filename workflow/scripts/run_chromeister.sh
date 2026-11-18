#!/usr/bin/env bash
# Run one comparison with isolated reference preparation and real PNG outputs.
set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo "Usage: run_chromeister.sh query ref matrix png score done stdout stderr allow_empty_hap" >&2
    exit 2
fi

query=$(realpath "$1")
ref=$(realpath "$2")
matrix=$3
png=$4
score=$5
done_file=$6
stdout_log=$7
stderr_log=$8
allow_empty_hap=$9

mkdir -p "$(dirname "$matrix")" "$(dirname "$stdout_log")" "$(dirname "$stderr_log")"
# Require artifacts from this invocation, including the native plotting filename
# when the declared PNG uses a different path.
rm -f -- "$done_file" "$matrix" "$png" "$score" "$matrix.filt.png"
: > "$stdout_log"
: > "$stderr_log"
test -s "$query"
test -s "$ref"

work_dir=$(mktemp -d "$(dirname "$matrix")/.chromeister.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

if [[ "$allow_empty_hap" == 1 ]] && grep -q '^>hap_empty\([[:space:]]\|$\)' "$query"; then
    # Optional haplotig comparisons have an explicit, visible skipped result.
    printf '# SKIPPED: no haplotigs (hap_empty)\n' > "$matrix"
    printf 'SKIPPED hap_empty\n' > "$score"
    Rscript -e 'png(commandArgs(TRUE)[1]); plot.new(); text(0.5, 0.5, "SKIPPED: no haplotigs"); invisible(dev.off())' \
        "$png" >> "$stdout_log" 2>> "$stderr_log"
else
    if [[ "$ref" == *.gz ]]; then
        gzip -cd -- "$ref" > "$work_dir/reference.fa"
        ref="$work_dir/reference.fa"
        test -s "$ref"
    fi
    CHROMEISTER -query "$ref" -db "$query" -out "$matrix" -dimension 2000 \
        > "$stdout_log" 2> "$stderr_log"
    if [[ ! -s "$matrix" ]]; then
        echo "CHROMEISTER did not produce a fresh matrix: $matrix" >> "$stderr_log"
        exit 1
    fi
    Rscript "$(command -v compute_score.R)" "$matrix" 2000 > "$score" 2>> "$stderr_log"
    # Chromeister's compute_score.R writes <matrix>.filt.png.
    generated_png="$matrix.filt.png"
    if [[ ! -s "$generated_png" ]]; then
        echo "Chromeister scoring did not produce a fresh PNG: $generated_png" >> "$stderr_log"
        exit 1
    fi
    if [[ "$generated_png" != "$png" ]]; then
        mv -- "$generated_png" "$png"
    fi
fi

test -s "$score"
test -s "$png"
signature=$(od -An -tx1 -N8 "$png" | tr -d '[:space:]')
if [[ "$signature" != 89504e470d0a1a0a ]]; then
    echo "Invalid PNG output: $png" >> "$stderr_log"
    exit 1
fi
touch "$done_file"
