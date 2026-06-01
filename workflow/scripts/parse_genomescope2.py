#!/usr/bin/env python3
import re
import sys
import csv
from pathlib import Path
from statistics import median

# -------------------------------
# Biological / heuristic filters
# -------------------------------
#
# These are designed to be species-agnostic and work for:
# - yeasts (~7–15 Mb genomes)
# - filamentous fungi (~25–35 Mb genomes)
#
# Instead of hard-coding absolute genome sizes, we:
#   1) compute the median genome_len across all models
#   2) keep only models within a reasonable factor of that median
#
# This reflects the GenomeScope2 behavior that the correct models
# tend to cluster near the true genome size, while artefactual
# models often give much smaller/larger sizes (e.g. 3–6 Mb or 50+ Mb).
#
# We also exclude models with unrealistically high heterozygosity
# (>20%), which is extremely unlikely for yeasts/fungi and usually
# indicates the model is explaining repeats/collapsed peaks as het.
#
# References:
# - Ranallo-Benavidez et al. 2020, Nat Commun 11:1432 (GenomeScope2)
# - Peter et al. 2018, Nature 556:339–344 (S. cerevisiae diversity)
# - Gladieux et al. 2014, Fungal Genet Biol 66:13–23 (fungal het ranges)

GENOME_LEN_FACTOR_MIN = 0.6   # keep models with len >= 60% of median
GENOME_LEN_FACTOR_MAX = 1.6   # and <= 160% of median
MAX_HET = 0.20                # 20% heterozygosity cap (very generous)

# -------------------------------
# Regexes for parsing
# -------------------------------

RE_HEADER = re.compile(r"GenomeScope analyzing .* p=(\d+)\s+k=(\d+)\s+outdir=(\S+)")
RE_MODEL  = re.compile(
    r"Model converged\s+het:(?P<het>[0-9.eE+-]+)\s+kcov:(?P<kcov>[0-9.eE+-]+)\s+err:(?P<err>[0-9.eE+-]+)\s+model fit:(?P<fit>[0-9.eE+-]+)\s+len:(?P<len>[0-9.eE+-]+)"
)
RE_HET_PCT = re.compile(r"heterozygosity:\s*([0-9.]+)\s*%")

def parse_one(log_path: Path):
    text = log_path.read_text(errors="replace").splitlines()
    rows = []
    current = None

    last_het_pct = None
    for line in text:
        m = RE_HET_PCT.search(line)
        if m:
            last_het_pct = float(m.group(1)) / 100.0

    for line in text:
        h = RE_HEADER.search(line)
        if h:
            current = {
                "ploidy": int(h.group(1)),
                "k": int(h.group(2)),
                "outdir": h.group(3),
                "het": None,
                "kcov": None,
                "err": None,
                "model_fit": None,
                "genome_len": None,
                "source_log": str(log_path),
            }
            continue

        m = RE_MODEL.search(line)
        if m and current:
            current["het"] = float(m.group("het"))
            current["kcov"] = float(m.group("kcov"))
            current["err"] = float(m.group("err"))
            current["model_fit"] = float(m.group("fit"))
            current["genome_len"] = float(m.group("len"))
            rows.append(current)
            current = None

    # If GenomeScope printed only "heterozygosity: X%" without "Model converged het:"
    # (rare), we still return the percent value as het if needed.
    for r in rows:
        if (r["het"] is None) and (last_het_pct is not None):
            r["het"] = last_het_pct

    return rows

def main():
    if len(sys.argv) < 3:
        print("Usage: parse_genomescope2.py <out.tsv> <log1> [log2 ...]", file=sys.stderr)
        sys.exit(2)

    out_tsv = Path(sys.argv[1])
    logs = [Path(p) for p in sys.argv[2:]]

    all_rows = []
    for lp in logs:
        all_rows.extend(parse_one(lp))

    # Filter to rows that have model_fit
    valid = [r for r in all_rows if r.get("model_fit") is not None]

    fields = ["ploidy","k","outdir","het","kcov","err","model_fit","genome_len","source_log","best"]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    if not valid:
        # Write empty TSV with header so Snakemake still has output
        with out_tsv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
        sys.exit(0)

    # --------------------------------------
    # Biological / heuristic filtering
    # --------------------------------------

    # Use median genome_len as a robust estimate of the "true" genome size
    # and discard models that are large outliers.
    genome_lens = [r["genome_len"] for r in valid if r.get("genome_len") is not None]
    if genome_lens:
        med_len = median(genome_lens)
    else:
        med_len = None

    def is_biologically_plausible(r):
        g = r.get("genome_len")
        h = r.get("het")
        if g is None or h is None:
            return False

        # Heterozygosity filter: discard models with absurdly high het
        if h > MAX_HET:
            return False

        # Genome length filter: relative to median, not absolute size
        if med_len is not None:
            if not (GENOME_LEN_FACTOR_MIN * med_len <= g <= GENOME_LEN_FACTOR_MAX * med_len):
                return False

        return True

    plausible = [r for r in valid if is_biologically_plausible(r)]

    # If we have at least one biologically plausible model, select best among those.
    # Otherwise, fall back to purely statistical choice (lowest model_fit).
    candidates = plausible if plausible else valid

    # Choose best = lowest model_fit among candidates
    best_row = min(candidates, key=lambda r: r["model_fit"])

    with out_tsv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(valid, key=lambda x: x["ploidy"]):
            r2 = dict(r)
            r2["best"] = "YES" if (r is best_row) else ""
            w.writerow(r2)

if __name__ == "__main__":
    main()