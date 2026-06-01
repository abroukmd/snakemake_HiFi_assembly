# rules/ploidyfrost.smk
import os
import csv

PROJECT_DIR = config["project_dir"]

# ----------------------------
# PloidyFrost install settings
# ----------------------------
PF_REPO   = config.get("ploidyfrost", {}).get("repo", "https://github.com/CMB-BNU/PloidyFrost")
PF_COMMIT = config.get("ploidyfrost", {}).get("commit", "")

PF_PREFIX = os.path.join(PROJECT_DIR, "resources", "ploidyfrost")
PF_SRC    = os.path.join(PF_PREFIX, "src")
PF_BUILD  = os.path.join(PF_PREFIX, "build")
PF_DONE   = os.path.join(PF_PREFIX, ".installed")

PF_BIN      = os.path.join(PF_BUILD, "PloidyFrost")
PF_BIFROST  = os.path.join(PF_SRC, "bifrost", "src", "Bifrost")
PF_FILTER_R = os.path.join(PF_SRC, "script", "Filter.R")
PF_DRAW_R   = os.path.join(PF_SRC, "script", "Drawfreq.R")

# ----------------------------
# Input reads from samples.tsv
# ----------------------------
SAMPLES_TSV = (
    config.get("sample_sheet")
    or config.get("samples_tsv")
    or os.path.join(PROJECT_DIR, "samples.tsv")
)

def load_samples_table(tsv):
    d = {}
    with open(tsv, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            s = row["sample"]
            d[s] = {"fq1": row["fq1"], "fq2": row["fq2"]}
    return d

SAMPLES_TABLE = load_samples_table(SAMPLES_TSV)

def fq1_of(wc):
    return SAMPLES_TABLE[wc.sample]["fq1"]

def fq2_of(wc):
    return SAMPLES_TABLE[wc.sample]["fq2"]

# ============================================================
# 1) INSTALL RULE
# ============================================================
rule install_ploidyfrost:
    output:
        done=PF_DONE
    params:
        repo=PF_REPO,
        commit=PF_COMMIT,
        prefix=PF_PREFIX,
        src=PF_SRC,
        build=PF_BUILD,
        bin=PF_BIN
    conda:
        "../envs/ploidyfrost_env.yaml"
    log:
        os.path.join(PF_PREFIX, "install.log")
    threads: 1
    shell:
        r"""
        set -euo pipefail
        set -x

        mkdir -p "{params.prefix}"
        touch "{log}"
        exec > >(tee -a "{log}") 2>&1

        echo "[INFO] repo={params.repo}"
        echo "[INFO] commit={params.commit}"
        echo "[INFO] prefix={params.prefix}"
        echo "[INFO] src={params.src}"
        echo "[INFO] build={params.build}"
        echo "[INFO] bin={params.bin}"
        echo "[INFO] CONDA_PREFIX=$CONDA_PREFIX"

        # If already installed, do nothing
        if [[ -f "{output.done}" ]]; then
            echo "[INFO] already installed: {output.done}"
            exit 0
        fi

        # Clone once (with submodules)
        if [[ ! -d "{params.src}/.git" ]]; then
            rm -rf "{params.src}" || true
            git clone --recursive "{params.repo}" "{params.src}"
        fi

        # Optional pin
        if [[ -n "{params.commit}" ]]; then
            (cd "{params.src}" && git fetch --all --tags)
            (cd "{params.src}" && git checkout -f "{params.commit}")
            (cd "{params.src}" && git submodule update --init --recursive)
        fi

        # Patch CMakeLists.txt safely (no fragile regex)
        python - <<'PY'
from pathlib import Path
import re

cml = Path(r"{params.src}") / "CMakeLists.txt"
lines = cml.read_text().splitlines(True)

out = []
for line in lines:
    # bump cmake_minimum_required to >= 3.5 (CMake >=3.28 requires this)
    m = re.match(r'\s*cmake_minimum_required\s*\(\s*VERSION\s+([0-9]+)\.([0-9]+)(\.[0-9]+)?\s*\)\s*', line)
    if m:
        out.append("cmake_minimum_required(VERSION 3.5)\n")
        continue

    # authors recommended replacements
    if re.match(r'^\s*find_library\s*\(\s*pthread\s+REQUIRED\s*\)\s*$', line):
        out.append("find_package(Threads)\n")
        continue

    if re.match(r'^\s*find_library\s*\(\s*z\s+REQUIRED\s*\)\s*$', line):
        out.append("find_package(ZLIB)\n")
        continue

    out.append(line)

txt = "".join(out)

# target_link_libraries line: pthread -> ${CMAKE_THREAD_LIBS_INIT}
# keep "z" as authors suggest; if this later fails at link time, we can switch to ${ZLIB_LIBRARIES}
txt = re.sub(
    r'target_link_libraries\s*\(\s*PloidyFrost\s+pthread\s+z\s+',
    'target_link_libraries(PloidyFrost ${{CMAKE_THREAD_LIBS_INIT}} z ',
    txt
)

cml.write_text(txt)
print("[INFO] Patched", cml)
PY

        mkdir -p "{params.build}"
        cd "{params.build}"

        rm -f CMakeCache.txt
        rm -rf CMakeFiles

        cmake -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
              "{params.src}"

        make -j 1

        test -x "{params.bin}" || (echo "[ERROR] Missing binary: {params.bin}" && exit 1)
        "{params.bin}" --help | head -n 20 || true

        touch "{output.done}"
        """

# ============================================================
# 2) RUN RULE
# ============================================================
rule ploidyfrost:
    input:
        tool=PF_DONE,
        fq1=fq1_of,
        fq2=fq2_of
    output:
        done=outpath("Assemblies/{sample}/QC/PloidyFrost/ploidyfrost.done"),
        ploidy=outpath("Assemblies/{sample}/QC/PloidyFrost/{sample}_ploidy.tsv")
    log:
        out=outpath("Assemblies/{sample}/QC/PloidyFrost/ploidyfrost.log")
    conda:
        "../envs/ploidyfrost_env.yaml"
    threads: 25
    shell:
        r"""
        set -euo pipefail
        set -x

        OUTDIR="$(dirname "{output.done}")"
        mkdir -p "$OUTDIR"
        cd "$OUTDIR"

        LIST="{wildcards.sample}.list"
        printf "%s\n" "{input.fq1}" "{input.fq2}" > "$LIST"

        while read -r f; do
            [[ -s "$f" ]] || (echo "[ERROR] Missing/empty FASTQ: $f" >> "{log.out}" && exit 1)
        done < "$LIST"

        mkdir -p kmc_tmp

        kmc -ci1 -cs10000 -k25 -t{threads} @"$LIST" "{wildcards.sample}_db" kmc_tmp >> "{log.out}" 2>&1
        kmc_tools transform "{wildcards.sample}_db" histogram hist >> "{log.out}" 2>&1

        lower_threshold=$("{PF_BIN}" cutoffL hist 2>> "{log.out}")
        echo "[INFO] cutoffL lower_threshold=$lower_threshold" >> "{log.out}"

        kmc_tools -t{threads} filter -hm "{wildcards.sample}_db" @"$LIST" -ci${{lower_threshold}} "{wildcards.sample}-filtered.fq" >> "{log.out}" 2>&1

        "{PF_BIFROST}" build -i -d -k 25 -v -r "{wildcards.sample}-filtered.fq" -o "{wildcards.sample}_dbg" -t {threads} >> "{log.out}" 2>&1

        "{PF_BIN}" -g "{wildcards.sample}_dbg.gfa" -d "{wildcards.sample}_db" -t {threads} -v -o "{wildcards.sample}" -h hist >> "{log.out}" 2>&1

        cd PloidyFrost_output

        Rscript "{PF_FILTER_R}" -i "{wildcards.sample}" -o "{wildcards.sample}_filtered" -n 6 -s 11 -q 0.05 >> "{log.out}" 2>&1
        "{PF_BIN}" model -g "{wildcards.sample}_filtered_allele_frequency.txt" -l 2 -u 10 -q 0.05 -o gmm >> "{log.out}" 2>&1
        Rscript "{PF_DRAW_R}" -f "{wildcards.sample}_filtered_allele_frequency.txt" -t "{wildcards.sample}_ploidy" -p 2 -o "{wildcards.sample}_histogram" >> "{log.out}" 2>&1

        # Placeholder output (replace when you decide canonical output)
        echo -e "{wildcards.sample}\tNA" > "{output.ploidy}"

        cd "$OUTDIR"
        touch "{output.done}"
        """
