import os

PROJECT_DIR = config["project_dir"]

# ----------------------------
# GenomeScope2 install settings
# ----------------------------
GS2_REPO   = config.get("genomescope2", {}).get("repo", "https://github.com/tbenavi1/genomescope2.0.git")
GS2_COMMIT = config.get("genomescope2", {}).get("commit", "")

GS2_PREFIX = os.path.join(PROJECT_DIR, "resources", "genomescope2")
GS2_SRC    = os.path.join(GS2_PREFIX, "src")
GS2_RLIB   = os.path.join(GS2_PREFIX, "Rlib")
GS2_DONE   = os.path.join(GS2_PREFIX, ".installed")

# ----------------------------
# nQuire install settings
# ----------------------------
NQ_REPO   = config.get("nquire", {}).get("repo", "https://github.com/clwgg/nQuire")
NQ_COMMIT = config.get("nquire", {}).get("commit", "")

NQ_PREFIX = os.path.join(PROJECT_DIR, "resources", "nQuire")
NQ_SRC    = os.path.join(NQ_PREFIX, "src")
NQ_DONE   = os.path.join(NQ_PREFIX, ".installed")
NQ_BIN    = os.path.join(NQ_SRC, "nQuire")


rule install_genomescope2:
    output:
        done=GS2_DONE
    params:
        repo=GS2_REPO,
        commit=GS2_COMMIT,
        prefix=GS2_PREFIX,
        src=GS2_SRC,
        rlib=GS2_RLIB,
    conda:
        "../envs/genomescope2_env.yaml"
    log:
        os.path.join(GS2_PREFIX, "install.log")
    shell:
        r"""
        set -euo pipefail
        set -x

        mkdir -p "{params.prefix}" "{params.rlib}"
        # DO NOT truncate the log blindly; keep history
        touch "{log}"
        exec > >(tee -a "{log}") 2>&1

        echo "[INFO] repo={params.repo}"
        echo "[INFO] commit={params.commit}"
        echo "[INFO] prefix={params.prefix}"
        echo "[INFO] rlib={params.rlib}"
        echo "[INFO] which git:     $(command -v git     || echo MISSING)"
        echo "[INFO] which Rscript: $(command -v Rscript  || echo MISSING)"
        echo "[INFO] which R:       $(command -v R        || echo MISSING)"

        # clone once
        if [[ ! -d "{params.src}/.git" ]]; then
            rm -rf "{params.src}" || true
            git clone --depth 1 "{params.repo}" "{params.src}"
        fi

        # optional pin
        if [[ -n "{params.commit}" ]]; then
            (cd "{params.src}" && git fetch --all --tags)
            (cd "{params.src}" && git checkout -f "{params.commit}")
        fi

        export R_LIBS_USER="{params.rlib}"

        Rscript -e 'cat("[INFO] R:", as.character(getRversion()), "\n"); cat("[INFO] R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n"); cat("[INFO] .libPaths():\n", paste(.libPaths(), collapse="\n"), "\n", sep="")'

        (cd "{params.src}" && Rscript install.R)

        # verify
        Rscript -e 'stopifnot(requireNamespace("genomescope", quietly=TRUE)); cat("[OK] genomescope loadable\n")'

        touch "{output.done}"
        """

