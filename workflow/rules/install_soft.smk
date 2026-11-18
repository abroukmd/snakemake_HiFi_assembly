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

        # Resolve the requested revision to its immutable commit and verify HEAD.
        if [[ -n "{params.commit}" ]]; then
            (cd "{params.src}" && git fetch --all --tags)
            (cd "{params.src}" && git checkout -f "{params.commit}")
            test "$(git -C {params.src:q} rev-parse HEAD)" = \
                 "$(git -C {params.src:q} rev-parse {params.commit:q}^{{commit}})"
        fi

        export R_LIBS_USER="{params.rlib}"

        Rscript -e 'cat("[INFO] R:", as.character(getRversion()), "\n"); cat("[INFO] R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n"); cat("[INFO] .libPaths():\n", paste(.libPaths(), collapse="\n"), "\n", sep="")'

        # The pinned DESCRIPTION imports only these three R packages. Upstream's
        # install.R reinstalls them with dependencies=TRUE, pulling unrelated
        # optional spatial packages and their system libraries from CRAN.
        Rscript -e 'required <- c("argparse", "jsonlite", "minpack.lm"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing required Conda R packages: ", paste(missing, collapse=", "))'
        R CMD INSTALL --library={params.rlib:q} {params.src:q}

        # Verify the package is loaded from this project's Rlib and its package
        # identity/version agree with the checked-out, pinned source.
        Rscript -e 'expected <- normalizePath(Sys.getenv("R_LIBS_USER")); actual <- normalizePath(find.package("genomescope", lib.loc=expected)); stopifnot(dirname(actual) == expected); installed <- packageDescription("genomescope", lib.loc=expected); source <- read.dcf(file.path(commandArgs(TRUE)[1], "DESCRIPTION")); stopifnot(installed$Package == source[1,"Package"], installed$Version == source[1,"Version"]); stopifnot(requireNamespace("genomescope", lib.loc=expected, quietly=TRUE)); cat("[OK] genomescope", installed$Version, "loaded from", actual, "\n")' {params.src:q}

        git -C {params.src:q} rev-parse HEAD > {output.done:q}
        """

