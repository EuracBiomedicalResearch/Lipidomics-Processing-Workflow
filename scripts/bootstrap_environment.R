#!/usr/bin/env Rscript

# Create the project-local R library and lockfile on the first run. On later
# runs (including fresh clones), restore the exact package versions recorded
# in renv.lock.

expected_conda_env <- "cembio_eurac"
expected_r <- "4.6.0"
bioconductor_version <- "3.23"
renv_version <- "1.2.4"
msbackend_revision <- "d606a8fecb7fc12f850426cc98c96fcc1cc3f133"
msbackend_source_url <- paste0(
    "https://github.com/RforMassSpectrometry/MsBackendMetaboLights/archive/",
    msbackend_revision,
    ".tar.gz"
)
msio_revision <- "a669fa1303a023161581b7a16caed3ed20a43299"
msio_source_url <- paste0(
    "https://github.com/RforMassSpectrometry/MsIO/archive/",
    msio_revision, ".tar.gz")

if (!file.exists("environment.yml") || !file.exists("DESCRIPTION")) {
    stop(
        "Run this script from the repository root:\n",
        "  Rscript scripts/bootstrap_environment.R",
        call. = FALSE
    )
}

active_conda_env <- Sys.getenv("CONDA_DEFAULT_ENV", unset = "")
if (!identical(active_conda_env, expected_conda_env)) {
    stop(
        "Expected Conda environment '", expected_conda_env,
        "', but CONDA_DEFAULT_ENV is '", active_conda_env, "'.\n",
        "Activate it with: conda activate ", expected_conda_env,
        call. = FALSE
    )
}

if (!identical(as.character(getRversion()), expected_r)) {
    stop(
        "Expected R ", expected_r, ", but found R ", getRversion(), ".\n",
        "Recreate the environment from environment.yml.",
        call. = FALSE
    )
}

Sys.setenv(DOWNLOAD_STATIC_LIBV8=1)
options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    Ncpus = max(1L, min(4L, parallel::detectCores(logical = FALSE)))
)

# V8 (via jsonvalidate and the alabaster packages) can use its upstream static
# Linux library. This avoids requiring distro-specific libnode/libv8 headers
# alongside the Conda toolchain. An explicitly supplied setting takes precedence.
if (identical(Sys.info()[["sysname"]], "Linux")) {
    Sys.setenv(DOWNLOAD_STATIC_LIBV8 = Sys.getenv(
        "DOWNLOAD_STATIC_LIBV8", unset = "1"
    ))
}

# Conda does not yet publish renv for R 4.6. Bootstrap the exact release from
# CRAN, after which renv manages the project library in the usual way.
has_expected_renv <- requireNamespace("renv", quietly = TRUE) &&
    identical(as.character(utils::packageVersion("renv")), renv_version)

if (!has_expected_renv) {
    renv_filename <- paste0("renv_", renv_version, ".tar.gz")
    renv_urls <- c(
        paste0("https://cloud.r-project.org/src/contrib/", renv_filename),
        paste0(
            "https://cloud.r-project.org/src/contrib/Archive/renv/",
            renv_filename
        )
    )
    message("Installing renv ", renv_version, " from CRAN...")
    for (renv_url in renv_urls) {
        install_result <- try(
            utils::install.packages(renv_url, repos = NULL, type = "source"),
            silent = TRUE
        )
        has_expected_renv <- requireNamespace("renv", quietly = TRUE) &&
            identical(
                as.character(utils::packageVersion("renv")),
                renv_version
            )
        if (has_expected_renv) break
    }
}

if (!has_expected_renv) {
    stop(
        "Could not install renv ", renv_version,
        " from either the current or archived CRAN location.",
        call. = FALSE
    )
}

if (file.exists("renv.lock")) {
    message("Restoring the R package library from renv.lock...")
    renv::restore(prompt = FALSE)
    message("Environment restored. Run: Rscript scripts/check_environment.R")
    quit(save = "no", status = 0L)
}

message(
    "No renv.lock found; initializing Bioconductor ",
    bioconductor_version,
    " and creating the first lockfile..."
)
if (!file.exists("renv/activate.R")) {
    renv::init(
        bare = TRUE,
        bioconductor = bioconductor_version,
        restart = FALSE
    )
} else {
    # Support a clean retry after an interrupted first installation.
    renv::settings$bioconductor.version(bioconductor_version)
}

# Install BiocManager from CRAN first, then expose the complete Bioconductor
# repository set before renv resolves any Bioconductor package names.
renv::install("BiocManager")
bioconductor_repos <- BiocManager::repositories(
    version = bioconductor_version
)
options(repos = bioconductor_repos)

cran_packages <- c(
    "dbplyr", "dplyr", "enviPat", "ggfortify", "ggplot2", "ggVennDiagram",
    "gridExtra", "knitr", "matrixStats", "pander", "pheatmap", "RColorBrewer",
    "Rcpp", "readxl", "rmarkdown", "RSQLite", "scam", "tictoc", "UpSetR",
    "vioplot", "writexl", "openxlsx", "quarto"
)

bioconductor_packages <- c(
    "alabaster.se", "AnnotationHub", "BiocFileCache", "BiocParallel",
    "CompoundDb", "limma", "MetaboAnnotation", "MetaboCoreUtils",
    "MsBackendSql", "MsExperiment", "Spectra", "SummarizedExperiment", "xcms"
)

renv::install(cran_packages, repos = bioconductor_repos)
renv::install(bioconductor_packages, repos = bioconductor_repos)
renv::install(msbackend_source_url)
renv::install(msio_source_url)
renv::snapshot(prompt = FALSE)

message(
    "Environment created and renv.lock written.\n",
    "Run: Rscript scripts/check_environment.R"
)
