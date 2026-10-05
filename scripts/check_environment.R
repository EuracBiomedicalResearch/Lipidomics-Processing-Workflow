#!/usr/bin/env Rscript

expected_conda_env <- "cembio_eurac"
expected_r <- "4.6.0"
expected_bioconductor <- "3.23"
expected_quarto <- "1.9.38"
expected_msstash <- "0.99.0"
expected_spectrastash <- "0.99.2"
expected_msexperimentstash <- "0.99.0"
expected_xcmsstash <- "0.97.2"

required_packages <- c(
    "alabaster.se", "AnnotationHub", "BiocFileCache", "BiocManager",
    "BiocParallel", "CompoundDb", "dbplyr", "dplyr", "enviPat",
    "ggfortify", "ggplot2", "ggVennDiagram", "gridExtra", "knitr", "limma",
    "matrixStats", "MetaboAnnotation", "MetaboCoreUtils",
    "MsBackendMetaboLights", "MsBackendSql", "MsExperiment",
    "MsExperimentStash", "MsStash", "pander", "pheatmap", "RColorBrewer",
    "Rcpp", "readxl", "rmarkdown", "RSQLite", "scam", "Spectra", "SpectraStash",
    "SummarizedExperiment", "tictoc", "UpSetR", "vioplot", "writexl", "xcms",
    "xcmsStash"
)

failures <- character()

if (!identical(Sys.getenv("CONDA_DEFAULT_ENV", unset = ""), expected_conda_env)) {
    failures <- c(failures, paste0("Conda environment is not ", expected_conda_env))
}

if (!identical(as.character(getRversion()), expected_r)) {
    failures <- c(
        failures,
        paste0("R is ", getRversion(), "; expected ", expected_r)
    )
}

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
    failures <- c(
        failures,
        paste("Missing R packages:", paste(missing_packages, collapse = ", "))
    )
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
    actual_bioconductor <- as.character(BiocManager::version())
    if (!identical(actual_bioconductor, expected_bioconductor)) {
        failures <- c(
            failures,
            paste0(
                "Bioconductor is ", actual_bioconductor,
                "; expected ", expected_bioconductor
            )
        )
    }
}

if (requireNamespace("MsIO", quietly = TRUE)) {
    actual_msio <- as.character(utils::packageVersion("MsIO"))
    if (!identical(actual_msio, expected_msio)) {
        failures <- c(
            failures,
            paste0("MsIO is ", actual_msio, "; expected ", expected_msio)
        )
    }
}

quarto_path <- Sys.which("quarto")
if (!nzchar(quarto_path)) {
    failures <- c(failures, "Quarto is not on PATH")
} else {
    actual_quarto <- system2(quarto_path, "--version", stdout = TRUE)
    if (!identical(actual_quarto[[1L]], expected_quarto)) {
        failures <- c(
            failures,
            paste0(
                "Quarto is ", actual_quarto[[1L]],
                "; expected ", expected_quarto
            )
        )
    }
}

if (!length(failures) && requireNamespace("Rcpp", quietly = TRUE)) {
    message("Compiling R/knn_impute.cpp to verify the native toolchain...")
    compile_result <- try(
        Rcpp::sourceCpp("R/knn_impute.cpp", rebuild = TRUE, verbose = FALSE),
        silent = TRUE
    )
    if (inherits(compile_result, "try-error")) {
        failures <- c(
            failures,
            paste("Rcpp compilation failed:", as.character(compile_result))
        )
    }
}

if (length(failures)) {
    message("Environment check failed:")
    for (failure in failures) message("  - ", failure)
    quit(save = "no", status = 1L)
}

message("Environment check passed:")
message("  Conda environment: ", expected_conda_env)
message("  R: ", R.version.string)
message("  Bioconductor: ", expected_bioconductor)
message("  Quarto: ", expected_quarto)
message("  MsIO: ", expected_msio)
message("  R packages: ", length(required_packages), " available")
message("  Native compilation: OK")
