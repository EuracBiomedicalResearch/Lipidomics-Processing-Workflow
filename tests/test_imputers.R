# Smoke-test and sanity-check for imputers.R
# Loads pilot POS data, runs all three imputers, asserts basic invariants,
# and prints a comparison summary table.

rm(list = ls())

suppressPackageStartupMessages({
    library(alabaster.se)
    library(SummarizedExperiment)
    library(methods)
})

source("R/imputers.R")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

res <- readObject(
    "applications/pilot_study/positive/objects/pilot_final_res_pos"
)
mat <- as.matrix(assay(res, "ISnorm_filled"))

cat("Input:", nrow(mat), "features x", ncol(mat), "samples,",
    sum(is.na(mat)), "NAs\n\n")

na_mask <- is.na(mat)
stopifnot(any(na_mask))   # test is only useful if there are NAs to fill

# ---------------------------------------------------------------------------
# Helper: assert invariants common to all imputers
# ---------------------------------------------------------------------------

check_imputer <- function(result, name, obs_tol = 0) {
    obs_err <- max(abs(result[!na_mask] - mat[!na_mask]))
    stopifnot(
        "no NAs in output"          = !any(is.na(result)),
        "same shape as input"       = identical(dim(result), dim(mat)),
        "observed values preserved" = obs_err <= obs_tol
    )
    imputed_vals <- result[na_mask]
    cat(sprintf(
        "%-25s  n_imputed=%d  mean=%.2f  min=%.2f  max=%.2f\n",
        name, sum(na_mask),
        mean(imputed_vals), min(imputed_vals), max(imputed_vals)
    ))
    invisible(result)
}

# ---------------------------------------------------------------------------
# impute_unidis
# ---------------------------------------------------------------------------

res_unidis <- t(apply(mat, MARGIN = 1, impute_unidis))
check_imputer(res_unidis, "impute_unidis")

# Each imputed value must lie in [min_i/2, min_i] for its feature
for (i in which(rowSums(na_mask) > 0)) {
    obs_min  <- min(mat[i, ], na.rm = TRUE)
    imputed  <- res_unidis[i, na_mask[i, ]]
    stopifnot(all(imputed >= obs_min / 2), all(imputed <= obs_min))
}
cat("  [impute_unidis] bound check passed\n\n")

# ---------------------------------------------------------------------------
# impute_block_knn
# ---------------------------------------------------------------------------

res_knn <- impute_block_knn(mat, a = 0.4, b = 0.7, k = 5)
check_imputer(res_knn, "impute_block_knn")
cat("\n")

# ---------------------------------------------------------------------------
# impute_block_knn_log
# ---------------------------------------------------------------------------

res_knn_log <- impute_block_knn_log(mat, a = 0.4, b = 0.7, k = 5)
check_imputer(res_knn_log, "impute_block_knn_log", obs_tol = 1e-6)

stopifnot("output non-negative (back-transform invariant)" = all(res_knn_log >= 0))
cat("  [impute_block_knn_log] positivity check passed\n")

# round-trip: matrix with no NAs should be (approximately) unchanged
mat_complete <- mat
mat_complete[na_mask] <- 1   # fill NAs with arbitrary positive value
rt <- impute_block_knn_log(mat_complete, a = 0.4, b = 0.7, k = 5)
stopifnot(
    "round-trip: no NAs in output"    = !any(is.na(rt)),
    "round-trip: values approx equal" = max(abs(rt - mat_complete)) < 1e-6
)
cat("  [impute_block_knn_log] round-trip check passed\n\n")

cat("All checks passed.\n")
