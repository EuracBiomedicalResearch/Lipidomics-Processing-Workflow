# Resolve script directory at source() time so .knn_impute_rcpp can find
# knn_impute.cpp regardless of the caller's working directory.
.script_dir <- local({
    for (i in seq_len(sys.nframe())) {
        f <- sys.frame(i)$ofile
        if (!is.null(f)) return(dirname(normalizePath(f)))
    }
    "."
})

# =============================================================================
# impute_block_knn — Three-tier block KNN imputation
# =============================================================================
# Adapted from MATLAB code by Santiago Angulo (reemplazarNaN3.m).
# Features are assigned to one of three blocks by proportion of missing values:
#   Block 1 (prop_na <= a):        direct KNN imputation
#   Block 2 (a < prop_na < b):     KNN with zero-filled context features
#   Block 3 (prop_na >= b):        replaced with zeros

#' @param mat Numeric matrix (features x samples), data.frame, or
#'   SummarizedExperiment.
#' @param a Maximum NA proportion for direct KNN (Block 1). 0 <= a < b.
#' @param b Minimum NA proportion for zero replacement (Block 3). b <= 1.
#' @param k Number of nearest neighbours.
#' @param ... For SummarizedExperiment: assay.type (integer or character,
#'   default 1L).
#' @return Object of same class as mat with NAs replaced.
setGeneric(
    "impute_block_knn",
    function(mat, a, b, k, ...) standardGeneric("impute_block_knn")
)

setMethod(
    "impute_block_knn", "matrix",
    function(mat, a, b, k, ...) .impute_block_knn_impl(mat, a, b, k)
)

setMethod(
    "impute_block_knn", "data.frame",
    function(mat, a, b, k, ...) {
        m <- as.matrix(mat)
        if (!is.numeric(m))
            stop("data.frame must contain only numeric columns")
        as.data.frame(.impute_block_knn_impl(m, a, b, k))
    }
)

setMethod(
    "impute_block_knn", "SummarizedExperiment",
    function(mat, a, b, k, ..., assay.type = 1L) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
            stop("SummarizedExperiment package is required.")
        m <- SummarizedExperiment::assay(mat, assay.type)
        imputed <- .impute_block_knn_impl(m, a, b, k)
        SummarizedExperiment::assay(
            mat, assay.type, withDimnames = FALSE
        ) <- imputed
        mat
    }
)

.impute_block_knn_impl <- function(mat, a, b, k) {
    if (a > b || a < 0 || b > 1)
        stop("Invalid thresholds: need 0 <= a < b <= 1 (a=", a, ", b=", b, ")")
    n <- nrow(mat)
    m <- ncol(mat)
    if (m == 1 && any(is.na(mat)))
        stop("KNN imputation requires at least 2 samples (columns).")
    orig_dimnames <- dimnames(mat)
    mat <- t(mat)              # m x n: samples in rows, features in cols
    mat_indexed <- rbind(seq_len(n), mat)
    prop_na <- colSums(is.na(mat)) / m
    block1 <- mat_indexed[, prop_na <= a, drop = FALSE]
    block2 <- mat_indexed[, prop_na > a & prop_na < b, drop = FALSE]
    block3 <- mat_indexed[, prop_na >= b, drop = FALSE]
    block1 <- t(block1)
    block2 <- t(block2)
    block3 <- t(block3)
    if (nrow(block1) > 0) {
        idx1 <- block1[, 1]
        block1 <- block1[, 2:(m + 1), drop = FALSE]
    } else {
        idx1 <- numeric(0)
        block1 <- matrix(nrow = 0, ncol = m)
    }
    if (nrow(block2) > 0) {
        idx2 <- block2[, 1]
        block2 <- block2[, 2:(m + 1), drop = FALSE]
    } else {
        idx2 <- numeric(0)
        block2 <- matrix(nrow = 0, ncol = m)
    }
    if (nrow(block3) > 0) {
        idx3 <- block3[, 1]
        block3 <- block3[, 2:(m + 1), drop = FALSE]
    } else {
        idx3 <- numeric(0)
        block3 <- matrix(nrow = 0, ncol = m)
    }
    block3_zeros <- block3; block3_zeros[is.na(block3_zeros)] <- 0
    block2_zeros <- block2; block2_zeros[is.na(block2_zeros)] <- 0
    mat_ctx <- rbind(block1, block2_zeros, block3_zeros)
    mat_ind <- rbind(block1, block2, block3_zeros)
    mat_imputed <- matrix(0, nrow = n, ncol = m)
    mat_work <- mat_ctx
    for (i in seq_len(m)) {
        if (nrow(mat_work) == 0) {
            mat_imputed[, i] <- mat_ind[, i]
            mat_imputed[is.na(mat_imputed[, i]), i] <- 0
            next
        }
        orig_col <- mat_work[, i]
        mat_work[, i] <- mat_ind[, i]
        if (!any(rowSums(is.na(mat_work)) == 0L))
            stop("All rows contain missing values. Unable to impute.")
        mat_imputed[, i] <- .knn_impute_rcpp(mat_work, k)[, i]
        mat_work[, i] <- orig_col
    }
    orig_order <- c(idx1, idx2, idx3)
    mat_ordered <- cbind(orig_order, mat_imputed)
    mat_ordered <- mat_ordered[order(orig_order), , drop = FALSE]
    result <- mat_ordered[, 2:(m + 1), drop = FALSE]
    dimnames(result) <- orig_dimnames
    result
}

.knn_impute_rcpp <- function(mat, k) {
    if (!requireNamespace("Rcpp", quietly = TRUE))
        stop("Rcpp package is required. Install with: install.packages('Rcpp')")
    cpp_file <- file.path(.script_dir, "knn_impute.cpp")
    if (!exists("knn_impute_rcpp", mode = "function", envir = .GlobalEnv))
        Rcpp::sourceCpp(cpp_file)
    # C++ expects samples in rows (n_samples x n_features)
    t(knn_impute_rcpp(t(mat), as.integer(k)))
}

# =============================================================================
# impute_block_knn_log — block KNN in log2(1+x) space
# =============================================================================

#' Impute with block KNN after log2(1+x) transformation
#'
#' @param mat Numeric matrix (features x samples).
#' @param a,b,k Passed to impute_block_knn.
#' @return Imputed matrix in the original (untransformed) scale.
impute_block_knn_log <- function(mat, a, b, k) {
    prop_na <- rowMeans(is.na(mat))
    block3_mask <- is.na(mat) & (prop_na >= b)
    imputed <- 2^impute_block_knn(log2(1 + mat), a, b, k) - 1
    for (i in which(rowSums(block3_mask) > 0)) {
        obs_min <- suppressWarnings(min(mat[i, ], na.rm = TRUE))
        if (!is.finite(obs_min)) next
        na_pos <- block3_mask[i, ]
        imputed[i, na_pos] <- runif(sum(na_pos), min = obs_min / 2, max = obs_min)
    }
    imputed
}

# =============================================================================
# impute_unidis — uniform distribution near per-feature minimum
# =============================================================================

#' Impute missing values using a uniform distribution near the feature minimum
#'
#' @param z Numeric vector with potential NA values.
#' @return Vector with NAs replaced by values drawn from runif(min/2, min).
impute_unidis <- function(z) {
    na <- is.na(z)
    if (any(na)) {
        min_val <- min(z, na.rm = TRUE)
        z[na] <- runif(sum(na), min = min_val / 2, max = min_val)
    }
    z
}
