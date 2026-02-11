#' imputeBlockKNN - Impute Missing Values Using Three-Tier Block KNN Strategy
#'
#' S4 generic function to replace missing values using a three-tier
#' k-nearest neighbors strategy based on the proportion of missing values
#' per variable (feature).
#'
#' @description
#' Adapted from MATLAB code written by Santiago Angulo
#' (original: reemplazarNaN3.m).
#'
#' The function assigns each variable to one of three blocks depending on
#' its proportion of missing values, then applies a different imputation
#' strategy per block:
#' \enumerate{
#'   \item If proportion of NA <= a: Apply KNN imputation directly
#'   \item If a < proportion of NA < b: Mixed strategy (KNN with
#'         zero-replaced context variables)
#'   \item If proportion of NA >= b: Replace all NA with zeros
#' }
#'
#' @param datos Input data. Can be a numeric \code{matrix}, a
#'   \code{data.frame}, or a \code{SummarizedExperiment}. For matrix and
#'   data.frame: rows = variables/features, columns = samples. For
#'   SummarizedExperiment: standard Bioconductor convention (features in
#'   rows, samples in columns).
#' @param a Scalar. Maximum proportion of NA for direct KNN (Block 1).
#'   Must be >= 0 and < b.
#' @param b Scalar. Minimum proportion of NA for zero replacement (Block 3).
#'   Must be <= 1 and > a.
#' @param k Integer. Number of nearest neighbors for KNN imputation.
#' @param ... Additional arguments passed to methods:
#'   \describe{
#'     \item{assay.type}{(SummarizedExperiment only) Integer or character
#'       selecting which assay to impute. Default \code{1L}.}
#'   }
#'
#' @return Object of same class as \code{datos} with NA values replaced.
#'
#' @author Original MATLAB code by Santiago Angulo. Translation by C.A. García
#' @export
#'
#' @examples
#' \dontrun{
#' # Matrix input
#' datos <- matrix(c(
#'     1.0, 2.0, NA, 4.0, 5.0,
#'     1.0, NA, NA, NA, 5.0,
#'     NA, NA, NA, NA, 5.0
#' ), nrow = 3, byrow = TRUE)
#'
#' result <- imputeBlockKNN(datos, a = 0.4, b = 0.7, k = 2)
#'
#' # SummarizedExperiment input
#' result_se <- imputeBlockKNN(se, a = 0.4, b = 0.7, k = 5)
#' }
# Resolve the directory of this script at source() time so that
# .knnimpute_rcpp can find knn_impute.cpp regardless of the caller's
# working directory.
.script_dir <- local({
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  "."
})

#' @rdname imputeBlockKNN
setGeneric("imputeBlockKNN", function(datos, a, b, k, ...) {
  standardGeneric("imputeBlockKNN")
})


# S4 methods --------------------------------------------------------------

#' @rdname imputeBlockKNN
#' @export
setMethod(
  "imputeBlockKNN",
  "matrix",
  function(datos, a, b, k, ...) {
    .imputeBlockKNN_impl(datos, a, b, k)
  }
)

#' @rdname imputeBlockKNN
#' @export
setMethod(
  "imputeBlockKNN",
  "data.frame",
  function(datos, a, b, k, ...) {
    mat <- as.matrix(datos)
    if (!is.numeric(mat)) {
      stop("data.frame must contain only numeric columns")
    }
    result <- .imputeBlockKNN_impl(mat, a, b, k)
    as.data.frame(result)
  }
)

#' @rdname imputeBlockKNN
#' @export
setMethod(
  "imputeBlockKNN",
  "SummarizedExperiment",
  function(datos, a, b, k, ..., assay.type = 1L) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop(
        "SummarizedExperiment package is required. ",
        "Install from Bioconductor."
      )
    }
    mat <- SummarizedExperiment::assay(datos, assay.type)
    imputed <- .imputeBlockKNN_impl(mat, a, b, k)
    SummarizedExperiment::assay(
      datos,
      assay.type,
      withDimnames = FALSE
    ) <- imputed
    datos
  }
)

# Internal workhorse ------------------------------------------------------

#' Core three-block KNN imputation on a numeric matrix
#'
#' @param datos Numeric matrix (variables in rows, samples in columns)
#' @param a Block 1 threshold
#' @param b Block 3 threshold
#' @param k Number of neighbors
#' @return Imputed numeric matrix
#' @keywords internal
.imputeBlockKNN_impl <- function(datos, a, b, k) {
  if (a > b || a < 0 || b > 1) {
    stop(
      "Invalid thresholds: need 0 <= a < b <= 1 (Found a=",
      a,
      ", b=",
      b,
      ")"
    )
  }

  n <- nrow(datos)
  m <- ncol(datos)

  if (m == 1 && any(is.na(datos))) {
    stop("KNN imputation requires at least 2 samples (columns).")
  }

  # Transpose: variables in columns (m x n)
  datos <- t(datos)

  # Track original variable indices
  datosindice <- rbind(1:n, datos)

  # Proportion of NA per variable
  proporcionNA <- colSums(is.na(datos)) / m

  # The blocks are generated (variables with a proportion <=a, variables with
  # a proportion between a and b, and variables with a proportion >=b). In addition, each
  # variable is accompanied by its index:
  bloque1 <- datosindice[, proporcionNA <= a, drop = FALSE]
  bloque2 <- datosindice[, proporcionNA > a & proporcionNA < b, drop = FALSE]
  bloque3 <- datosindice[, proporcionNA >= b, drop = FALSE]

  # Transpose blocks back to variable-in-rows format
  bloque1 <- t(bloque1)
  bloque2 <- t(bloque2)
  bloque3 <- t(bloque3)

  # Extract indices and data
  if (nrow(bloque1) > 0) {
    indice1 <- bloque1[, 1]
    bloque1 <- bloque1[, 2:(m + 1), drop = FALSE]
  } else {
    indice1 <- numeric(0)
    bloque1 <- matrix(nrow = 0, ncol = m)
  }

  if (nrow(bloque2) > 0) {
    indice2 <- bloque2[, 1]
    bloque2 <- bloque2[, 2:(m + 1), drop = FALSE]
  } else {
    indice2 <- numeric(0)
    bloque2 <- matrix(nrow = 0, ncol = m)
  }

  if (nrow(bloque3) > 0) {
    indice3 <- bloque3[, 1]
    bloque3 <- bloque3[, 2:(m + 1), drop = FALSE]
  } else {
    indice3 <- numeric(0)
    bloque3 <- matrix(nrow = 0, ncol = m)
  }

  # Zero-fill blocks 2 and 3
  bloque3mod <- bloque3
  bloque3mod[is.na(bloque3mod)] <- 0

  bloque2mod <- bloque2
  bloque2mod[is.na(bloque2mod)] <- 0

  # For KNN, we use the matrix formed with bloque1, bloque2mod, and
  # bloque3mod, but we estimate the NAs for each individual (column) from
  # bloque1, bloque 2 y bloque3mod:
  matrizKNN <- rbind(bloque1, bloque2mod, bloque3mod)
  individuosKNN <- rbind(bloque1, bloque2, bloque3mod)

  # Impute each individual (column) separately
  datossinNaNKNN <- matrix(0, nrow = n, ncol = m)

  # Working copy of matrizKNN: we swap column i in-place rather than cbind-ing
  # a new matrix each iteration, avoiding m full-matrix allocations.
  matrizWork <- matrizKNN

  for (i in seq_len(m)) {
    if (nrow(matrizWork) == 0) {
      datossinNaNKNN[, i] <- individuosKNN[, i]
      datossinNaNKNN[is.na(datossinNaNKNN[, i]), i] <- 0
      next
    }

    # Swap column i with the version that carries real NAs for this sample
    orig_col <- matrizWork[, i]
    matrizWork[, i] <- individuosKNN[, i]

    if (!any(rowSums(is.na(matrizWork)) == 0L)) {
      stop("All rows contain missing values. Unable to impute.")
    }

    sinNaN <- .knnimpute_rcpp(matrizWork, k)
    datossinNaNKNN[, i] <- sinNaN[, i]

    # Restore column i for the next iteration
    matrizWork[, i] <- orig_col
  }

  # Restore original variable order
  posicion <- c(indice1, indice2, indice3)
  datossinNaNKNNindice <- cbind(posicion, datossinNaNKNN)
  datossinNaNKNNord <- datossinNaNKNNindice[order(posicion), , drop = FALSE]
  datossinNaNKNNord[, 2:(m + 1), drop = FALSE]
}


# knnimpute translation ---------------------------------------------------

#' KNN Imputation Using Rcpp (matching MATLAB's knnimpute)
#'
#' Compiles and calls the C++ implementation in \code{knn_impute.cpp}.
#' Uses globally-complete features for Euclidean distances, includes ties
#' at the k-th boundary, and falls back to 1-NN (nearest valid neighbour)
#' when none of the k nearest have a value for the target feature.
#'
#' @param data Numeric matrix (variables in rows, samples in columns)
#' @param k Number of nearest neighbors
#' @return Imputed matrix
#' @keywords internal
.knnimpute_rcpp <- function(data, k) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp package is required. Install with: install.packages('Rcpp')")
  }

  cpp_file <- file.path(.script_dir, "knn_impute.cpp")

  if (!exists("knn_impute_rcpp", mode = "function", envir = .GlobalEnv)) {
    Rcpp::sourceCpp(cpp_file)
  }

  # C++ expects samples-in-rows (n_samples x n_features)
  X <- t(data)
  result <- knn_impute_rcpp(X, as.integer(k))
  t(result)
}
