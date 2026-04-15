# =============================================================================
# Regression snapshot helpers
# =============================================================================
#
# save_snapshot() serialises a lightweight R object to an RDS file so that
# two runs of the same workflow can be compared by test_regression*.R.
#
# File naming convention:
#   {step_name}_{polarity}_{smoke|full}.rds
#
# Typical usage:
#   save_snapshot(mtched_data[, snap_cols], "annot_rank1", "pos")
#   save_snapshot(mtched_data[, snap_cols], "annot_rank1", "neg")
#
# The caller must have SNAPSHOT_DIR defined (set in the workflow config section)
# before calling save_snapshot().

#' Save a lightweight regression snapshot
#'
#' Serialises \code{object} to an RDS file under \code{snap_dir}. The filename
#' encodes \code{step_name}, \code{polarity}, and a \code{"smoke"}/\code{"full"}
#' tag so that different run modes can be compared independently.
#'
#' Basic sanity checks are applied before saving and any issues are reported as
#' warnings (the file is still written so the run is not blocked):
#' \itemize{
#'   \item \code{data.frame}: \code{nrow >= min_rows}; warns on columns >80\% NA.
#'   \item Structured \code{list}: checks \code{r_squared >= min_r2},
#'     \code{n_features >= 1}, and flags all-NA numeric elements.
#'   \item List of numeric vectors (e.g. \code{adjustedRtime}): checks for
#'     zero-length or all-NA elements.
#' }
#'
#' @param object R object to serialise. Keep it small (data frames, vectors,
#'   compact named lists).
#' @param step_name Short identifier used in the filename,
#'   e.g. \code{"annot_rank1"}.
#' @param polarity \code{"pos"} or \code{"neg"}.
#' @param smoke Logical. \code{TRUE} tags the file as \code{"smoke"};
#'   \code{FALSE} (default) tags it as \code{"full"}.
#' @param min_rows For data frames: minimum acceptable row count (default 1).
#' @param min_r2 For structured lists with an \code{r_squared} element:
#'   minimum acceptable R\u00b2 (default 0.5).
#' @param snap_dir Output directory. Defaults to \code{SNAPSHOT_DIR} looked up
#'   in the calling environment.
#' @return The file path, invisibly.
save_snapshot <- function(object, step_name, polarity, smoke = FALSE,
                          min_rows = 1L, min_r2 = 0.5,
                          snap_dir = SNAPSHOT_DIR) {
    if (!dir.exists(snap_dir)) dir.create(snap_dir, recursive = TRUE)
    mode_tag <- if (isTRUE(smoke)) "smoke" else "full"
    fname <- file.path(snap_dir,
                       paste0(step_name, "_", polarity, "_", mode_tag, ".rds"))

    issues <- character(0)

    if (is.data.frame(object)) {
        if (nrow(object) < min_rows)
            issues <- c(issues, sprintf("data frame has %d rows (min %d)",
                                        nrow(object), min_rows))
        na_frac  <- colMeans(is.na(object))
        bad_cols <- names(na_frac[na_frac > 0.8])
        if (length(bad_cols))
            issues <- c(issues, sprintf("column(s) >80%% NA: %s",
                                        paste(bad_cols, collapse = ", ")))

    } else if (is.list(object) && !is.data.frame(object)) {
        if (all(vapply(object, is.numeric, logical(1)))) {
            # List of per-sample numeric vectors (e.g. adjustedRtime)
            empty_els <- sum(vapply(object, function(x) length(x) == 0L, logical(1)))
            if (empty_els > 0)
                issues <- c(issues, sprintf("%d element(s) have length 0",
                                            empty_els))
            allna_els <- sum(vapply(object, function(x) all(is.na(x)), logical(1)))
            if (allna_els > 0)
                issues <- c(issues, sprintf("%d element(s) are all-NA",
                                            allna_els))
        } else {
            # Structured list (RT correction, gap-fill stats, etc.)
            if (!is.null(object$r_squared) && object$r_squared < min_r2)
                issues <- c(issues, sprintf("R\u00b2 = %.3f < %.1f",
                                            object$r_squared, min_r2))
            if (!is.null(object$n_features) && object$n_features < 1L)
                issues <- c(issues, "n_features = 0")
            na_nums <- Filter(function(x) is.numeric(x) && all(is.na(x)), object)
            if (length(na_nums))
                issues <- c(issues, sprintf("all-NA element(s): %s",
                                            paste(names(na_nums), collapse = ", ")))
        }
    }

    if (length(issues))
        warning("Snapshot '", step_name, "' (", mode_tag, "/", polarity, "): ",
                paste(issues, collapse = "; "),
                " \u2014 snapshot saved but may be degenerate for regression use")

    saveRDS(object, fname)
    message(sprintf("Snapshot saved [%s/%s]: %s  (%s)",
                    mode_tag, polarity, basename(fname),
                    if (length(issues)) "WARNINGS \u2014 see above" else "OK"))
    invisible(fname)
}
