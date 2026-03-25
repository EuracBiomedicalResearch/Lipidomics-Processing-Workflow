library(methods)
library(xcms)

# =============================================================================
# PosNegExp — dual-mode XcmsExperiment wrapper
# =============================================================================
#
# Holds one XcmsExperiment per ionization mode and delegates xcms method calls
# to both, optionally with mode-specific parameters via perMode().
#
# Usage:
#   data <- PosNegExp(mse_pos, mse_neg)
#
#   # Same parameters for both modes
#   data <- findChromPeaks(data, param = CentWaveParam(peakwidth = c(4, 8)))
#
#   # Different parameters per mode
#   data <- findChromPeaks(data, param = perMode(
#       pos = CentWaveParam(peakwidth = c(4, 8)),
#       neg = CentWaveParam(peakwidth = c(3, 7))
#   ))

# -----------------------------------------------------------------------------
# perMode() — per-mode parameter container (S3)
#
# Intentionally S3: it is a lightweight tagged container, not a domain object
# requiring S4 inheritance or slot validation. Mixing S3 containers with S4
# dispatch is idiomatic in Bioconductor.
# -----------------------------------------------------------------------------

#' Specify different parameters for positive and negative ionization modes.
#'
#' When passed as the \code{param} argument to a method dispatched on a
#' \code{PosNegExp} object, the positive-mode and negative-mode components are
#' routed to their respective \code{XcmsExperiment} automatically.
#'
#' \code{param} is used as the interception point because xcms processing
#' functions consistently use that argument name for their \code{*Param}
#' objects (e.g. \code{CentWaveParam}, \code{PeakDensityParam}).
#'
#' If only \code{pos} is supplied, \code{neg} defaults to the same value so
#' symmetric cases require no repetition.
#'
#' @param pos Parameter object for positive mode.
#' @param neg Parameter object for negative mode. Defaults to \code{pos}.
#'
#' @return An S3 object of class \code{PerModeParam} (named list with class
#'   tag).
#'
#' @examples
#' perMode(
#'     pos = CentWaveParam(peakwidth = c(4, 8)),
#'     neg = CentWaveParam(peakwidth = c(3, 7))
#' )
perMode <- function(pos, neg = pos) {
    stopifnot("pos must not be missing" = !missing(pos))
    structure(list(pos = pos, neg = neg), class = "PerModeParam")
}

is_per_mode <- function(x) inherits(x, "PerModeParam")

# -----------------------------------------------------------------------------
# PosNegExp S4 class
# -----------------------------------------------------------------------------

setClass("PosNegExp",
    slots = c(
        pos = "XcmsExperiment",
        neg = "XcmsExperiment"
    )
)

#' Create a dual-mode experiment object.
#'
#' @param pos An \code{XcmsExperiment} for positive ionization mode.
#' @param neg An \code{XcmsExperiment} for negative ionization mode.
#'
#' @return A \code{PosNegExp} object.
PosNegExp <- function(pos, neg) {
    new("PosNegExp", pos = pos, neg = neg)
}

setGeneric("posExp", function(x) standardGeneric("posExp"))
setMethod("posExp", "PosNegExp", function(x) x@pos)

setGeneric("negExp", function(x) standardGeneric("negExp"))
setMethod("negExp", "PosNegExp", function(x) x@neg)

# -----------------------------------------------------------------------------
# Delegating methods — return type is XcmsExperiment, so results are wrapped
# back into a PosNegExp to keep the pipeline flowing.
# -----------------------------------------------------------------------------

# Factory returning a method body that splits PerModeParam and delegates to
# each mode, or applies the same argument to both when params are identical.
#
# force(fname) is required to correctly capture fname in the closure when
# called from a loop — without it all methods would reference the last value
# of the loop variable at call time (classic R late-binding bug).
makeDelegatingMethod <- function(fname) {
    force(fname)
    function(object, ...) {
        args <- list(...)

        # Find whichever argument (by any name) is a PerModeParam.
        # This covers both the common xcms convention (param = perMode(...))
        # and functions that use other argument names (e.g. rt = perMode(...)).
        pm_idx <- which(vapply(args, is_per_mode, logical(1L)))

        if (length(pm_idx) == 1L) {
            arg_name          <- names(args)[pm_idx]
            pos_args          <- args
            neg_args          <- args
            pos_args[[arg_name]] <- args[[pm_idx]]$pos
            neg_args[[arg_name]] <- args[[pm_idx]]$neg
        } else {
            pos_args <- args
            neg_args <- args
        }

        pos_result <- do.call(fname, c(list(object@pos), pos_args))
        neg_result <- do.call(fname, c(list(object@neg), neg_args))

        PosNegExp(pos_result, neg_result)
    }
}

# xcms generics whose return type is XcmsExperiment; extend as needed.
.XCMS_DELEGATING <- c(
    "findChromPeaks",
    "refineChromPeaks",
    "adjustRtime",
    "applyAdjustedRtime",
    "dropAdjustedRtime",
    "groupChromPeaks",
    "fillChromPeaks",
    "filterRt",
    "filterMzRange",
    "filterFile"
)

for (.gen in .XCMS_DELEGATING) {
    if (isGeneric(.gen))
        setMethod(.gen, "PosNegExp", makeDelegatingMethod(.gen))
}
rm(.gen)

# -----------------------------------------------------------------------------
# Accessor methods — return type is data (matrix, DataFrame, logical, …), not
# XcmsExperiment. Results are returned as a named list(pos = …, neg = …) so
# callers can do chromPeaks(data)$pos or iterate over both modes.
#
# Note: quantify() returns a SummarizedExperiment, so it belongs here rather
# than in the delegating group above.
# -----------------------------------------------------------------------------

.makeAccessorMethod <- function(fname) {
    force(fname)
    function(object, ...) {
        list(
            pos = do.call(fname, c(list(object@pos), list(...))),
            neg = do.call(fname, c(list(object@neg), list(...)))
        )
    }
}

.XCMS_ACCESSORS <- c(
    "chromPeaks",
    "chromPeakData",
    "featureDefinitions",
    "featureValues",
    "adjustedRtime",
    "hasAdjustedRtime",
    "hasChromPeaks",
    "hasFeatures",
    "quantify"
)

for (.acc in .XCMS_ACCESSORS) {
    if (isGeneric(.acc))
        setMethod(.acc, "PosNegExp", .makeAccessorMethod(.acc))
}
rm(.acc)

# -----------------------------------------------------------------------------
# show method
# -----------------------------------------------------------------------------

setMethod("show", "PosNegExp", function(object) {
    cat("PosNegExp — dual ionization mode experiment\n")
    cat("  Positive:", length(object@pos), "sample(s)")
    if (hasChromPeaks(object@pos))
        cat(",", nrow(chromPeaks(object@pos)), "peaks")
    if (hasFeatures(object@pos))
        cat(",", nrow(featureDefinitions(object@pos)), "features")
    cat("\n")
    cat("  Negative:", length(object@neg), "sample(s)")
    if (hasChromPeaks(object@neg))
        cat(",", nrow(chromPeaks(object@neg)), "peaks")
    if (hasFeatures(object@neg))
        cat(",", nrow(featureDefinitions(object@neg)), "features")
    cat("\n")
})
