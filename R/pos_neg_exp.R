library(methods)
library(xcms)

# =============================================================================
# Dual-mode experiment wrappers
# =============================================================================
#
# Two S4 classes represent the pipeline at different stages:
#
#   PosNegMsExp   — holds two MsExperiment objects (pre-peak detection)
#   PosNegXcmsExp — holds two XcmsExperiment objects (post-peak detection)
#
# findChromPeaks() on a PosNegMsExp is the "promotion" step that returns a
# PosNegXcmsExp. All subsequent xcms operations work on PosNegXcmsExp.
#
# Usage:
#   data <- PosNegMsExp(mse_pos, mse_neg)
#   data <- filterRt(data, rt = c(10, 800))        # still PosNegMsExp
#   data <- findChromPeaks(data, param = ...)       # → PosNegXcmsExp
#   data <- adjustRtime(data, param = ...)          # still PosNegXcmsExp
#
#   # Mode-specific parameters via perMode():
#   data <- filterRt(data, rt = perMode(
#       pos = c(10, 950),
#       neg = c(10, 800)
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
#' When passed as any named argument to a method dispatched on a
#' \code{PosNegMsExp} or \code{PosNegXcmsExp} object, the positive-mode and
#' negative-mode components are routed to their respective experiment
#' automatically. Must always be passed as a named argument (e.g.
#' \code{rt = perMode(...)}, \code{param = perMode(...)}); positional use
#' is not supported.
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
# split_per_mode() — argument splitting helper
#
# Takes a named list of arguments, finds all PerModeParam values, and returns
# two lists: one with all PerModeParam values replaced by their $pos component,
# one by their $neg component. Non-PerModeParam arguments are kept as-is in
# both lists.
#
# Handles 0, 1, or N perMode() arguments uniformly.
# Errors loudly if a PerModeParam is passed without a name (positional use),
# which would otherwise cause silent corruption of the argument list.
# -----------------------------------------------------------------------------

split_per_mode <- function(args) {
    pm_idx <- which(vapply(args, is_per_mode, logical(1L)))
    if (length(pm_idx) == 0L) return(list(pos = args, neg = args))
    pos_args <- args
    neg_args <- args
    for (i in pm_idx) {
        nm <- names(args)[i]
        if (is.null(nm) || !nzchar(nm))
            stop("perMode() must be passed as a named argument, ",
                 "e.g. rt = perMode(...) or param = perMode(...)")
        pos_args[[nm]] <- args[[i]]$pos
        neg_args[[nm]] <- args[[i]]$neg
    }
    list(pos = pos_args, neg = neg_args)
}

# -----------------------------------------------------------------------------
# S4 classes
# -----------------------------------------------------------------------------

#' Dual-mode experiment wrapping two MsExperiment objects (pre-peak detection).
setClass("PosNegMsExp",
    slots = c(
        pos = "MsExperiment",
        neg = "MsExperiment"
    )
)

#' Dual-mode experiment wrapping two XcmsExperiment objects (post-peak
#' detection).
setClass("PosNegXcmsExp",
    slots = c(
        pos = "XcmsExperiment",
        neg = "XcmsExperiment"
    )
)

#' Create a pre-peak-detection dual-mode experiment.
#'
#' @param pos An \code{MsExperiment} for positive ionization mode.
#' @param neg An \code{MsExperiment} for negative ionization mode.
#' @return A \code{PosNegMsExp} object.
PosNegMsExp <- function(pos, neg) new("PosNegMsExp", pos = pos, neg = neg)

#' Create a post-peak-detection dual-mode experiment.
#'
#' @param pos An \code{XcmsExperiment} for positive ionization mode.
#' @param neg An \code{XcmsExperiment} for negative ionization mode.
#' @return A \code{PosNegXcmsExp} object.
PosNegXcmsExp <- function(pos, neg) new("PosNegXcmsExp", pos = pos, neg = neg)

# Accessors — same generics registered for both classes.
setGeneric("posExp", function(x) standardGeneric("posExp"))
setMethod("posExp", "PosNegMsExp",   function(x) x@pos)
setMethod("posExp", "PosNegXcmsExp", function(x) x@pos)

setGeneric("negExp", function(x) standardGeneric("negExp"))
setMethod("negExp", "PosNegMsExp",   function(x) x@neg)
setMethod("negExp", "PosNegXcmsExp", function(x) x@neg)

# =============================================================================
# PosNegMsExp methods  (pre-peak detection)
# =============================================================================
#
# All generics below were verified to use "object" as the first formal in
# xcms 4.8.0.
#
# Filter generics (filterRt, filterFile, filterMzRange) have signature
# function(object, ...) — no explicit named formals beyond object — so
# all user arguments land in ... and split_per_mode() handles them correctly.
#
# findChromPeaks has signature function(object, param, ...) — param is an
# explicit formal. S4 binds the caller's param= value directly to param
# inside the method, so we handle it separately before splitting the rest
# of ... for any additional perMode() arguments.

setMethod("filterRt", "PosNegMsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegMsExp(
        do.call("filterRt", c(list(object@pos), s$pos)),
        do.call("filterRt", c(list(object@neg), s$neg))
    )
})

setMethod("filterFile", "PosNegMsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegMsExp(
        do.call("filterFile", c(list(object@pos), s$pos)),
        do.call("filterFile", c(list(object@neg), s$neg))
    )
})

setMethod("filterMzRange", "PosNegMsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegMsExp(
        do.call("filterMzRange", c(list(object@pos), s$pos)),
        do.call("filterMzRange", c(list(object@neg), s$neg))
    )
})

# Promotion step: PosNegMsExp → PosNegXcmsExp.
# param may itself be a perMode() object (different CentWaveParam per mode).
setMethod("findChromPeaks", "PosNegMsExp", function(object, param, ...) {
    s <- split_per_mode(c(list(param = param), list(...)))
    PosNegXcmsExp(
        do.call("findChromPeaks", c(list(object@pos), s$pos)),
        do.call("findChromPeaks", c(list(object@neg), s$neg))
    )
})

# =============================================================================
# PosNegXcmsExp methods  (post-peak detection)
# =============================================================================
#
# Processing generics keep the object as PosNegXcmsExp.
# Filter generics are repeated here so post-peak filtering still works.
#
# Generics with explicit param formal (verified in xcms 4.8.0):
#   refineChromPeaks, adjustRtime, groupChromPeaks, fillChromPeaks
#
# Generics with only ... beyond object:
#   filterRt, filterFile, filterMzRange, dropAdjustedRtime
#
# Note: applyAdjustedRtime is NOT an S4 generic in xcms 4.8.0 — it is a
# plain function. It cannot be overloaded with setMethod(). Call it on each
# mode separately:
#   data@pos <- applyAdjustedRtime(posExp(data))
#   data@neg <- applyAdjustedRtime(negExp(data))
# or use dropAdjustedRtime / adjustRtime which are proper generics.

setMethod("filterRt", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegXcmsExp(
        do.call("filterRt", c(list(object@pos), s$pos)),
        do.call("filterRt", c(list(object@neg), s$neg))
    )
})

setMethod("filterFile", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegXcmsExp(
        do.call("filterFile", c(list(object@pos), s$pos)),
        do.call("filterFile", c(list(object@neg), s$neg))
    )
})

setMethod("filterMzRange", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegXcmsExp(
        do.call("filterMzRange", c(list(object@pos), s$pos)),
        do.call("filterMzRange", c(list(object@neg), s$neg))
    )
})

setMethod("dropAdjustedRtime", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    PosNegXcmsExp(
        do.call("dropAdjustedRtime", c(list(object@pos), s$pos)),
        do.call("dropAdjustedRtime", c(list(object@neg), s$neg))
    )
})

setMethod("refineChromPeaks", "PosNegXcmsExp", function(object, param, ...) {
    s <- split_per_mode(c(list(param = param), list(...)))
    PosNegXcmsExp(
        do.call("refineChromPeaks", c(list(object@pos), s$pos)),
        do.call("refineChromPeaks", c(list(object@neg), s$neg))
    )
})

setMethod("adjustRtime", "PosNegXcmsExp", function(object, param, ...) {
    s <- split_per_mode(c(list(param = param), list(...)))
    PosNegXcmsExp(
        do.call("adjustRtime", c(list(object@pos), s$pos)),
        do.call("adjustRtime", c(list(object@neg), s$neg))
    )
})

setMethod("groupChromPeaks", "PosNegXcmsExp", function(object, param, ...) {
    s <- split_per_mode(c(list(param = param), list(...)))
    PosNegXcmsExp(
        do.call("groupChromPeaks", c(list(object@pos), s$pos)),
        do.call("groupChromPeaks", c(list(object@neg), s$neg))
    )
})

setMethod("fillChromPeaks", "PosNegXcmsExp", function(object, param, ...) {
    s <- split_per_mode(c(list(param = param), list(...)))
    PosNegXcmsExp(
        do.call("fillChromPeaks", c(list(object@pos), s$pos)),
        do.call("fillChromPeaks", c(list(object@neg), s$neg))
    )
})

# =============================================================================
# Accessor methods on PosNegXcmsExp
# =============================================================================
#
# These return data (not an experiment object), wrapped as list(pos=…, neg=…).
# All have signature function(object, ...) in xcms 4.8.0.
# perMode() splitting is applied so callers can filter accessors by mode
# if the underlying generic supports it.
#
# quantify() returns a SummarizedExperiment per mode — belongs here, not above.

setMethod("chromPeaks", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("chromPeaks", c(list(object@pos), s$pos)),
         neg = do.call("chromPeaks", c(list(object@neg), s$neg)))
})

setMethod("chromPeakData", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("chromPeakData", c(list(object@pos), s$pos)),
         neg = do.call("chromPeakData", c(list(object@neg), s$neg)))
})

setMethod("featureDefinitions", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("featureDefinitions", c(list(object@pos), s$pos)),
         neg = do.call("featureDefinitions", c(list(object@neg), s$neg)))
})

setMethod("featureValues", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("featureValues", c(list(object@pos), s$pos)),
         neg = do.call("featureValues", c(list(object@neg), s$neg)))
})

setMethod("adjustedRtime", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("adjustedRtime", c(list(object@pos), s$pos)),
         neg = do.call("adjustedRtime", c(list(object@neg), s$neg)))
})

setMethod("hasAdjustedRtime", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("hasAdjustedRtime", c(list(object@pos), s$pos)),
         neg = do.call("hasAdjustedRtime", c(list(object@neg), s$neg)))
})

setMethod("hasChromPeaks", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("hasChromPeaks", c(list(object@pos), s$pos)),
         neg = do.call("hasChromPeaks", c(list(object@neg), s$neg)))
})

setMethod("hasFeatures", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("hasFeatures", c(list(object@pos), s$pos)),
         neg = do.call("hasFeatures", c(list(object@neg), s$neg)))
})

setMethod("quantify", "PosNegXcmsExp", function(object, ...) {
    s <- split_per_mode(list(...))
    list(pos = do.call("quantify", c(list(object@pos), s$pos)),
         neg = do.call("quantify", c(list(object@neg), s$neg)))
})

# -----------------------------------------------------------------------------
# show methods
# -----------------------------------------------------------------------------

setMethod("show", "PosNegMsExp", function(object) {
    cat("PosNegMsExp — dual ionization mode experiment (pre-peak detection)\n")
    cat("  Positive:", length(object@pos), "sample(s)\n")
    cat("  Negative:", length(object@neg), "sample(s)\n")
})

setMethod("show", "PosNegXcmsExp", function(object) {
    cat("PosNegXcmsExp — dual ionization mode experiment (post-peak detection)\n")
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

# =============================================================================
# PosNegSumExp — dual-mode quantitative result container
# =============================================================================
#
# Wraps two SummarizedExperiment objects (one per ionization mode).
# Functions that accept a SummarizedExperiment can accept a PosNegSumExp
# via inherits() guards, applying the operation to both modes and
# returning an updated PosNegSumExp.
#
# Usage:
#   res_both <- PosNegSumExp(res_pos, res_neg)
#   res_both <- apply_volume_correction(res_both, VOLUME_FACTORS, ...)
#   res_both <- normalize_by_is(res_both, ...)
#   res_both <- filter_by_qc_rsd(res_both, threshold = 0.3)

setClass("PosNegSumExp",
    slots = c(pos = "SummarizedExperiment", neg = "SummarizedExperiment"))

#' Create a dual-mode quantitative result container.
#'
#' @param pos A \code{SummarizedExperiment} for positive ionization mode.
#' @param neg A \code{SummarizedExperiment} for negative ionization mode.
#' @return A \code{PosNegSumExp} object.
PosNegSumExp <- function(pos, neg)
    new("PosNegSumExp", pos = pos, neg = neg)

# Accessors
setGeneric("posRes", function(x) standardGeneric("posRes"))
setGeneric("negRes", function(x) standardGeneric("negRes"))

setMethod("posRes", "PosNegSumExp", function(x) x@pos)
setMethod("negRes", "PosNegSumExp", function(x) x@neg)

setGeneric("posRes<-", function(x, value) standardGeneric("posRes<-"))
setGeneric("negRes<-", function(x, value) standardGeneric("negRes<-"))

setReplaceMethod("posRes", "PosNegSumExp", function(x, value) {
    x@pos <- value
    x
})
setReplaceMethod("negRes", "PosNegSumExp", function(x, value) {
    x@neg <- value
    x
})

setMethod("show", "PosNegSumExp", function(object) {
    cat("PosNegSumExp\n")
    cat("  Positive:", nrow(object@pos), "features,",
        ncol(object@pos), "samples\n")
    cat("  Negative:", nrow(object@neg), "features,",
        ncol(object@neg), "samples\n")
})
