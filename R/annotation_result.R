library(methods)

# =============================================================================
# AnnotationResult — annotation state container
# =============================================================================
#
# Bundles the three objects that always travel together through the annotation
# pipeline after rank-1 matching:
#
#   matches  — evolving matched data frame (mtched_data)
#   query    — experimental features extracted from the SummarizedExperiment
#   database — full lipid database for this ionization mode (all ranks)
#
# The database is stored in full (not rank-1 only) so that match_adducts()
# can filter for rank > 1 entries internally without needing an external arg.
#
# Usage:
#   annotation <- match_features_to_database(res, lipid_database, ...)
#   # returns AnnotationResult
#
#   annotation <- match_adducts(annotation, ppm = 20, rt_tol = 5)
#   annotation <- resolve_annotation_ambiguities(annotation)
#   annotation <- resolve_sm_isomers(annotation)

#' @slot matches  data.frame. Matched lipid-feature pairs (evolves at each step).
#' @slot query    data.frame. Experimental feature table from the SE rowData.
#' @slot database data.frame. Full lipid database (all ranks) for this mode.
setClass("AnnotationResult",
    slots = c(
        matches  = "data.frame",
        query    = "data.frame",
        database = "data.frame"
    )
)

#' Create an AnnotationResult object
#'
#' @param matches  data.frame (or tibble) of matched lipid-feature pairs.
#' @param query    data.frame (or tibble) of experimental features.
#' @param database data.frame (or tibble) of the full lipid database (all ranks).
#' @return An \code{AnnotationResult} object.
AnnotationResult <- function(matches, query, database)
    new("AnnotationResult", matches = matches, query = query, database = database)

# -----------------------------------------------------------------------------
# Accessors
# -----------------------------------------------------------------------------

setGeneric("annotMatches",  function(x) standardGeneric("annotMatches"))
setGeneric("annotQuery",    function(x) standardGeneric("annotQuery"))
setGeneric("annotDatabase", function(x) standardGeneric("annotDatabase"))

setMethod("annotMatches",  "AnnotationResult", function(x) x@matches)
setMethod("annotQuery",    "AnnotationResult", function(x) x@query)
setMethod("annotDatabase", "AnnotationResult", function(x) x@database)

setGeneric("annotMatches<-",
    function(x, value) standardGeneric("annotMatches<-"))
setReplaceMethod("annotMatches", "AnnotationResult", function(x, value) {
    x@matches <- value
    x
})

# -----------------------------------------------------------------------------
# show
# -----------------------------------------------------------------------------

setMethod("show", "AnnotationResult", function(object) {
    cat("AnnotationResult\n")
    cat("  Matches: ", nrow(object@matches),  " rows\n",   sep = "")
    cat("  Features:", nrow(object@query),    "\n",         sep = "")
    cat("  Database:", nrow(object@database), " entries\n", sep = "")
})

# =============================================================================
# PosNegAnnotation — dual-mode annotation state container
# =============================================================================
#
# Wraps one AnnotationResult per ionization mode. Functions that accept an
# AnnotationResult automatically accept a PosNegAnnotation via the inherits()
# guards in lipid_helpers.R, applying the operation to both modes and
# returning an updated PosNegAnnotation.
#
# Usage:
#   annotation <- PosNegAnnotation(pos = annot_pos, neg = annot_neg)
#   annotation <- match_adducts(annotation, ppm = 20, rt_tol = 5)
#   annotation <- resolve_annotation_ambiguities(annotation)
#   annotation <- resolve_sm_isomers(annotation)

setClass("PosNegAnnotation",
    slots = c(pos = "AnnotationResult", neg = "AnnotationResult"))

#' Create a dual-mode annotation container
#'
#' @param pos An \code{AnnotationResult} for positive ionization mode.
#' @param neg An \code{AnnotationResult} for negative ionization mode.
#' @return A \code{PosNegAnnotation} object.
PosNegAnnotation <- function(pos, neg)
    new("PosNegAnnotation", pos = pos, neg = neg)

# -----------------------------------------------------------------------------
# Accessors
# -----------------------------------------------------------------------------

setGeneric("posAnnot", function(x) standardGeneric("posAnnot"))
setGeneric("negAnnot", function(x) standardGeneric("negAnnot"))

setMethod("posAnnot", "PosNegAnnotation", function(x) x@pos)
setMethod("negAnnot", "PosNegAnnotation", function(x) x@neg)

setGeneric("posAnnot<-", function(x, value) standardGeneric("posAnnot<-"))
setGeneric("negAnnot<-", function(x, value) standardGeneric("negAnnot<-"))

setReplaceMethod("posAnnot", "PosNegAnnotation", function(x, value) {
    x@pos <- value
    x
})
setReplaceMethod("negAnnot", "PosNegAnnotation", function(x, value) {
    x@neg <- value
    x
})

# -----------------------------------------------------------------------------
# show
# -----------------------------------------------------------------------------

setMethod("show", "PosNegAnnotation", function(object) {
    cat("PosNegAnnotation\n")
    cat("  Positive:", nrow(object@pos@matches), "matches,",
        nrow(object@pos@database), "db entries\n")
    cat("  Negative:", nrow(object@neg@matches), "matches,",
        nrow(object@neg@database), "db entries\n")
})
