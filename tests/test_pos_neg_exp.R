library(testthat)
library(xcms)
library(MsExperiment)
source("../R/pos_neg_exp.R")

# =============================================================================
# Helpers
# =============================================================================

# Load the first n mzML files from a data directory.
# Returns NULL (and the calling test will skip) if no files are found.
load_mse <- function(data_dir, n = 1) {
    # data_dir is relative to project root; tests/ is the working dir
    abs_dir <- file.path("..", data_dir)
    files <- list.files(abs_dir, pattern = "\\.mzML$", full.names = TRUE)
    if (length(files) == 0L) return(NULL)
    readMsExperiment(spectraFiles = head(files, n))
}

skip_no_data <- function() {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
}

# =============================================================================
# perMode()
# =============================================================================

test_that("perMode() creates a PerModeParam object", {
    p <- perMode(pos = "A", neg = "B")
    expect_s3_class(p, "PerModeParam")
    expect_equal(p$pos, "A")
    expect_equal(p$neg, "B")
})

test_that("perMode() neg defaults to pos when omitted", {
    param <- CentWaveParam(peakwidth = c(4, 8))
    p <- perMode(pos = param)
    expect_identical(p$pos, p$neg)
})

test_that("perMode() errors when pos is missing", {
    expect_error(perMode(), "pos must not be missing")
})

test_that("is_per_mode() returns TRUE only for PerModeParam", {
    expect_true(is_per_mode(perMode(pos = "X")))
    expect_false(is_per_mode(list(pos = "X", neg = "X")))
    expect_false(is_per_mode(CentWaveParam()))
    expect_false(is_per_mode(NULL))
})

# =============================================================================
# PosNegMsExp construction (from MsExperiment — no coercion needed)
# =============================================================================

test_that("PosNegMsExp() constructs from two MsExperiment objects", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg), "could not load mzML files")

    pn <- PosNegMsExp(mse_pos, mse_neg)
    expect_s4_class(pn, "PosNegMsExp")
})

test_that("posMode() and negMode() work on PosNegMsExp", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegMsExp(mse_pos, mse_neg)
    expect_s4_class(posMode(pn), "MsExperiment")
    expect_s4_class(negMode(pn), "MsExperiment")
})

test_that("show() on PosNegMsExp runs without error", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegMsExp(mse_pos, mse_neg)
    expect_output(show(pn), "PosNegMsExp")
    expect_output(show(pn), "Positive")
    expect_output(show(pn), "Negative")
})

# =============================================================================
# PosNegXcmsExp construction (from XcmsExperiment)
# =============================================================================

test_that("PosNegXcmsExp() constructs from two XcmsExperiment objects", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg), "could not load mzML files")

    xcms_pos <- as(mse_pos, "XcmsExperiment")
    xcms_neg <- as(mse_neg, "XcmsExperiment")

    pn <- PosNegXcmsExp(xcms_pos, xcms_neg)
    expect_s4_class(pn, "PosNegXcmsExp")
})

test_that("posMode() and negMode() work on PosNegXcmsExp", {
    skip_no_data()
    xcms_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    xcms_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(xcms_pos) || is.null(xcms_neg))

    pn <- PosNegXcmsExp(xcms_pos, xcms_neg)
    expect_identical(posMode(pn), xcms_pos)
    expect_identical(negMode(pn), xcms_neg)
})

test_that("show() on PosNegXcmsExp runs without error", {
    skip_no_data()
    xcms_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    xcms_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(xcms_pos) || is.null(xcms_neg))

    pn <- PosNegXcmsExp(xcms_pos, xcms_neg)
    expect_output(show(pn), "PosNegXcmsExp")
    expect_output(show(pn), "Positive")
    expect_output(show(pn), "Negative")
})

# =============================================================================
# Delegating methods on PosNegMsExp
# =============================================================================

test_that("filterRt() on PosNegMsExp returns a PosNegMsExp", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegMsExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, c(100, 500))
    expect_s4_class(pn_flt, "PosNegMsExp")
})

test_that("filterRt() on PosNegMsExp applies same RT range to both modes", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegMsExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, c(100, 500))

    rt_pos <- range(rtime(spectra(posMode(pn_flt))))
    rt_neg <- range(rtime(spectra(negMode(pn_flt))))
    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 500)
    expect_gte(rt_neg[1], 100); expect_lte(rt_neg[2], 500)
})

# =============================================================================
# Delegating methods on PosNegXcmsExp
# =============================================================================

test_that("filterRt() on PosNegXcmsExp returns a PosNegXcmsExp", {
    skip_no_data()
    xcms_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    xcms_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(xcms_pos) || is.null(xcms_neg))

    pn     <- PosNegXcmsExp(xcms_pos, xcms_neg)
    pn_flt <- filterRt(pn, c(100, 500))
    expect_s4_class(pn_flt, "PosNegXcmsExp")
})

test_that("filterRt() on PosNegXcmsExp applies same RT range to both modes", {
    skip_no_data()
    xcms_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    xcms_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(xcms_pos) || is.null(xcms_neg))

    pn     <- PosNegXcmsExp(xcms_pos, xcms_neg)
    pn_flt <- filterRt(pn, c(100, 500))

    rt_pos <- range(rtime(spectra(posMode(pn_flt))))
    rt_neg <- range(rtime(spectra(negMode(pn_flt))))
    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 500)
    expect_gte(rt_neg[1], 100); expect_lte(rt_neg[2], 500)
})

# =============================================================================
# perMode() routing
# =============================================================================

test_that("perMode() routes distinct params to the correct mode (PosNegMsExp)", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegMsExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, rt = perMode(
        pos = c(100, 600),
        neg = c(200, 700)
    ))

    rt_pos <- range(rtime(spectra(posMode(pn_flt))))
    rt_neg <- range(rtime(spectra(negMode(pn_flt))))

    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 600)
    expect_gte(rt_neg[1], 200); expect_lte(rt_neg[2], 700)
})

test_that("perMode(pos = X) applies X to both modes identically (PosNegMsExp)", {
    skip_no_data()
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegMsExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, rt = perMode(pos = c(100, 500)))

    rt_pos <- range(rtime(spectra(posMode(pn_flt))))
    rt_neg <- range(rtime(spectra(negMode(pn_flt))))

    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 500)
    expect_gte(rt_neg[1], 100); expect_lte(rt_neg[2], 500)
})

# =============================================================================
# Accessor methods — return type must be named list(pos, neg)
# =============================================================================

test_that("hasChromPeaks() on PosNegXcmsExp returns a named list", {
    skip_no_data()
    xcms_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    xcms_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(xcms_pos) || is.null(xcms_neg))

    pn     <- PosNegXcmsExp(xcms_pos, xcms_neg)
    result <- hasChromPeaks(pn)

    expect_type(result, "list")
    expect_named(result, c("pos", "neg"))
    expect_false(result$pos)
    expect_false(result$neg)
})
