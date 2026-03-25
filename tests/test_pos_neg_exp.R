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
# PosNegExp construction
# =============================================================================

test_that("PosNegExp() constructs from two XcmsExperiment objects", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- load_mse("POS_data")
    mse_neg <- load_mse("NEG_data")
    skip_if(is.null(mse_pos) || is.null(mse_neg), "could not load mzML files")

    mse_pos <- as(mse_pos, "XcmsExperiment")
    mse_neg <- as(mse_neg, "XcmsExperiment")

    pn <- PosNegExp(mse_pos, mse_neg)
    expect_s4_class(pn, "PosNegExp")
})

test_that("posExp() and negExp() return the correct slots", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegExp(mse_pos, mse_neg)
    expect_identical(posExp(pn), mse_pos)
    expect_identical(negExp(pn), mse_neg)
})

test_that("show() runs without error", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegExp(mse_pos, mse_neg)
    expect_output(show(pn), "PosNegExp")
    expect_output(show(pn), "Positive")
    expect_output(show(pn), "Negative")
})

# =============================================================================
# Delegating methods — return type must be PosNegExp
# =============================================================================

# filterRt is the cheapest delegating method: no peak detection required.
test_that("filterRt() on PosNegExp returns a PosNegExp", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, c(100, 500))
    expect_s4_class(pn_flt, "PosNegExp")
})

test_that("filterRt() applies same RT range to both modes", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, c(100, 500))

    rt_pos <- range(rtime(spectra(posExp(pn_flt))))
    rt_neg <- range(rtime(spectra(negExp(pn_flt))))
    expect_gte(rt_pos[1], 100)
    expect_lte(rt_pos[2], 500)
    expect_gte(rt_neg[1], 100)
    expect_lte(rt_neg[2], 500)
})

# =============================================================================
# perMode() routing — correct param reaches the correct mode
# =============================================================================

test_that("perMode() routes distinct params to the correct mode", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn <- PosNegExp(mse_pos, mse_neg)

    # Use RT ranges that are intentionally different per mode
    pn_flt <- filterRt(pn, rt = perMode(
        pos = c(100, 600),
        neg = c(200, 700)
    ))

    rt_pos <- range(rtime(spectra(posExp(pn_flt))))
    rt_neg <- range(rtime(spectra(negExp(pn_flt))))

    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 600)
    expect_gte(rt_neg[1], 200); expect_lte(rt_neg[2], 700)
})

test_that("perMode(pos = X) applies X to both modes identically", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegExp(mse_pos, mse_neg)
    pn_flt <- filterRt(pn, rt = perMode(pos = c(100, 500)))

    rt_pos <- range(rtime(spectra(posExp(pn_flt))))
    rt_neg <- range(rtime(spectra(negExp(pn_flt))))

    expect_gte(rt_pos[1], 100); expect_lte(rt_pos[2], 500)
    expect_gte(rt_neg[1], 100); expect_lte(rt_neg[2], 500)
})

# =============================================================================
# Accessor methods — return type must be named list(pos, neg)
# =============================================================================

test_that("hasChromPeaks() returns a named list with pos and neg", {
    skip_if_not(
        file.exists("../POS_data") && file.exists("../NEG_data"),
        "POS_data / NEG_data not available"
    )
    mse_pos <- as(load_mse("POS_data"), "XcmsExperiment")
    mse_neg <- as(load_mse("NEG_data"), "XcmsExperiment")
    skip_if(is.null(mse_pos) || is.null(mse_neg))

    pn     <- PosNegExp(mse_pos, mse_neg)
    result <- hasChromPeaks(pn)

    expect_type(result, "list")
    expect_named(result, c("pos", "neg"))
    expect_false(result$pos)
    expect_false(result$neg)
})
