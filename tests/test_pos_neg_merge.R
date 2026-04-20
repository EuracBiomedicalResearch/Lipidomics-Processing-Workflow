library(testthat)
library(SummarizedExperiment)

source("../R/pos_neg_merge.R")

make_result <- function(mode, lipids, rt, values,
                        sample_ids = c("S1", "S2"),
                        sample_order = seq_along(sample_ids)) {
    mode_suffix <- if (mode == "positive") "pos" else "neg"
    sample_ids <- sample_ids[sample_order]
    sample_names <- paste0(sample_ids, "_", mode_suffix)

    se <- SummarizedExperiment(
        assays = list(ISnorm_filled_imputed = values[, sample_order, drop = FALSE]),
        rowData = DataFrame(
            target_lipid_name_unique = lipids,
            rtmed = rt
        ),
        colData = DataFrame(
            file_name = paste0(sample_names, ".mzML"),
            spectraOrigin = paste0(sample_names, ".mzML"),
            polarity = mode
        )
    )
    colnames(se) <- sample_names
    se
}

test_that("strip_mode_suffix() removes mode and mzML suffixes", {
    expect_equal(
        strip_mode_suffix(c("S1_pos", "S2_neg.mzML", "S3_POS.MZML")),
        c("S1", "S2", "S3")
    )
})

test_that("align_neg_to_pos_samples() reorders negative samples", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(1, 2), nrow = 1)
    )
    res_neg <- make_result(
        "negative", "B", 110,
        matrix(c(3, 4), nrow = 1),
        sample_order = c(2, 1)
    )

    aligned <- align_neg_to_pos_samples(res_pos, res_neg)

    expect_equal(aligned$neg_reorder_idx, c(2, 1))
    expect_equal(strip_mode_suffix(colnames(aligned$res_neg)), c("S1", "S2"))
})

test_that("align_neg_to_pos_samples() fails on missing negative samples", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(1, 2), nrow = 1),
        sample_ids = c("S1", "S2")
    )
    res_neg <- make_result(
        "negative", "B", 110,
        matrix(c(3, 4), nrow = 1),
        sample_ids = c("S1", "S3")
    )

    expect_error(
        align_neg_to_pos_samples(res_pos, res_neg),
        "Could not align negative-mode samples"
    )
})

test_that("merge_pos_neg_results() removes close duplicates with lower abundance", {
    res_pos <- make_result(
        "positive",
        lipids = c("A", "B", "C"),
        rt = c(100, 200, 400),
        values = matrix(c(
            100, 100,
            10, 10,
            5, 5
        ), nrow = 3, byrow = TRUE)
    )
    res_neg <- make_result(
        "negative",
        lipids = c("A", "B", "D"),
        rt = c(110, 205, 500),
        values = matrix(c(
            20, 20,
            100, 100,
            7, 7
        ), nrow = 3, byrow = TRUE)
    )

    result <- merge_pos_neg_results(res_pos, res_neg, rt_diff_threshold = 30)
    combined <- result$res_combined

    expect_equal(result$pos_remove_idx, 2)
    expect_equal(result$neg_remove_idx, 1)
    expect_equal(rowData(combined)$target_lipid_name_unique,
                 c("A", "C", "B", "D"))
    expect_equal(rowData(combined)$ionization_mode,
                 c("positive", "positive", "negative", "negative"))
    expect_equal(dim(combined), c(4, 2))
})

test_that("merge_pos_neg_results() keeps duplicate lipid names when RT differs", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(100, 100), nrow = 1)
    )
    res_neg <- make_result(
        "negative", "A", 180,
        matrix(c(1, 1), nrow = 1)
    )

    result <- merge_pos_neg_results(res_pos, res_neg, rt_diff_threshold = 30)

    expect_length(result$pos_remove_idx, 0)
    expect_length(result$neg_remove_idx, 0)
    expect_equal(rowData(result$res_combined)$target_lipid_name_unique,
                 c("A", "A"))
    expect_equal(rowData(result$res_combined)$ionization_mode,
                 c("positive", "negative"))
})

test_that("merge_pos_neg_results() handles zero common lipids", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(100, 100), nrow = 1)
    )
    res_neg <- make_result(
        "negative", "B", 110,
        matrix(c(1, 1), nrow = 1)
    )

    result <- merge_pos_neg_results(res_pos, res_neg)

    expect_equal(nrow(result$comparison_df), 0)
    expect_length(result$pos_remove_idx, 0)
    expect_length(result$neg_remove_idx, 0)
    expect_equal(dim(result$res_combined), c(2, 2))
})

test_that("merge_pos_neg_results() strips sample suffixes and drops mode-specific colData", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(100, 100), nrow = 1)
    )
    res_neg <- make_result(
        "negative", "B", 110,
        matrix(c(1, 1), nrow = 1),
        sample_order = c(2, 1)
    )

    combined <- merge_pos_neg_results(res_pos, res_neg)$res_combined

    expect_equal(colnames(combined), c("S1", "S2"))
    expect_equal(colData(combined)$file_name, c("S1", "S2"))
    expect_false("spectraOrigin" %in% colnames(colData(combined)))
    expect_false("polarity" %in% colnames(colData(combined)))
})

test_that("merge_pos_neg_results() keeps positive mode on equal abundance ties", {
    res_pos <- make_result(
        "positive", "A", 100,
        matrix(c(10, 10), nrow = 1)
    )
    res_neg <- make_result(
        "negative", "A", 105,
        matrix(c(10, 10), nrow = 1)
    )

    result <- merge_pos_neg_results(res_pos, res_neg)

    expect_length(result$pos_remove_idx, 0)
    expect_equal(result$neg_remove_idx, 1)
    expect_equal(rowData(result$res_combined)$ionization_mode, "positive")
})
