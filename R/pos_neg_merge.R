# Utilities for combining positive and negative ionization mode annotation
# results into a single SummarizedExperiment.

strip_mode_suffix <- function(x) {
    x <- sub("\\.mzML$", "", x, ignore.case = TRUE)
    sub("_(pos|neg)$", "", x, ignore.case = TRUE)
}

align_neg_to_pos_samples <- function(res_pos, res_neg,
                                     strip_fun = strip_mode_suffix) {
    pos_order <- strip_fun(colnames(res_pos))
    neg_order <- strip_fun(colnames(res_neg))

    neg_reorder_idx <- match(pos_order, neg_order)
    if (anyNA(neg_reorder_idx)) {
        stop("Could not align negative-mode samples to positive-mode samples. ",
             "Missing in NEG: ",
             paste(pos_order[is.na(neg_reorder_idx)], collapse = ", "),
             call. = FALSE)
    }

    res_neg <- res_neg[, neg_reorder_idx]

    if (!identical(pos_order, strip_fun(colnames(res_neg)))) {
        stop("Sample order mismatch remains after reordering negative mode.",
             call. = FALSE)
    }

    list(
        res_neg = res_neg,
        pos_order = pos_order,
        neg_reorder_idx = neg_reorder_idx
    )
}

find_pos_neg_duplicates <- function(res_pos, res_neg,
                                    rt_diff_threshold = 30,
                                    abundance_assay = "ISnorm_filled_imputed",
                                    lipid_col = "target_lipid_name_unique",
                                    rt_col = "rtmed") {
    row_pos <- SummarizedExperiment::rowData(res_pos)
    row_neg <- SummarizedExperiment::rowData(res_neg)

    required_pos <- c(lipid_col, rt_col)
    missing_pos <- setdiff(required_pos, colnames(row_pos))
    missing_neg <- setdiff(required_pos, colnames(row_neg))
    if (length(missing_pos)) {
        stop("Positive result is missing rowData column(s): ",
             paste(missing_pos, collapse = ", "), call. = FALSE)
    }
    if (length(missing_neg)) {
        stop("Negative result is missing rowData column(s): ",
             paste(missing_neg, collapse = ", "), call. = FALSE)
    }
    if (!abundance_assay %in% SummarizedExperiment::assayNames(res_pos)) {
        stop("Positive result is missing assay: ", abundance_assay,
             call. = FALSE)
    }
    if (!abundance_assay %in% SummarizedExperiment::assayNames(res_neg)) {
        stop("Negative result is missing assay: ", abundance_assay,
             call. = FALSE)
    }

    common_lipids <- intersect(row_pos[[lipid_col]], row_neg[[lipid_col]])

    comparison_df <- data.frame(
        lipid_name_unique = common_lipids,
        pos_idx = match(common_lipids, row_pos[[lipid_col]]),
        neg_idx = match(common_lipids, row_neg[[lipid_col]]),
        stringsAsFactors = FALSE
    )

    if (!nrow(comparison_df)) {
        comparison_df$pos_mean_abundance <- numeric(0)
        comparison_df$neg_mean_abundance <- numeric(0)
        comparison_df$pos_rt <- numeric(0)
        comparison_df$neg_rt <- numeric(0)
        comparison_df$rt_diff <- numeric(0)
        comparison_df$rt_similar <- logical(0)
        comparison_df$keep_mode <- character(0)
        return(comparison_df)
    }

    pos_assay <- SummarizedExperiment::assay(res_pos, abundance_assay)
    neg_assay <- SummarizedExperiment::assay(res_neg, abundance_assay)

    comparison_df$pos_mean_abundance <- vapply(
        comparison_df$pos_idx,
        function(i) mean(pos_assay[i, ], na.rm = TRUE),
        numeric(1)
    )
    comparison_df$neg_mean_abundance <- vapply(
        comparison_df$neg_idx,
        function(i) mean(neg_assay[i, ], na.rm = TRUE),
        numeric(1)
    )

    comparison_df$pos_rt <- row_pos[[rt_col]][comparison_df$pos_idx]
    comparison_df$neg_rt <- row_neg[[rt_col]][comparison_df$neg_idx]
    comparison_df$rt_diff <- abs(comparison_df$pos_rt - comparison_df$neg_rt)
    comparison_df$rt_similar <- comparison_df$rt_diff < rt_diff_threshold

    comparison_df$keep_mode <- ifelse(
        comparison_df$pos_mean_abundance >= comparison_df$neg_mean_abundance,
        "positive", "negative"
    )

    comparison_df
}

drop_col_data_columns <- function(res, columns) {
    col_data <- SummarizedExperiment::colData(res)
    keep <- setdiff(colnames(col_data), columns)
    SummarizedExperiment::colData(res) <- col_data[, keep, drop = FALSE]
    res
}

merge_pos_neg_results <- function(res_pos, res_neg,
                                  rt_diff_threshold = 30,
                                  abundance_assay = "ISnorm_filled_imputed",
                                  lipid_col = "target_lipid_name_unique",
                                  rt_col = "rtmed",
                                  strip_fun = strip_mode_suffix) {
    row_pos <- SummarizedExperiment::rowData(res_pos)
    row_neg <- SummarizedExperiment::rowData(res_neg)
    row_pos$ionization_mode <- "positive"
    row_neg$ionization_mode <- "negative"
    SummarizedExperiment::rowData(res_pos) <- row_pos
    SummarizedExperiment::rowData(res_neg) <- row_neg

    aligned <- align_neg_to_pos_samples(res_pos, res_neg, strip_fun = strip_fun)
    res_neg <- aligned$res_neg

    comparison_df <- find_pos_neg_duplicates(
        res_pos = res_pos,
        res_neg = res_neg,
        rt_diff_threshold = rt_diff_threshold,
        abundance_assay = abundance_assay,
        lipid_col = lipid_col,
        rt_col = rt_col
    )

    pos_remove_idx <- comparison_df$pos_idx[
        comparison_df$keep_mode == "negative" & comparison_df$rt_similar
    ]
    neg_remove_idx <- comparison_df$neg_idx[
        comparison_df$keep_mode == "positive" & comparison_df$rt_similar
    ]

    res_pos_filtered <- if (length(pos_remove_idx)) {
        res_pos[-pos_remove_idx, ]
    } else {
        res_pos
    }
    res_neg_filtered <- if (length(neg_remove_idx)) {
        res_neg[-neg_remove_idx, ]
    } else {
        res_neg
    }

    colnames(res_pos_filtered) <- strip_fun(colnames(res_pos_filtered))
    colnames(res_neg_filtered) <- strip_fun(colnames(res_neg_filtered))

    if ("file_name" %in% colnames(SummarizedExperiment::colData(res_pos_filtered))) {
        SummarizedExperiment::colData(res_pos_filtered)$file_name <-
            strip_fun(SummarizedExperiment::colData(res_pos_filtered)$file_name)
    }
    if ("file_name" %in% colnames(SummarizedExperiment::colData(res_neg_filtered))) {
        SummarizedExperiment::colData(res_neg_filtered)$file_name <-
            strip_fun(SummarizedExperiment::colData(res_neg_filtered)$file_name)
    }

    res_pos_filtered <- drop_col_data_columns(
        res_pos_filtered, c("spectraOrigin", "polarity")
    )
    res_neg_filtered <- drop_col_data_columns(
        res_neg_filtered, c("spectraOrigin", "polarity")
    )

    list(
        res_combined = rbind(res_pos_filtered, res_neg_filtered),
        comparison_df = comparison_df,
        pos_remove_idx = pos_remove_idx,
        neg_remove_idx = neg_remove_idx,
        res_pos_filtered = res_pos_filtered,
        res_neg_filtered = res_neg_filtered,
        neg_reorder_idx = aligned$neg_reorder_idx
    )
}
