#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
study_id <- if (length(args) >= 1) args[[1]] else "posNegExp"

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("applications/pilot_study/generate_annotation_metrics.R")
}
study_dir <- dirname(script_path)
project_root <- normalizePath(file.path(study_dir, "../.."))

source(file.path(project_root, "R/lipid_helpers.R"))
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(MsIO)
  library(alabaster.se)
  library(SummarizedExperiment)
  library(xcms)
  library(Spectra)
  library(MetaboCoreUtils)
  library(MetaboAnnotation)
  library(enviPat)
  library(ProtGenerics)
})

MATCH_PPM <- 20
MATCH_RT_TOL <- 20
ISOPEAK_SIM_THRESHOLD <- 0.78
PPM <- 10

first_present <- function(x, candidates) {
  candidates[candidates %in% colnames(x)][1]
}

bind_rows_aligned <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  rows <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (!length(rows)) {
    return(data.frame())
  }

  all_cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing_cols <- setdiff(all_cols, colnames(x))
    for (col in missing_cols) x[[col]] <- NA
    x[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

expand_annotation_source <- function(x) {
  if (is.null(x) || !nrow(x)) {
    return(data.frame())
  }

  lipid_col <- first_present(
    x,
    c("target_lipid_name_unique", "lipid_name_unique",
      "target_lipid.name", "lipid.name", "target_lipid_name", "lipid_name")
  )
  if (is.na(lipid_col) || !"feature_id" %in% colnames(x)) {
    return(data.frame())
  }

  expanded <- vector("list", nrow(x))
  for (i in seq_len(nrow(x))) {
    lipids <- as.character(x[[lipid_col]][i])
    feature <- as.character(x$feature_id[i])
    if (is.na(lipids) || !nzchar(lipids) || is.na(feature) || !nzchar(feature)) {
      next
    }
    lipids <- unlist(strsplit(lipids, "\\s*;\\s*"), use.names = FALSE)
    lipids <- lipids[nzchar(lipids)]
    if (!length(lipids)) next
    expanded[[i]] <- data.frame(
      source_row = i,
      feature_id = feature,
      lipid_name = lipids,
      stringsAsFactors = FALSE
    )
  }

  expanded <- do.call(rbind, Filter(Negate(is.null), expanded))
  if (is.null(expanded)) {
    return(data.frame())
  }
  expanded$pair_key <- paste(expanded$feature_id, expanded$lipid_name, sep = "||")
  expanded
}

pair_metadata <- function(x) {
  expanded <- expand_annotation_source(x)
  if (!nrow(expanded)) {
    return(data.frame())
  }

  keys <- unique(expanded$pair_key)
  first_rows <- expanded$source_row[match(keys, expanded$pair_key)]
  pair_counts <- as.integer(table(factor(expanded$pair_key, levels = keys)))

  metadata_cols <- intersect(
    c("target_lipid.name", "target_lipid_name_unique",
      "target_LIPID.CATEGORY..ABBREV.", "target_LIPID.SUBCLASS..ABBREV.",
      "target_mz", "mzmed", "target_rt_adjusted", "target_rt_sd",
      "rtmed", "ppm_error", "score_rt", "target_rank", "score",
      "isopeak_count", "isopeak_sim", "adduct_count", "adduct_ratio",
      "nb_annotations", "target_IS_norm"),
    colnames(x)
  )

  out <- expanded[match(keys, expanded$pair_key),
                  c("feature_id", "lipid_name", "pair_key"), drop = FALSE]
  out$source_row_count <- pair_counts

  if (length(metadata_cols)) {
    out <- cbind(out, x[first_rows, metadata_cols, drop = FALSE])
  }

  adduct_col <- first_present(x, c("target_adduct", "target_Adduct"))
  if (!is.na(adduct_col)) {
    adducts <- tapply(
      as.character(x[[adduct_col]][expanded$source_row]),
      expanded$pair_key,
      function(z) paste(unique(z[!is.na(z) & nzchar(z)]), collapse = "; ")
    )
    out$adducts <- unname(adducts[out$pair_key])
  }

  rownames(out) <- NULL
  out
}

make_phase_detail <- function(candidate_annotations, truth, phase, polarity) {
  candidates <- pair_metadata(candidate_annotations)
  truth_meta <- pair_metadata(truth)

  truth_keys <- unique(truth_meta$pair_key)
  truth_features <- unique(truth_meta$feature_id)
  candidate_keys <- unique(candidates$pair_key)

  if (nrow(candidates)) {
    candidates$annotation_status <- ifelse(
      candidates$pair_key %in% truth_keys,
      "true",
      ifelse(candidates$feature_id %in% truth_features, "false", "unknown")
    )
    candidates$interpretation <- ifelse(
      candidates$annotation_status == "true",
      "Recovered the manually curated annotation.",
      ifelse(candidates$annotation_status == "false",
             "Candidate for a curated feature, but a different annotation was selected during manual curation.",
             "Candidate for a feature absent from the curated truth set; not counted as wrong.")
    )
    candidates$candidate_generated <- "yes"
    candidates$feature_scope <- ifelse(
      candidates$feature_id %in% truth_features,
      "curated_feature",
      "uncurated_feature"
    )
    candidates$included_in_curated_metrics <- ifelse(
      candidates$feature_scope == "curated_feature", "yes", "no"
    )
  }

  missed <- truth_meta[!truth_meta$pair_key %in% candidate_keys, , drop = FALSE]
  if (nrow(missed)) {
    missed$annotation_status <- "missed"
    missed$interpretation <- "Manually curated annotation not recovered in this phase."
    missed$candidate_generated <- "no"
    missed$feature_scope <- "curated_feature"
    missed$included_in_curated_metrics <- "yes"
  }

  out <- bind_rows_aligned(list(candidates, missed))
  if (!nrow(out)) {
    out <- data.frame(feature_id = character(), lipid_name = character(),
                      pair_key = character(), stringsAsFactors = FALSE)
  }

  out$polarity <- polarity
  out$phase <- phase

  rename_map <- c(
    target_LIPID.CATEGORY..ABBREV. = "lipid_category",
    target_LIPID.SUBCLASS..ABBREV. = "lipid_subclass",
    target_mz = "database_mz",
    mzmed = "feature_mz",
    target_rt_adjusted = "database_rt_adjusted",
    target_rt_sd = "database_rt_uncorrected",
    rtmed = "feature_rt",
    score_rt = "rt_error_sec",
    isopeak_sim = "isotope_similarity",
    isopeak_count = "isotope_peak_count"
  )
  for (old in names(rename_map)) {
    if (old %in% colnames(out)) {
      colnames(out)[colnames(out) == old] <- rename_map[[old]]
    }
  }

  keep_cols <- c(
    "polarity", "annotation_status", "interpretation",
    "candidate_generated", "feature_scope", "included_in_curated_metrics",
    "feature_id", "lipid_name", "lipid_category", "lipid_subclass",
    "adducts", "feature_mz", "database_mz", "ppm_error", "feature_rt",
    "database_rt_adjusted", "database_rt_uncorrected", "rt_error_sec",
    "isotope_similarity", "isotope_peak_count", "adduct_ratio",
    "adduct_count", "nb_annotations"
  )
  keep_cols <- keep_cols[keep_cols %in% colnames(out)]
  out <- out[, keep_cols, drop = FALSE]

  numeric_cols <- intersect(
    c("feature_mz", "database_mz", "ppm_error", "feature_rt",
      "database_rt_adjusted", "database_rt_uncorrected", "rt_error_sec",
      "isotope_similarity", "adduct_ratio"),
    colnames(out)
  )
  for (col in numeric_cols) out[[col]] <- round(as.numeric(out[[col]]), 5)

  rownames(out) <- NULL
  out
}

evaluate_publication_phase <- function(candidate_annotations, truth, phase) {
  cand <- pair_metadata(candidate_annotations)
  truth_meta <- pair_metadata(truth)

  cand_keys <- unique(cand$pair_key)
  truth_keys <- unique(truth_meta$pair_key)
  truth_features <- unique(truth_meta$feature_id)

  curated_candidates <- cand[cand$feature_id %in% truth_features, , drop = FALSE]
  curated_candidate_keys <- unique(curated_candidates$pair_key)
  uncurated_candidates <- cand[!cand$feature_id %in% truth_features, , drop = FALSE]

  curated_matches <- intersect(cand_keys, truth_keys)
  alternative_curated <- setdiff(curated_candidate_keys, truth_keys)
  missed_curated <- setdiff(truth_keys, cand_keys)

  curated_precision <- if (length(curated_candidate_keys) > 0) {
    length(curated_matches) / length(curated_candidate_keys)
  } else {
    NA_real_
  }
  curated_recall <- if (length(truth_keys) > 0) {
    length(curated_matches) / length(truth_keys)
  } else {
    NA_real_
  }
  curated_f1 <- if (is.na(curated_precision) || is.na(curated_recall) ||
                    (curated_precision + curated_recall) == 0) {
    NA_real_
  } else {
    2 * curated_precision * curated_recall /
      (curated_precision + curated_recall)
  }
  all_feature_curated_match_rate <- if (length(cand_keys) > 0) {
    length(curated_matches) / length(cand_keys)
  } else {
    NA_real_
  }

  by_feature <- table(cand$feature_id)
  curated_by_feature <- table(curated_candidates$feature_id)

  data.frame(
    phase = phase,
    candidate_annotations_all_features = length(cand_keys),
    candidate_features_all = length(unique(cand$feature_id)),
    candidate_lipids_all_features = length(unique(cand$lipid_name)),
    curated_truth_annotations = length(truth_keys),
    curated_truth_features = length(truth_features),
    curated_matches = length(curated_matches),
    alternative_candidates_on_curated_features = length(alternative_curated),
    candidates_on_uncurated_features = length(unique(uncurated_candidates$pair_key)),
    uncurated_features_with_candidates = length(unique(uncurated_candidates$feature_id)),
    missed_curated_annotations = length(missed_curated),
    curated_feature_precision = curated_precision,
    curated_feature_recall = curated_recall,
    curated_feature_f1 = curated_f1,
    all_feature_curated_match_rate = all_feature_curated_match_rate,
    ambiguous_candidate_features_all = sum(by_feature > 1),
    ambiguous_candidate_features_curated = sum(curated_by_feature > 1),
    mean_candidates_per_feature_all = if (length(by_feature)) {
      mean(as.numeric(by_feature))
    } else {
      0
    },
    max_candidates_per_feature_all = if (length(by_feature)) {
      max(as.numeric(by_feature))
    } else {
      0
    },
    stringsAsFactors = FALSE
  )
}

safe_sheet_name <- function(phase) {
  lookup <- c(
    mz_only_all_ranks = "mz_only",
    mz_rt_uncorrected_all_ranks = "mz_rt_uncorrected",
    mz_rt_corrected_all_ranks = "mz_rt_corrected",
    rank1_mz_only = "rank1_mz_only",
    rank1_mz_rt_corrected = "rank1_mz_rt",
    rank1_isotope_filter = "isotope_filter",
    rank1_isotope_adduct_scored = "adduct_scored",
    rank1_isotope_with_adduct_support = "adduct_required",
    automatic_ambiguity_resolution = "ambiguity_auto",
    manual_curation = "manual_truth",
    qc_rsd_filtered = "qc_rsd_filtered"
  )
  out <- unname(lookup[phase])
  if (is.na(out)) out <- gsub("[^A-Za-z0-9_]", "_", phase)
  substr(out, 1, 31)
}

phase_rt_fit_method <- function(phase) {
  ifelse(
    phase %in% c("mz_rt_corrected_all_ranks",
                 "rank1_mz_rt_corrected",
                 "rank1_isotope_filter",
                 "rank1_isotope_adduct_scored",
                 "rank1_isotope_with_adduct_support",
                 "automatic_ambiguity_resolution",
                 "manual_curation",
                 "qc_rsd_filtered"),
    "scam_monotone_pspline",
    "not_used"
  )
}

rt_fit_r_squared <- function(fit) {
  if (inherits(fit, "scam")) {
    summary(fit)$r.sq
  } else {
    summary(fit)$r.squared
  }
}

rt_fit_diagnostics <- function(rt_fit, lipid_database, polarity) {
  detected <- !is.na(rt_fit$ref_rt) & !is.na(rt_fit$exp_rt)
  pred <- rep(NA_real_, length(rt_fit$ref_rt))
  if (any(detected)) {
    pred[detected] <- as.numeric(predict(
      rt_fit$fit,
      newdata = data.frame(ref_rt = rt_fit$ref_rt[detected])
    ))
  }
  residual <- rt_fit$exp_rt - pred
  abs_residual <- abs(residual[detected])
  rng <- attr(rt_fit$fit, "ref_rt_range")
  extrapolated <- lipid_database$rt_sd < rng[1] | lipid_database$rt_sd > rng[2]

  data.frame(
    polarity = if (polarity == "pos") "positive" else "negative",
    rt_fit_method = if (inherits(rt_fit$fit, "scam")) {
      "scam_monotone_pspline"
    } else {
      "poly_lm"
    },
    reference_lipids_total = length(rt_fit$ref_rt),
    reference_lipids_detected = sum(detected),
    reference_lipids_missing = sum(!detected),
    r_squared = rt_fit_r_squared(rt_fit$fit),
    mean_abs_residual_sec = if (length(abs_residual)) {
      mean(abs_residual)
    } else {
      NA_real_
    },
    median_abs_residual_sec = if (length(abs_residual)) {
      median(abs_residual)
    } else {
      NA_real_
    },
    max_abs_residual_sec = if (length(abs_residual)) {
      max(abs_residual)
    } else {
      NA_real_
    },
    calibration_min_rt_sec = rng[1],
    calibration_max_rt_sec = rng[2],
    database_entries = nrow(lipid_database),
    database_entries_extrapolated = sum(extrapolated, na.rm = TRUE),
    database_entries_extrapolated_pct =
      sum(extrapolated, na.rm = TRUE) / nrow(lipid_database),
    stringsAsFactors = FALSE
  )
}

rt_reference_details <- function(rt_fit, intern_standard, polarity) {
  label_col <- first_present(
    intern_standard,
    c("short_name", "lipid.name", "lipid_name", "Lipid", "Name")
  )
  labels <- if (!is.na(label_col)) {
    as.character(intern_standard[[label_col]])
  } else {
    paste0("reference_", seq_along(rt_fit$ref_rt))
  }
  pred <- rep(NA_real_, length(rt_fit$ref_rt))
  detected <- !is.na(rt_fit$ref_rt) & !is.na(rt_fit$exp_rt)
  if (any(detected)) {
    pred[detected] <- as.numeric(predict(
      rt_fit$fit,
      newdata = data.frame(ref_rt = rt_fit$ref_rt[detected])
    ))
  }

  data.frame(
    polarity = if (polarity == "pos") "positive" else "negative",
    reference_lipid = labels,
    database_rt_sec = rt_fit$ref_rt,
    observed_rt_sec = rt_fit$exp_rt,
    fitted_rt_sec = pred,
    residual_sec = rt_fit$exp_rt - pred,
    detected = detected,
    stringsAsFactors = FALSE
  )
}

evaluate_polarity <- function(polarity) {
  mode_label <- if (polarity == "pos") "positive" else "negative"
  polarity_dir <- file.path(study_dir, mode_label)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(polarity_dir)

  message("\nEvaluating ", mode_label, " mode for study_id = ", study_id)

  folders <- setup_folders(polarity = polarity)
  is_file <- if (polarity == "pos") {
    "pos_lipid_reference_set.xlsx"
  } else {
    "neg_lipid_reference_set.xlsx"
  }
  db_sheet <- if (polarity == "pos") 4L else 5L
  rt_col <- if (polarity == "pos") {
    "RT ESI(+) (min)"
  } else {
    paste0("RT ESI(", intToUtf8(0x2013), ") (min)")
  }

  mse_save_path <- file.path(
    folders$objects,
    paste0(study_id, "_preprocessed_mse_", polarity)
  )
  res_save_path <- file.path(
    folders$objects,
    paste0(study_id, "_preprocessed_res_", polarity)
  )

  mse <- readMsObject(XcmsExperiment(), AlabasterParam(path = mse_save_path))
  res <- readObject(res_save_path)

  intern_standard <- prepare_reference_lipids(
    is_file,
    rt_window_left = 30,
    rt_window_right = ifelse(polarity == "pos", 15, 20)
  )
  eic_is <- extract_is_eics(mse, intern_standard)

  group <- eic_is$sample_type
  group[group == "QC"] <- NA
  group <- as.factor(group)
  param_group <- PeakDensityParam(
    sampleGroups = group,
    minFraction = 0.5,
    binSize = 0.01,
    ppm = PPM,
    bw = 2
  )

  rt_fit <- fit_rt_correction(
    eic_is,
    intern_standard,
    param_group,
    output_dir = folders$ref_lipid
  )

  lipid_database <- prepare_lipid_database(
    db_path = "../LipidDatabase_R.xlsx",
    sheet = db_sheet,
    polarity = polarity,
    rt_col = rt_col,
    rt_fit = rt_fit
  )

  match_result <- match_features_to_database(
    res,
    lipid_database,
    ppm = MATCH_PPM,
    rt_tol = MATCH_RT_TOL
  )

  mtched_data <- match_result$mtched_data
  query <- match_result$query

  candidate_annotations_by_phase <- list(
    mz_only_all_ranks = match_lipid_candidates(
      res, lipid_database,
      ppm = MATCH_PPM,
      use_rt = FALSE,
      phase_name = "mz_only_all_ranks"
    ),
    mz_rt_uncorrected_all_ranks = match_lipid_candidates(
      res, lipid_database,
      ppm = MATCH_PPM,
      rt_tol = MATCH_RT_TOL,
      rt_col = "rt_sd",
      use_rt = TRUE,
      phase_name = "mz_rt_uncorrected_all_ranks"
    ),
    mz_rt_corrected_all_ranks = match_lipid_candidates(
      res, lipid_database,
      ppm = MATCH_PPM,
      rt_tol = MATCH_RT_TOL,
      rt_col = "rt_adjusted",
      use_rt = TRUE,
      phase_name = "mz_rt_corrected_all_ranks"
    ),
    rank1_mz_only = match_lipid_candidates(
      res, lipid_database,
      ppm = MATCH_PPM,
      rank_filter = 1,
      use_rt = FALSE,
      phase_name = "rank1_mz_only"
    ),
    rank1_mz_rt_corrected = mtched_data
  )

  iso_result <- calculate_isotope_similarity(
    mse,
    mtched_data,
    polarity = polarity,
    isopeak_threshold = 2,
    similarity_threshold = ISOPEAK_SIM_THRESHOLD
  )
  mtched_data <- iso_result$mtched_data
  candidate_annotations_by_phase$rank1_isotope_filter <- mtched_data

  mtched_data <- match_adducts(
    mtched_data,
    lipid_database,
    query,
    ppm = MATCH_PPM,
    rt_tol = 5
  )
  candidate_annotations_by_phase$rank1_isotope_adduct_scored <- mtched_data
  candidate_annotations_by_phase$rank1_isotope_with_adduct_support <-
    mtched_data[mtched_data$adduct_ratio > 0, ]

  key <- paste(mtched_data$target_lipid_name_unique,
               mtched_data$target_adduct, sep = "_")
  rm_vector <- c()
  for (i in key[duplicated(key)]) {
    idxs <- which(key == i)
    count <- mtched_data$adduct_ratio[idxs]
    if (which.max(count) == which.min(count)) next
    best_idx <- idxs[which.max(count)]
    rm_vector <- c(rm_vector, setdiff(idxs, best_idx))
  }
  if (length(rm_vector) > 0) {
    mtched_data <- mtched_data[-rm_vector, ]
  }

  key <- paste(mtched_data$target_lipid_name_unique,
               mtched_data$target_adduct, sep = "_")
  dupe_keys <- unique(key[duplicated(key)])
  rows_to_keep <- rep(TRUE, nrow(mtched_data))
  for (k in dupe_keys) {
    idxs <- which(key == k)
    scores <- abs(mtched_data$score_rt[idxs])
    best_score <- min(scores, na.rm = TRUE)
    to_remove <- idxs[(scores - best_score) > 4]
    if (length(to_remove) > 0) rows_to_keep[to_remove] <- FALSE
  }
  mtched_data <- mtched_data[rows_to_keep, ]
  mtched_data <- resolve_sm_isomers(mtched_data)
  candidate_annotations_by_phase$automatic_ambiguity_resolution <- mtched_data

  lipid_amb_res <- readxl::read_excel("lipid_ambiguity_resolution.xlsx")
  lipid_amb_res$keep_row <- as.logical(lipid_amb_res$keep_row)
  mtch_to_be_rm <- lipid_amb_res$ntch_idx[!lipid_amb_res$keep_row]
  mtched_data <- mtched_data[!mtched_data$ntch_idx %in% mtch_to_be_rm, ]

  feature_amb_res <- readxl::read_excel("feature_ambiguity_resolution.xlsx")
  feature_amb_res$keep_row <- as.logical(feature_amb_res$keep_row)
  mtch_to_be_rm_ft <- feature_amb_res$ntch_idx[!feature_amb_res$keep_row]
  mtched_data <- mtched_data[!mtched_data$ntch_idx %in% mtch_to_be_rm_ft, ]

  mtched_data$nb_annotations <- 1
  fids_ambiguous <- mtched_data$feature_id[duplicated(mtched_data$feature_id)]
  for (i in unique(fids_ambiguous)) {
    idxs <- which(mtched_data$feature_id == i)
    if (length(idxs) > 1) {
      combined_names <- paste(unique(mtched_data$target_lipid_name_unique[idxs]),
                              collapse = "; ")
      combined_adducts <- paste(unique(mtched_data$target_adduct[idxs]),
                                collapse = "; ")
      mtched_data$target_lipid_name_unique[idxs[1]] <- combined_names
      mtched_data$target_adduct[idxs[1]] <- combined_adducts
      mtched_data$nb_annotations[idxs[1]] <- length(idxs)
      mtched_data <- mtched_data[-idxs[-1], ]
    }
  }

  truth_annotations <- mtched_data
  candidate_annotations_by_phase$manual_curation <- truth_annotations
  istd_filter <- filter_internal_standards_for_metrics(
    candidate_annotations_by_phase,
    truth_annotations
  )
  candidate_annotations_by_phase <- istd_filter$candidates
  truth_annotations <- istd_filter$truth

  metrics_path <- file.path(
    folders$objects,
    paste0(study_id, "_annotation_metrics_", polarity, ".csv")
  )
  truth_path <- file.path(
    folders$objects,
    paste0(study_id, "_truth_annotations_manual_", polarity, ".csv")
  )

  metrics <- do.call(
    rbind,
    Map(
      f = function(candidate_annotations, phase) {
        evaluate_publication_phase(candidate_annotations, truth_annotations, phase)
      },
      candidate_annotations_by_phase,
      names(candidate_annotations_by_phase)
    )
  )
  rownames(metrics) <- NULL
  metrics$rt_fit_method <- phase_rt_fit_method(metrics$phase)
  write.csv(metrics, metrics_path, row.names = FALSE, na = "")
  export_annotation_truth(truth_annotations, truth_path, polarity = mode_label)

  phase_details <- Map(
    f = function(candidate_annotations, phase) {
      make_phase_detail(
        candidate_annotations,
        truth_annotations,
        phase,
        mode_label
      )
    },
    candidate_annotations_by_phase,
    names(candidate_annotations_by_phase)
  )

  qc_metrics <- NULL
  if (nrow(truth_annotations) > 0) {
    volume_factors <- list("W" = 2.5, "P" = 1.08)
    res_qc <- apply_volume_correction(
      res,
      volume_factors,
      assay_name = "raw",
      new_assay_name = "raw_corr"
    )
    res_qc <- apply_volume_correction(
      res_qc,
      volume_factors,
      assay_name = "raw_filled",
      new_assay_name = "corr_filled"
    )

    annotated_res <- res_qc[rownames(res_qc) %in% truth_annotations$feature_id, ]
    annotation_cols <- intersect(
      c("target_lipid.name", "target_lipid_name_unique",
        "target_LIPID.CATEGORY..ABBREV.",
        "target_LIPID.SUBCLASS..ABBREV.", "target_Adduct",
        "target_IS_norm"),
      colnames(truth_annotations)
    )
    rowData(annotated_res) <- cbind(
      rowData(annotated_res),
      truth_annotations[
        match(rownames(annotated_res), truth_annotations$feature_id),
        annotation_cols,
        drop = FALSE
      ]
    )
    annotated_res <- normalize_by_is(
      annotated_res,
      input_assay = "corr_filled",
      output_assay = "ISnorm_filled",
      is_col = "target_IS_norm",
      lipid_col = "target_lipid.name"
    )
    annotated_res <- filter_by_qc_rsd(annotated_res, threshold = 0.3)

    qc_filtered_annotations <- data.frame(
      feature_id = rownames(annotated_res),
      as.data.frame(rowData(annotated_res)),
      stringsAsFactors = FALSE
    )
    qc_metrics <- evaluate_publication_phase(
      qc_filtered_annotations,
      truth_annotations,
      phase = "qc_rsd_filtered"
    )
    qc_metrics$rt_fit_method <- phase_rt_fit_method(qc_metrics$phase)
    phase_details$qc_rsd_filtered <- make_phase_detail(
      qc_filtered_annotations,
      truth_annotations,
      "qc_rsd_filtered",
      mode_label
    )
    write.csv(
      qc_metrics,
      file.path(folders$objects,
                paste0(study_id, "_qc_retention_metrics_", polarity, ".csv")),
      row.names = FALSE,
      na = ""
    )
  }

  message("Saved metrics: ", metrics_path)
  message("Saved manual truth: ", truth_path)
  invisible(list(
    metrics = metrics,
    qc_metrics = qc_metrics,
    details = phase_details,
    rt_diagnostics = rt_fit_diagnostics(rt_fit, lipid_database, polarity),
    rt_references = rt_reference_details(rt_fit, intern_standard, polarity),
    istd_filter_summary = data.frame(
      polarity = mode_label,
      istd_filter$summary,
      stringsAsFactors = FALSE
    )
  ))
}

result_pos <- evaluate_polarity("pos")
result_neg <- evaluate_polarity("neg")

metrics_pos <- result_pos$metrics
metrics_neg <- result_neg$metrics
metrics_pos$polarity <- "positive"
metrics_neg$polarity <- "negative"
combined_metrics <- rbind(metrics_pos, metrics_neg)
combined_metrics <- combined_metrics[
  c("polarity", setdiff(colnames(combined_metrics), "polarity"))
]

rt_diagnostics <- bind_rows_aligned(list(
  result_pos$rt_diagnostics,
  result_neg$rt_diagnostics
))
rt_references <- bind_rows_aligned(list(
  result_pos$rt_references,
  result_neg$rt_references
))
istd_filter_summary <- bind_rows_aligned(list(
  result_pos$istd_filter_summary,
  result_neg$istd_filter_summary
))

qc_summary <- rbind(
  cbind(polarity = "positive", result_pos$qc_metrics),
  cbind(polarity = "negative", result_neg$qc_metrics)
)
summary_all <- bind_rows_aligned(list(combined_metrics, qc_summary))
summary_all$phase <- vapply(
  summary_all$phase,
  safe_sheet_name,
  character(1)
)

summary_cols <- c(
  "polarity", "phase", "rt_fit_method",
  "candidate_annotations_all_features", "candidate_features_all",
  "candidate_lipids_all_features", "curated_truth_annotations",
  "curated_truth_features", "curated_matches",
  "alternative_candidates_on_curated_features",
  "candidates_on_uncurated_features", "uncurated_features_with_candidates",
  "missed_curated_annotations", "all_feature_curated_match_rate",
  "ambiguous_candidate_features_all", "mean_candidates_per_feature_all",
  "max_candidates_per_feature_all"
)
summary_cols <- summary_cols[summary_cols %in% colnames(summary_all)]
summary_all <- summary_all[, summary_cols, drop = FALSE]

db_only_cols <- c(
  "polarity", "phase", "rt_fit_method",
  "curated_truth_annotations", "curated_truth_features",
  "curated_matches", "alternative_candidates_on_curated_features",
  "missed_curated_annotations", "curated_feature_precision",
  "curated_feature_recall", "curated_feature_f1",
  "ambiguous_candidate_features_curated"
)
db_only_cols <- db_only_cols[db_only_cols %in% colnames(summary_all) |
                               db_only_cols %in% colnames(bind_rows_aligned(list(combined_metrics, qc_summary)))]
db_only_metrics <- bind_rows_aligned(list(combined_metrics, qc_summary))
db_only_metrics$phase <- vapply(
  db_only_metrics$phase,
  safe_sheet_name,
  character(1)
)
db_only_metrics <- db_only_metrics[, db_only_cols, drop = FALSE]

combined_metrics_path <- file.path(
  study_dir,
  "objects",
  paste0(study_id, "_annotation_metrics_pos_neg.csv")
)
dir.create(dirname(combined_metrics_path), recursive = TRUE, showWarnings = FALSE)
write.csv(combined_metrics, combined_metrics_path, row.names = FALSE, na = "")
message("Saved combined POS/NEG metrics: ", combined_metrics_path)

merge_metrics_path <- file.path(
  study_dir,
  "objects",
  paste0(study_id, "_annotation_merge_metrics.csv")
)
if (file.exists(merge_metrics_path)) {
  merge_metrics <- read.csv(
    merge_metrics_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  warning(
    "Merge metrics were not found at ", merge_metrics_path,
    ". Render POS_NEG_merge.qmd before generating the metrics workbook."
  )
  merge_metrics <- data.frame(
    scope = "merged",
    phase = "not_available",
    status = "Render POS_NEG_merge.qmd before generating metrics.",
    stringsAsFactors = FALSE
  )
}

merge_details_path <- file.path(
  study_dir,
  "objects",
  paste0(study_id, "_annotation_merge_details.csv")
)
if (file.exists(merge_details_path)) {
  merge_details <- read.csv(
    merge_details_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  warning(
    "Merged compound details were not found at ", merge_details_path,
    ". Render POS_NEG_merge.qmd before generating the metrics workbook."
  )
  merge_details <- data.frame(
    scope = "merged",
    phase = "not_available",
    status = "Render POS_NEG_merge.qmd before generating metrics.",
    stringsAsFactors = FALSE
  )
}

rt_diagnostics_path <- file.path(
  study_dir,
  "objects",
  paste0(study_id, "_rt_fit_diagnostics_pos_neg.csv")
)
write.csv(rt_diagnostics, rt_diagnostics_path, row.names = FALSE, na = "")
message("Saved RT fit diagnostics: ", rt_diagnostics_path)

detail_names <- unique(c(names(result_pos$details), names(result_neg$details)))
detail_sheets <- lapply(detail_names, function(phase) {
  rows <- list(result_pos$details[[phase]], result_neg$details[[phase]])
  bind_rows_aligned(rows)
})
names(detail_sheets) <- vapply(detail_names, safe_sheet_name, character(1))

readme <- data.frame(
  field = c("truth_definition", "annotation_status_true",
            "annotation_status_false",
            "annotation_status_unknown",
            "annotation_status_missed",
            "Summary",
            "DB_only_metrics",
            "Merge_metrics",
            "merged",
            "RT_fit_diagnostics",
            "ISTD_filter",
            "all_feature_curated_match_rate",
            "why_unknown_is_not_wrong", "adduct_required_note",
            "polarity_scope"),
  description = c(
    "Truth is the final manually curated feature-lipid annotation table after ambiguity curation.",
    "The candidate annotation matches the manually curated feature-lipid pair.",
    "The candidate is on a curated feature, but this lipid assignment was not selected in the manual truth set.",
    "The candidate is on a feature absent from the curated truth set; it is not counted as wrong.",
    "The manually curated annotation was not recovered in this phase.",
    "All-candidate overview. Positive and negative ionization modes are kept separate, and unknown candidates outside the curated set are retained.",
    "Precision, recall and F1 after excluding candidates on features outside the curated truth set.",
    "Operational merge counts, including positive-only, negative-only, overlap and union annotation totals for a POS/NEG Venn diagram. These are retention metrics, not annotation precision/recall metrics.",
    "One row per unique lipid annotation. Overlapping POS/NEG annotations are listed first with source details side by side, followed by positive-only and negative-only annotations.",
    "Diagnostics for the merged main RT correction method: monotone scam P-spline fitting, including reference-lipid residuals and database extrapolation counts.",
    "Counts of injected-standard annotations/features removed before computing annotation metrics.",
    "Curated matches divided by all candidate annotations, including uncurated features. This is a candidate-burden indicator, not precision.",
    "The absence of a curated annotation for a feature does not prove the candidate is wrong; it may simply not have been manually resolved.",
    "The adduct_required sheet is a strict comparison phase requiring adduct_ratio > 0; default workflow uses adduct_ratio as support/ambiguity evidence, not as a hard filter.",
    "Summary and DB_only_metrics keep positive and negative ionization modes as separate rows; Merge_metrics reports the later combined scope."
  ),
  stringsAsFactors = FALSE
)

workbook_path <- file.path(
  study_dir,
  "objects",
  paste0(study_id, "_annotation_metrics_publication.xlsx")
)
write_xlsx(
  c(list(README = readme, Summary = summary_all,
         DB_only_metrics = db_only_metrics,
         Merge_metrics = merge_metrics,
         merged = merge_details,
         RT_fit_diagnostics = rt_diagnostics,
         ISTD_filter = istd_filter_summary,
         RT_reference_residuals = rt_references), detail_sheets),
  path = workbook_path
)
message("Saved publication workbook: ", workbook_path)
