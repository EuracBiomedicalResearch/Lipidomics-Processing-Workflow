# Run from the repository root: Rscript tests/test_annotation_reporting.R
source("R/annotation_reporting.R")
make_input <- function(ids, labels, standards = rep(FALSE, length(ids))) {
  data.frame(feature_id = ids, target_lipid_name_unique = labels,
    target_lipid.name = ifelse(standards, "IS(d7)", labels),
    target_IS_norm = ifelse(standards, "IS(d7)", "another IS"),
    mzmed = seq_along(ids), rtmed = seq_along(ids) * 10)
}
make_report <- function(mode) {
  report <- new_annotation_report("TEST", mode)
  candidates <- make_input(c("FT1", "FT1", "FT2", "FT3", "IS1", "IS1"),
    c("A_1", "B_1", "C_1", "D_1", "IS_1", "alternative_1"),
    c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  truth <- make_input(c("FT1", "FT2", "IS1"), c("A_1", "C_1; E_1", "IS_1"),
    c(FALSE, FALSE, TRUE))
  report <- capture_annotation_phase(report, "preprocessed_features",
    make_input(c("FT1", "FT2", "FT3", "IS1", "FT4"), rep(NA_character_, 5)), annotations = FALSE)
  report <- capture_annotation_phase(report, "rank1_mz",
    rbind(candidates, make_input("FT4", "F_1")))
  report <- insert_rank1_unadjusted_phase(report,
    make_input(c("FT1", "FT3", "FT4"), c("B_1", "D_1", "F_1")))
  for (phase in c("rank1_mz_rt", "isotope_filter", "adduct_scored", "ambiguity_auto"))
    report <- capture_annotation_phase(report, phase, candidates)
  report <- capture_annotation_phase(report, "manual_curation", truth)
  report <- capture_annotation_phase(report, "normalization", truth)
  report <- capture_annotation_phase(report, "qc_rsd_filtered", truth[c(1, 3), ])
  report <- capture_annotation_phase(report, "within_mode_adduct_resolution", truth[c(1, 3), ])
  capture_annotation_phase(report, "qc_samples_removed", truth[c(1, 3), ])
}
positive <- make_report("positive")
negative <- make_report("negative")
merged <- new_annotation_report("TEST", "merged")
merged$phases$merge_input <- list(data = rbind(positive$phases$qc_samples_removed$data,
  negative$phases$qc_samples_removed$data), annotations = TRUE)
merged$phases$resolved_annotation_preference <- merged$phases$merge_input
merged$phases$duplicate_resolution <- list(data = merged$phases$merge_input$data[c(1, 2, 3), ],
  annotations = TRUE)
merged$phases$internal_standard_removal <- list(
  data = merged$phases$duplicate_resolution$data[c(1, 3), ], annotations = TRUE)
original <- list(positive, negative, merged)
tables <- annotation_metrics_tables(positive, negative, merged)
s <- tables$Summary
stopifnot(identical(original, list(positive, negative, merged)),
  identical(names(s), c("polarity", "phase", "features", "removed_features")),
  s$features[s$phase == "ambiguity_resolution"][1] == 2L,
  !any(s$phase %in% c("manual_curation", "rank1_mz_rt_unadjusted", "qc_samples_removed", "adduct_scored")),
  tail(s$features, 1) == 2L,
  tail(s$removed_features, 1) == 0L,
  !any(c("README", "Merge_metrics", "ISTD_filter", "manual_curation",
         "Curated_reference_comparison") %in% names(tables)),
  nrow(tables$final_annotations) == 2L,
  all(vapply(tables[-c(1, 2)], function(d) !any(d$feature_id == "IS1"), logical(1))))
summary <- annotation_report_summary(positive)
stopifnot(identical(summary$features[1:3], c(4L, 4L, 3L)),
  summary$annotation_pairs[summary$phase == "manual_curation"] == 3L,
  summary$ambiguous_features[summary$phase == "manual_curation"] == 1L,
  "adduct_scored" %in% names(tables), "rank1_mz_rt_unadjusted" %in% names(tables))
# The parallel comparison must not change sequential removed-feature counts.
without_comparison <- positive
without_comparison$phases$rank1_mz_rt_unadjusted <- NULL
stopifnot(identical(annotation_report_summary(without_comparison), summary))
comparison <- tables$Final_reference_comparison
stopifnot(identical(names(comparison), c("polarity", "phase", "matching_pairs",
  "alternative_pairs", "outside_reference_pairs", "reference_retention",
  "nonreference_reduction", "within_reference_discrimination")),
  identical(comparison$phase[1:7], c("rank1_mz", "rank1_mz_rt_unadjusted", "rank1_mz_rt",
    "isotope_filter", "ambiguity_auto", "ambiguity_resolution", "qc_rsd_filtered")))
rt <- comparison[comparison$phase == "rank1_mz_rt", ][1, ]
unadjusted <- comparison[comparison$phase == "rank1_mz_rt_unadjusted", ][1, ]
stopifnot(rt$matching_pairs == 1L, rt$alternative_pairs == 1L,
  rt$outside_reference_pairs == 2L, rt$reference_retention == 1,
  abs(rt$nonreference_reduction - 1/3) < 1e-12,
  rt$within_reference_discrimination == 0.5,
  unadjusted$reference_retention == 0, unadjusted$within_reference_discrimination == 0,
  all(comparison$nonreference_reduction[comparison$phase == "rank1_mz"] == 0),
  all(comparison$reference_retention[comparison$phase == "qc_rsd_filtered"] == 1),
  all(comparison$within_reference_discrimination[comparison$phase == "qc_rsd_filtered"] == 1))
qc_without_standards <- positive
qc_without_standards$phases$qc_rsd_filtered$data <-
  qc_without_standards$phases$qc_rsd_filtered$data[!positive$phases$qc_rsd_filtered$data$.is_standard, ]
stopifnot(identical(annotation_reference_comparison(qc_without_standards),
                    annotation_reference_comparison(positive)))
# Distinct pairs, not features or duplicate rows, define T_s, C_s and U_s.
ambiguous <- new_annotation_report("TEST", "positive")
candidates <- make_input(c("FT1", "FT1", "FT2", "FT3"), c("A_1; B_1; Z_1", "A_1", "D_1", "E_1"))
for (phase in c("rank1_mz", "rank1_mz_rt"))
  ambiguous <- capture_annotation_phase(ambiguous, phase, candidates)
for (phase in c("manual_curation", "qc_rsd_filtered"))
  ambiguous <- capture_annotation_phase(ambiguous, phase,
    make_input(c("FT1", "FT2"), c("A_1; B_1", "C_1")))
ambiguous <- capture_annotation_phase(ambiguous, "within_mode_adduct_resolution", make_input("FT1", "A_1; B_1"))
a <- annotation_reference_comparison(ambiguous)
stopifnot(a$matching_pairs[1] == 2L, a$alternative_pairs[1] == 2L,
  a$outside_reference_pairs[1] == 1L,
  abs(a$reference_retention[1] - 2/3) < 1e-12,
  a$within_reference_discrimination[1] == 0.5,
  attr(a, "reference_pairs") == 3L, attr(a, "reference_features") == 2L,
  attr(a, "baseline_outside_pairs") == 1L)
empty <- new_annotation_report("TEST", "positive")
for (phase in c("rank1_mz", "rank1_mz_rt", "manual_curation", "qc_rsd_filtered"))
  empty <- capture_annotation_phase(empty, phase, make_input(character(), character()))
rates <- c("reference_retention", "nonreference_reduction", "within_reference_discrimination")
stopifnot(all(annotation_report_summary(empty)$features == 0L),
  all(is.na(as.matrix(annotation_reference_comparison(empty)[rates]))))
no_outside <- positive
no_outside$phases$rank1_mz$data <- no_outside$phases$qc_rsd_filtered$data
stopifnot(all(is.na(annotation_reference_comparison(no_outside)$nonreference_reduction)))
missing <- positive
missing$phases$qc_rsd_filtered <- NULL
stopifnot(inherits(try(annotation_reference_comparison(missing), silent = TRUE), "try-error"))
missing <- positive
missing$phases$rank1_mz <- NULL
stopifnot(inherits(try(annotation_reference_comparison(missing), silent = TRUE), "try-error"))
payload <- annotation_metrics_payload(tables)
stopifnot(length(payload$sheets) == 14L,
  grepl("qc_rsd_filtered", payload$sheets[[2]]$note, fixed = TRUE),
  grepl("|G| = 1", payload$sheets[[2]]$note, fixed = TRUE),
  grepl("U0 = 3", payload$sheets[[2]]$note, fixed = TRUE),
  all(vapply(payload$sheets, function(x) length(x$columns) == length(x$descriptions) &&
    all(nzchar(x$descriptions)), logical(1))))
if (requireNamespace("jsonlite", quietly = TRUE) && length(commandArgs(TRUE)))
  jsonlite::write_json(payload, commandArgs(TRUE)[1], auto_unbox = TRUE, na = "null", null = "null")
cat("PASS: ISTD-free counts, pair partition and formulas, phase labels, parallel comparisons, missing and empty inputs.\n")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  path <- tempfile(fileext = ".xlsx")
  write_annotation_metrics_workbook(positive, negative, merged, path)
  stopifnot(identical(openxlsx::getSheetNames(path), names(tables)))
  summary <- openxlsx::read.xlsx(path, sheet = "Summary", startRow = 3, colNames = FALSE)
  headers <- openxlsx::read.xlsx(path, sheet = "Final_reference_comparison", rows = 3, colNames = FALSE)
  descriptions <- openxlsx::read.xlsx(path, sheet = "Final_reference_comparison", rows = 4, colNames = FALSE)
  exported <- openxlsx::read.xlsx(path, sheet = "Final_reference_comparison", startRow = 5, colNames = FALSE)
  stopifnot(identical(unname(unlist(headers)), names(comparison)),
    length(unlist(descriptions)) == 8L, all(nzchar(unlist(descriptions))),
    ncol(summary) == 4L, ncol(exported) == 8L, is.numeric(summary[[3]]),
    identical(as.integer(summary[[3]]), tables$Summary$features),
    identical(as.integer(summary[[4]]), tables$Summary$removed_features),
    isTRUE(all.equal(exported[[6]], comparison$reference_retention)),
    isTRUE(all.equal(exported[[7]], comparison$nonreference_reduction)),
    isTRUE(all.equal(exported[[8]], comparison$within_reference_discrimination)))
  # Undefined percentages are explicit N/A text, not errors or empty cells.
  empty$polarity <- "negative"
  write_annotation_metrics_workbook(positive, empty, merged, path)
  exported <- openxlsx::read.xlsx(path, sheet = "Final_reference_comparison", startRow = 5, colNames = FALSE)
  stopifnot(all(as.matrix(exported[exported[[1]] == "negative", 6:8]) == "N/A"))
  unlink(path)
  cat("PASS: XLSX numeric values, header descriptions, labels and zero-denominator N/A.\n")
  # Legacy keys and the absence of a retired within-mode stage remain valid.
  temporary_root <- tempfile()
  folder <- file.path(temporary_root, "applications", "TEST_study")
  for (mode in c("positive", "negative")) {
    dir.create(file.path(folder, mode, "objects"), recursive = TRUE)
    report <- if (mode == "positive") positive else negative
    report$phases$within_mode_adduct_resolution <- NULL
    suffix <- if (mode == "positive") "pos" else "neg"
    saveRDS(report, file.path(folder, mode, "objects", paste0("TEST_annotation_report_", suffix, ".rds")))
  }
  dir.create(file.path(folder, "objects"))
  saveRDS(merged, file.path(folder, "objects", "TEST_annotation_report_merged.rds"))
  regenerate_annotation_metrics(folder, "TEST")
  path <- file.path(folder, "objects", "TEST_annotation_metrics.xlsx")
  stopifnot(file.exists(path), !"within_mode_adduct_resolution" %in% openxlsx::getSheetNames(path))
  report$phases$qc_rsd_filtered <- NULL
  saveRDS(report, file.path(folder, "negative", "objects", "TEST_annotation_report_neg.rds"))
  stopifnot(inherits(try(regenerate_annotation_metrics(folder, "TEST"), silent = TRUE), "try-error"))
  unlink(temporary_root, recursive = TRUE)
  cat("PASS: regeneration accepts legacy keys and rejects missing QC snapshots.\n")
}
