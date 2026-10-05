# Run from the repository root: Rscript tests/test_annotation_reporting.R
source("R/annotation_reporting.R")
make_input <- function(ids, names, standards = rep(FALSE, length(ids))) {
  data.frame(feature_id = ids, target_lipid_name_unique = names,
             target_lipid.name = ifelse(standards, "IS(d7)", names),
             target_IS_norm = ifelse(standards, "IS(d7)", "another IS"),
             mzmed = seq_along(ids), rtmed = seq_along(ids) * 10)
}
make_report <- function(mode) {
  report <- new_annotation_report("TEST", mode)
  candidates <- make_input(c("FT1", "FT1", "FT2", "FT3", "IS1"),
                           c("A_1", "B_1", "C_1", "D_1", "IS_1"),
                           c(FALSE, FALSE, FALSE, FALSE, TRUE))
  truth <- make_input(c("FT1", "FT2", "IS1"), c("A_1", "C_1; E_1", "IS_1"),
                      c(FALSE, FALSE, TRUE))
  report <- capture_annotation_phase(report, "preprocessed_features",
    make_input(c("FT1", "FT2", "FT3", "IS1", "FT4"), rep(NA_character_, 5)), annotations = FALSE)
  mz_candidates <- rbind(candidates, make_input("FT4", "F_1"))
  report <- capture_annotation_phase(report, "rank1_mz", mz_candidates)
  for (phase in c("rank1_mz_rt", "isotope_filter", "adduct_scored", "ambiguity_auto"))
    report <- capture_annotation_phase(report, phase, candidates)
  report <- capture_annotation_phase(report, "manual_curation", truth)
  report <- capture_annotation_phase(report, "normalization", truth)
  report <- capture_annotation_phase(report, "qc_rsd_filtered", truth[c(1, 3), ])
  report <- capture_annotation_phase(report, "within_mode_adduct_resolution", truth[c(1, 3), ])
  report <- capture_annotation_phase(report, "qc_samples_removed", truth[c(1, 3), ])
  report
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
tables <- annotation_metrics_tables(positive, negative, merged)
s <- tables$Summary
stopifnot(identical(names(s), c("polarity", "phase", "features", "removed_features")),
          s$features[s$phase == "manual_curation"][1] == 3L,
          tail(s$features, 1) == 2L, # Identical IDs in different polarities are distinct.
          tail(s$removed_features, 1) == 1L,
          !any(c("README", "Merge_metrics", "ISTD_filter") %in% names(tables)),
          nrow(tables$final_annotations) == 2L)
# Internal counts remain available; only the workbook presentation is reduced.
internal_summary <- annotation_report_summary(positive)
stopifnot(internal_summary$annotation_pairs[internal_summary$phase == "manual_curation"] == 4L,
          internal_summary$ambiguous_features[internal_summary$phase == "manual_curation"] == 1L)
stopifnot(identical(s$phase[1:3], c("preprocessed_features", "rank1_mz", "rank1_mz_rt")),
          s$features[2] == 5L, s$features[3] == 4L, s$removed_features[3] == 1L,
          "rank1_mz" %in% names(tables))
comparison <- tables$Curated_reference_comparison
first <- comparison[comparison$phase == "rank1_mz_rt", ][1, ]
stopifnot(identical(names(comparison), c("polarity", "phase", "curated_matches",
          "uncurated_assignments", "missed_curated_pairs", "precision", "recall", "f1")),
          first$curated_matches == 2L,
          first$uncurated_assignments == 1L, first$missed_curated_pairs == 1L,
          abs(first$precision - 2/3) < 1e-12,
          abs(first$recall - 2/3) < 1e-12,
          all(comparison$precision[comparison$phase == "manual_curation"] == 1),
          all(comparison$recall[comparison$phase == "manual_curation"] == 1))
empty <- make_input(character(), character())
empty_report <- new_annotation_report("TEST", "positive")
empty_report <- capture_annotation_phase(empty_report, "rank1_mz_rt", empty)
empty_report <- capture_annotation_phase(empty_report, "manual_curation", empty)
stopifnot(all(annotation_report_summary(empty_report)$features == 0L),
          all(is.na(annotation_reference_comparison(empty_report)$precision)))
payload <- annotation_metrics_payload(tables)
stopifnot(length(payload$sheets) == 13L,
          all(vapply(payload$sheets, function(x) length(x$columns) == length(x$descriptions), logical(1))))
if (requireNamespace("jsonlite", quietly = TRUE) && length(commandArgs(TRUE)))
  jsonlite::write_json(payload, commandArgs(TRUE)[1], auto_unbox = TRUE, na = "null", null = "null")
cat("PASS: counts, ambiguities, cross-polarity IDs, standard exclusion, reference metrics and empty inputs.\n")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  path <- tempfile(fileext = ".xlsx")
  write_annotation_metrics_workbook(positive, negative, merged, path)
  stopifnot(identical(openxlsx::getSheetNames(path), names(tables)))
  summary <- openxlsx::read.xlsx(path, sheet = "Summary", startRow = 3, colNames = FALSE)
  summary_headers <- openxlsx::read.xlsx(path, sheet = "Summary", rows = 1, colNames = FALSE)
  comparison_headers <- openxlsx::read.xlsx(path, sheet = "Curated_reference_comparison",
                                            rows = 3, colNames = FALSE)
  exported_comparison <- openxlsx::read.xlsx(path, sheet = "Curated_reference_comparison",
                                             startRow = 5, colNames = FALSE)
  stopifnot(identical(unname(unlist(summary_headers)), names(tables$Summary)),
            identical(unname(unlist(comparison_headers)), names(comparison)),
            ncol(summary) == 4L, ncol(exported_comparison) == 8L,
            is.numeric(summary[[3]]),
            identical(as.integer(summary[[4]]), tables$Summary$removed_features),
            isTRUE(all.equal(exported_comparison[[6]], comparison$precision)),
            isTRUE(all.equal(exported_comparison[[7]], comparison$recall)),
            isTRUE(all.equal(exported_comparison[[8]], comparison$f1)),
            identical(as.integer(summary[[3]]), tables$Summary$features))
  final <- openxlsx::read.xlsx(path, sheet = "final_annotations", startRow = 3, colNames = FALSE)
  stopifnot(nrow(final) == 2L)
  unlink(path)
  cat("PASS: workflow XLSX export preserves numeric counts, sheet order and final rows.\n")
}
