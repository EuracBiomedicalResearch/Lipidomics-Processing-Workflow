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
  candidates <- rbind(candidates, make_input("IS1", "endogenous_candidate_1"))
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
          !any(s$phase %in% c("qc_samples_removed", "adduct_scored")),
          tail(s$features, 1) == 2L, # Identical IDs in different polarities are distinct.
          tail(s$removed_features, 1) == 1L,
          !any(c("README", "Merge_metrics", "ISTD_filter") %in% names(tables)),
          nrow(tables$final_annotations) == 2L)
# Internal counts remain available; only the workbook presentation is reduced.
internal_summary <- annotation_report_summary(positive)
stopifnot(internal_summary$annotation_pairs[internal_summary$phase == "manual_curation"] == 4L,
          internal_summary$ambiguous_features[internal_summary$phase == "manual_curation"] == 1L,
          all(c("qc_samples_removed", "adduct_scored") %in% internal_summary$phase))
stopifnot(identical(s$phase[1:3], c("preprocessed_features", "rank1_mz", "rank1_mz_rt")),
          s$features[2] == 5L, s$features[3] == 4L, s$removed_features[3] == 1L,
          "rank1_mz" %in% names(tables))
comparison <- tables$Curated_reference_comparison
stopifnot("adduct_scored" %in% comparison$phase, "adduct_scored" %in% names(tables))
first <- comparison[comparison$phase == "rank1_mz_rt", ][1, ]
stopifnot(identical(names(comparison), c("polarity", "phase", "curated_matches",
          "uncurated_assignments", "missed_curated_pairs", "precision", "recall", "f1")),
          first$curated_matches == 1L,
          first$uncurated_assignments == 2L, first$missed_curated_pairs == 0L,
          first$precision == 0.25, first$recall == 1,
          abs(first$f1 - 0.4) < 1e-12,
          all(comparison$curated_matches[comparison$phase == "manual_curation"] == 1L),
          all(comparison$precision[comparison$phase == "manual_curation"] == 1/3),
          all(comparison$recall[comparison$phase == "manual_curation"] == 1))
qc_without_standards <- positive
qc_without_standards$phases$qc_rsd_filtered$data <-
  qc_without_standards$phases$qc_rsd_filtered$data[!positive$phases$qc_rsd_filtered$data$.is_standard, ]
stopifnot(identical(annotation_reference_comparison(qc_without_standards),
                    annotation_reference_comparison(positive)))
# Multiple matching semicolon assignments and duplicate rows count as one
# feature. A feature with only a nonmatching assignment does not count.
ambiguous_report <- new_annotation_report("TEST", "positive")
ambiguous_report <- capture_annotation_phase(ambiguous_report, "rank1_mz_rt",
  make_input(c("FT1", "FT1", "FT2", "FT3"), c("A_1; B_1; Z_1", "A_1", "D_1", "E_1")))
ambiguous_report <- capture_annotation_phase(ambiguous_report, "manual_curation",
  make_input(c("FT1", "FT2"), c("A_1; B_1", "C_1")))
ambiguous_report <- capture_annotation_phase(ambiguous_report, "qc_rsd_filtered",
  make_input(c("FT1", "FT2"), c("A_1; B_1", "C_1")))
ambiguous_report <- capture_annotation_phase(ambiguous_report, "within_mode_adduct_resolution",
  make_input("FT1", "A_1; B_1"))
ambiguous_comparison <- annotation_reference_comparison(ambiguous_report)
stopifnot(identical(ambiguous_comparison$curated_matches, c(1L, 2L)),
          ambiguous_comparison$precision[1] == 0.4,
          abs(ambiguous_comparison$recall[1] - 2/3) < 1e-12,
          abs(ambiguous_comparison$f1[1] - 0.5) < 1e-12)
empty <- make_input(character(), character())
empty_report <- new_annotation_report("TEST", "positive")
empty_report <- capture_annotation_phase(empty_report, "rank1_mz_rt", empty)
empty_report <- capture_annotation_phase(empty_report, "manual_curation", empty)
empty_report <- capture_annotation_phase(empty_report, "qc_rsd_filtered", empty)
stopifnot(all(annotation_report_summary(empty_report)$features == 0L),
          all(is.na(annotation_reference_comparison(empty_report)$precision)))
missing_qc <- positive
missing_qc$phases$qc_rsd_filtered <- NULL
stopifnot(inherits(try(annotation_reference_comparison(missing_qc), silent = TRUE), "try-error"))
payload <- annotation_metrics_payload(tables)
stopifnot(length(payload$sheets) == 13L,
          grepl("Reference snapshot: qc_rsd_filtered", payload$sheets[[2]]$note, fixed = TRUE),
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

# Current annotations resolve adduct duplicates during ambiguity handling;
# metrics regeneration must accept snapshots without the retired later stage.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  temporary_root <- tempfile()
  folder <- file.path(temporary_root, "applications", "TEST_study")
  dir.create(folder, recursive = TRUE)
  dir.create(file.path(temporary_root, "R"))
  file.copy("R/lipid_helpers.R", file.path(temporary_root, "R", "lipid_helpers.R"))
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
  cat("PASS: metrics regeneration accepts current stages and rejects missing QC snapshots.\n")
}
