# Run from the repository root: Rscript tests/test_annotation_istd_exclusion.R
source("R/annotation_reporting.R")
frame <- function(ids, labels, standard = rep(FALSE, length(ids))) {
  data.frame(feature_id = ids, polarity = "positive", annotation = labels,
             .is_standard = standard)
}
report <- new_annotation_report("TEST", "positive")
report$phases$preprocessed_features <- list(
  data = frame(c("IS1", "IS_FAILED_QC", "FT1"), rep(NA_character_, 3)), annotations = FALSE)
report$phases$rank1_mz <- list(data = frame(
  c("IS1", "IS1", "IS_FAILED_QC", "FT1", "FT1"),
  c("standard", "alternative", "standard", "endogenous", "standard"),
  c(TRUE, FALSE, TRUE, FALSE, TRUE)), annotations = TRUE)
report$phases$manual_curation <- list(data = frame(
  c("IS1", "IS_FAILED_QC", "FT1"), c("standard", "standard", "endogenous"),
  c(TRUE, TRUE, FALSE)), annotations = TRUE)
report$phases$qc_rsd_filtered <- list(data = frame(
  c("IS1", "FT1"), c("standard", "endogenous"), c(TRUE, FALSE)), annotations = TRUE)
original <- report
filtered <- exclude_annotation_istds(report)
stopifnot(identical(report, original),
          all(vapply(filtered$phases, function(s) identical(s$data$feature_id, "FT1"), logical(1))),
          identical(filtered$phases$rank1_mz$data$annotation, "endogenous"),
          all(annotation_report_summary(report)$features == 1L))
# Identical feature IDs in the other polarity are not excluded accidentally.
other <- new_annotation_report("TEST", "negative")
other$phases$preprocessed_features <- list(data = transform(
  frame("IS1", NA_character_), polarity = "negative"), annotations = FALSE)
stopifnot(nrow(exclude_annotation_istds(other, annotation_istd_keys(list(report)))$
  phases$preprocessed_features$data) == 1L)
cat("PASS: ISTDs excluded from every phase, alternative labels, preprocessing and failed-QC rows; endogenous candidates and polarity IDs preserved.\n")
