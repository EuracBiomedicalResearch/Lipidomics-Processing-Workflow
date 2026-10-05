# Run from the repository root: Rscript tests/test_rank1_mz_reporting.R
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MetaboAnnotation))
source("R/lipid_helpers.R")
source("R/annotation_reporting.R")
res <- SummarizedExperiment(assays = list(raw = matrix(1, 3, 1,
  dimnames = list(c("FT1", "FT2", "FT3"), "sample"))))
rowData(res) <- S4Vectors::DataFrame(mzmed = c(100, 100, 200), rtmed = c(10, 200, 10))
database <- data.frame(rank = c(1L, 2L), mz = c(100, 200), rt_adjusted = c(10, 10),
  lipid_name_unique = c("A_1", "B_1"), lipid.name = c("A", "B"))
mz <- match_features_mz_only(res, database, ppm = 20)
rt <- match_features_to_database(res, database, ppm = 20, rt_tol = 20, verbose = FALSE)$mtched_data
stopifnot(identical(mz$feature_id, c("FT1", "FT2")),
          identical(rt$feature_id, "FT1"), !"FT3" %in% mz$feature_id)
# Removing RT columns must have no effect on m/z-only candidates.
no_rt <- database[, setdiff(names(database), "rt_adjusted")]
stopifnot(identical(mz[, c("feature_id", "target_lipid_name_unique")],
  match_features_mz_only(res, no_rt)[, c("feature_id", "target_lipid_name_unique")]))
report <- new_annotation_report("TEST", "positive")
report <- capture_annotation_phase(report, "preprocessed_features", res, annotations = FALSE)
report <- capture_annotation_phase(report, "rank1_mz_rt", rt)
before <- report$phases
report <- insert_rank1_mz_phase(report, mz)
stopifnot(identical(report$phases[names(before)], before),
          identical(names(report$phases), c("preprocessed_features", "rank1_mz", "rank1_mz_rt")))
summary <- annotation_report_summary(report)
stopifnot(identical(summary$features, c(3L, 2L, 1L)), summary$removed_features[3] == 1L)
bad <- mz[mz$feature_id != "FT1", ]
old <- new_annotation_report("TEST", "positive"); old$phases <- before
stopifnot(inherits(try(insert_rank1_mz_phase(old, bad), silent = TRUE), "try-error"))
empty <- match_features_mz_only(res, transform(database, mz = mz + 1000))
stopifnot(nrow(empty) == 0L)
cat("PASS: mass-only ignores RT, excludes rank-2 entries, preserves saved phases, and checks candidate containment\n")
