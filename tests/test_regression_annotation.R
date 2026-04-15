library(testthat)

# =============================================================================
# Regression tests: compare annotation dev snapshots against golden baselines.
#
# Workflow:
#   1. Run applications/pilot_study/Annotation.qmd once with
#      SNAPSHOT_DIR = "tests/snapshots_annotation_dev" to produce the baseline.
#      Promote to golden:
#        file.copy(
#            list.files("tests/snapshots_annotation_dev", full.names = TRUE),
#            "tests/snapshots_annotation/"
#        )
#   2. After each refactor step, re-run Annotation.qmd (same SNAPSHOT_DIR).
#   3. Run this script to verify outputs match the golden baselines:
#        Rscript -e "setwd('/path/to/project')
#                    testthat::test_file('tests/test_regression_annotation.R')"
#
# Snapshot checkpoints (both "pos" and "neg" per step):
#   annot_rank1    — after match_features_to_database
#   annot_isotope  — after calculate_isotope_similarity
#   annot_adducts  — after match_adducts
#   annot_resolved — after resolve_annotation_ambiguities + resolve_sm_isomers
#   annot_isnorm   — after normalize_by_is
#   annot_final    — after RSD filter + QC removal
# =============================================================================

GOLDEN_DIR <- "../tests/snapshots_annotation"
DEV_DIR    <- "../tests/snapshots_annotation_dev"

if (!dir.exists(DEV_DIR) || length(list.files(DEV_DIR)) == 0L) {
    cat("No dev snapshots found in", DEV_DIR,
        "— run Annotation.qmd with SNAPSHOT_DIR pointing there first.\n")
    testthat::skip("dev snapshot directory is empty")
}

golden_files <- list.files(GOLDEN_DIR, pattern = "\\.rds$", full.names = FALSE)

for (f in golden_files) {
    golden_path <- file.path(GOLDEN_DIR, f)
    dev_path    <- file.path(DEV_DIR,    f)

    test_that(paste("annotation snapshot matches golden:", f), {
        skip_if_not(file.exists(dev_path),
                    paste("dev snapshot missing:", f,
                          "— re-run Annotation.qmd for this polarity"))

        golden <- readRDS(golden_path)
        dev    <- readRDS(dev_path)

        expect_equal(dev, golden)
    })
}
