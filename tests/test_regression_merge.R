library(testthat)

# =============================================================================
# Regression tests: compare merge dev snapshots against golden baselines.
#
# Workflow:
#   1. Run applications/pilot_study/POS_NEG_merge.qmd once with
#      SNAPSHOT_DIR = "tests/snapshots_merge_dev" to produce the baseline.
#      Promote to golden:
#        dir.create("tests/snapshots_merge", showWarnings = FALSE)
#        file.copy(
#            list.files("tests/snapshots_merge_dev", full.names = TRUE),
#            "tests/snapshots_merge/"
#        )
#   2. After each refactor step, re-run POS_NEG_merge.qmd (same SNAPSHOT_DIR).
#   3. Run this script to verify outputs match the golden baselines:
#        Rscript -e "setwd('/path/to/project')
#                    testthat::test_file('tests/test_regression_merge.R')"
#
# Snapshot checkpoints (polarity "both" — each captures data from both modes):
#   merge_aligned    — after sample alignment (pos vs neg colnames)
#   merge_duplicates — after duplicate resolution (comparison_df)
#   merge_combined   — after rbind (ionization_mode + lipid identifiers)
# =============================================================================

GOLDEN_DIR <- "../tests/snapshots_merge"
DEV_DIR    <- "../tests/snapshots_merge_dev"

if (!dir.exists(DEV_DIR) || length(list.files(DEV_DIR)) == 0L) {
    cat("No dev snapshots found in", DEV_DIR,
        "— run POS_NEG_merge.qmd with SNAPSHOT_DIR pointing there first.\n")
    testthat::skip("dev snapshot directory is empty")
}

golden_files <- list.files(GOLDEN_DIR, pattern = "\\.rds$", full.names = FALSE)

for (f in golden_files) {
    golden_path <- file.path(GOLDEN_DIR, f)
    dev_path    <- file.path(DEV_DIR,    f)

    test_that(paste("merge snapshot matches golden:", f), {
        skip_if_not(file.exists(dev_path),
                    paste("dev snapshot missing:", f,
                          "— re-run POS_NEG_merge.qmd"))

        golden <- readRDS(golden_path)
        dev    <- readRDS(dev_path)

        expect_equal(dev, golden)
    })
}
