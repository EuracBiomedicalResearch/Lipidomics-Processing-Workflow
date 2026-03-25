library(testthat)

# =============================================================================
# Regression tests: compare dev snapshots against committed golden baselines.
#
# Workflow:
#   1. Run Lipidomics_workflow.qmd with SNAPSHOT_DIR = "tests/snapshots_dev"
#      (for both POLARITY values and both SMOKE modes as needed)
#   2. Run this script to verify outputs match the golden baselines in
#      tests/snapshots/
#
# Usage (from project root):
#   Rscript -e "setwd('/path/to/project'); testthat::test_file('tests/test_regression.R')"
# =============================================================================

GOLDEN_DIR  <- "../tests/snapshots"
DEV_DIR     <- "../tests/snapshots_dev"

# Skip the whole file if the dev directory hasn't been populated yet
if (!dir.exists(DEV_DIR) || length(list.files(DEV_DIR)) == 0L) {
    cat("No dev snapshots found in", DEV_DIR,
        "— run the workflow with SNAPSHOT_DIR = 'tests/snapshots_dev' first.\n")
    testthat::skip("dev snapshot directory is empty")
}

golden_files <- list.files(GOLDEN_DIR, pattern = "\\.rds$", full.names = FALSE)

for (f in golden_files) {
    golden_path <- file.path(GOLDEN_DIR, f)
    dev_path    <- file.path(DEV_DIR,    f)

    test_that(paste("snapshot matches golden:", f), {
        skip_if_not(file.exists(dev_path),
                    paste("dev snapshot missing:", f,
                          "— re-run workflow for this polarity/mode"))

        golden <- readRDS(golden_path)
        dev    <- readRDS(dev_path)

        expect_equal(dev, golden)
    })
}
