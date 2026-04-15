#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# test_annotation.sh
#
# Render applications/pilot_study/Annotation.qmd (which writes fresh snapshots
# to tests/snapshots_annotation_dev/) then run the regression test suite to
# compare against the golden baselines in tests/snapshots_annotation/.
#
# Usage (from project root):
#   bash test_annotation.sh
#
# Exit codes:
#   0 — render succeeded and all regression tests passed
#   1 — render failed (tests are skipped)
#   2 — render succeeded but one or more regression tests failed
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QMD="$SCRIPT_DIR/applications/pilot_study/Annotation.qmd"
TEST="$SCRIPT_DIR/tests/test_regression_annotation.R"

# Locate quarto — RStudio bundles its own binary outside the system PATH
QUARTO_BIN=$(command -v quarto 2>/dev/null \
    || echo "/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto")
if [ ! -x "$QUARTO_BIN" ]; then
    echo "ERROR: quarto binary not found. Set PATH or install quarto." >&2
    exit 1
fi
echo "Using quarto: $QUARTO_BIN"

echo "=== Step 1: Render Annotation.qmd ==="
if ! "$QUARTO_BIN" render "$QMD" --to html; then
    echo "ERROR: Quarto render failed. Regression tests skipped." >&2
    exit 1
fi

echo ""
echo "=== Step 2: Run regression tests ==="
Rscript -e "
setwd('$SCRIPT_DIR')
results <- testthat::test_file('$TEST')
if (any(as.data.frame(results)[['failed']] > 0)) quit(status = 2)
"
