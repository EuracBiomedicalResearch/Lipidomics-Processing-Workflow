# Validation script for R translation of impute_block_knn
# This script compares R implementation outputs against MATLAB reference data
#
# Prerequisites:
# 1. Use MATLAB to generate reference data

# Clear workspace
rm(list = ls())

# Load required libraries
if (!require("here")) {
    install.packages("here")
}
library(here)

# SummarizedExperiment needed for SE integration tests
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    library(SummarizedExperiment)
}

source("R/imputers.R")

cat("===== R Translation Validation Suite =====\n\n")
t_total_start <- proc.time()

# Initialize results tracking
validation_results <- data.frame(
    test_num = integer(),
    test_id = character(),
    status = character(),
    max_abs_diff = numeric(),
    rmse = numeric(),
    correlation = numeric(),
    elapsed_sec = numeric(),
    notes = character(),
    stringsAsFactors = FALSE
)

# Helper function to load test case
load_test_case <- function(test_num, test_id) {
    prefix <- sprintf("tests/reference_data/test_%03d_%s", test_num, test_id)

    # Check if test exists
    if (!file.exists(paste0(prefix, "_input.csv"))) {
        return(list(exists = FALSE))
    }

    # Load input
    input_data <- as.matrix(read.csv(
        paste0(prefix, "_input.csv"),
        header = FALSE
    ))

    # Load parameters
    params_file <- paste0(prefix, "_params.txt")
    params_lines <- readLines(params_file)
    a <- as.numeric(sub("a = ", "", params_lines[1]))
    b <- as.numeric(sub("b = ", "", params_lines[2]))
    k <- as.integer(sub("k = ", "", params_lines[3]))

    # Load MATLAB output (if exists)
    matlab_output <- NULL
    error_msg <- NULL

    if (file.exists(paste0(prefix, "_output.csv"))) {
        matlab_output <- as.matrix(read.csv(
            paste0(prefix, "_output.csv"),
            header = FALSE
        ))
    } else if (file.exists(paste0(prefix, "_error.txt"))) {
        error_lines <- readLines(paste0(prefix, "_error.txt"))
        error_msg <- paste(error_lines, collapse = "\n")
    }

    return(list(
        exists = TRUE,
        test_num = test_num,
        test_id = test_id,
        input = input_data,
        a = a,
        b = b,
        k = k,
        matlab_output = matlab_output,
        matlab_error = error_msg
    ))
}

# Helper function to compare matrices
compare_outputs <- function(matlab_output, r_output, tolerance = 6e-5) {
    if (is.null(matlab_output) || is.null(r_output)) {
        return(list(
            match = FALSE,
            max_abs_diff = NA,
            rmse = NA,
            correlation = NA,
            notes = "One or both outputs are NULL"
        ))
    }

    # Check dimensions
    if (!all(dim(matlab_output) == dim(r_output))) {
        return(list(
            match = FALSE,
            max_abs_diff = NA,
            rmse = NA,
            correlation = NA,
            notes = sprintf(
                "Dimension mismatch: MATLAB %dx%d vs R %dx%d",
                nrow(matlab_output),
                ncol(matlab_output),
                nrow(r_output),
                ncol(r_output)
            )
        ))
    }

    # Calculate differences
    diff_matrix <- matlab_output - r_output
    max_abs_diff <- max(abs(diff_matrix))
    rmse <- sqrt(mean(diff_matrix^2))

    # Calculate correlation (only for non-constant matrices)
    correlation <- NA
    if (sd(as.vector(matlab_output)) > 0 && sd(as.vector(r_output)) > 0) {
        correlation <- cor(as.vector(matlab_output), as.vector(r_output))
    }

    # Determine if match
    match <- max_abs_diff < tolerance

    notes <- if (match) {
        sprintf("Perfect match (max diff: %.2e)", max_abs_diff)
    } else {
        sprintf("Differences found (max diff: %.2e)", max_abs_diff)
    }

    return(list(
        match = match,
        max_abs_diff = max_abs_diff,
        rmse = rmse,
        correlation = correlation,
        notes = notes
    ))
}

# Helper function to run single validation test
validate_test <- function(test_num, test_id) {
    cat(sprintf("Test %d (%s): ", test_num, test_id))

    # Load test case
    test_case <- load_test_case(test_num, test_id)

    if (!test_case$exists) {
        cat("SKIPPED (not found)\n")
        return(data.frame(
            test_num = test_num,
            test_id = test_id,
            status = "SKIPPED",
            max_abs_diff = NA,
            rmse = NA,
            correlation = NA,
            elapsed_sec = NA,
            notes = "Test file not found"
        ))
    }

    # Check if MATLAB had an error
    if (!is.null(test_case$matlab_error)) {
        # R should also produce an error or handle gracefully
        cat("MATLAB_ERROR - validating R handles same case\n")

        r_error <- NULL
        r_output <- tryCatch(
            {
                impute_block_knn(
                    test_case$input,
                    test_case$a,
                    test_case$b,
                    test_case$k
                )
            },
            error = function(e) {
                r_error <<- e$message
                NULL
            }
        )

        status <- if (!is.null(r_error)) "MATCH_ERROR" else "R_NO_ERROR"
        notes <- sprintf(
            "MATLAB error: %s; R error: %s",
            substr(test_case$matlab_error, 1, 50),
            if (!is.null(r_error)) substr(r_error, 1, 50) else "None"
        )

        cat(sprintf("  %s\n", status))

        return(data.frame(
            test_num = test_num,
            test_id = test_id,
            status = status,
            max_abs_diff = NA,
            rmse = NA,
            correlation = NA,
            elapsed_sec = NA,
            notes = notes
        ))
    }

    # Run R implementation
    r_output <- NULL
    r_error <- NULL
    elapsed <- NA

    t0 <- proc.time()
    r_output <- tryCatch(
        {
            impute_block_knn(
                test_case$input,
                test_case$a,
                test_case$b,
                test_case$k
            )
        },
        error = function(e) {
            r_error <<- e$message
            NULL
        }
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]

    if (!is.null(r_error)) {
        cat(sprintf("R_ERROR: %s\n", substr(r_error, 1, 60)))
        return(data.frame(
            test_num = test_num,
            test_id = test_id,
            status = "R_ERROR",
            max_abs_diff = NA,
            rmse = NA,
            correlation = NA,
            elapsed_sec = elapsed,
            notes = r_error
        ))
    }

    # Compare outputs
    comparison <- compare_outputs(test_case$matlab_output, r_output)

    status <- if (comparison$match) "PASS" else "FAIL"
    cat(sprintf("%s (%.2fs) - %s\n", status, elapsed, comparison$notes))

    return(data.frame(
        test_num = test_num,
        test_id = test_id,
        status = status,
        max_abs_diff = comparison$max_abs_diff,
        rmse = comparison$rmse,
        correlation = comparison$correlation,
        elapsed_sec = elapsed,
        notes = comparison$notes
    ))
}

# Main validation loop
cat("Starting validation...\n\n")

# Run all tests
# The test numbers and IDs must match those in test_knn_imputation.m

test_specs <- list(
    # Input validation tests
    list(1, "VAL-01"),
    list(2, "VAL-02"),
    list(3, "VAL-03"),
    list(4, "VAL-04"),
    list(5, "VAL-05"),
    list(6, "VAL-06"),

    # Block assignment tests
    list(7, "BLK-01"),
    list(8, "BLK-02"),
    list(9, "BLK-03"),
    list(10, "BLK-04"),
    list(11, "BLK-05"),
    list(12, "BLK-06"),

    # Functional tests
    list(13, "B1-01"),
    list(14, "B1-02-k1"),
    list(15, "B1-02-k3"),
    list(16, "B2-01"),
    list(17, "B2-02"),
    list(18, "B3-01"),
    list(19, "B3-02"),

    # Edge cases
    list(20, "EDG-01"),
    list(21, "EDG-02"),
    list(22, "EDG-03"),
    list(23, "EDG-04"),
    list(24, "EDG-05"),
    list(25, "EDG-06"),

    # Comprehensive test cases
    list(26, "TC-01"),
    list(27, "TC-02"),
    list(28, "TC-03"),
    list(29, "TC-04"),

    # fallback tests
    list(34, "FALLBACK-B2-basic"),
    list(35, "FALLBACK-B2-multi"),
    list(36, "FALLBACK-B2-mixed")
)

for (spec in test_specs) {
    result <- validate_test(spec[[1]], spec[[2]])
    validation_results <- rbind(validation_results, result)
}

# ============================================================================
# SummarizedExperiment integration tests
# ============================================================================
cat("\n--- SummarizedExperiment Tests ---\n")

se_test_specs <- list(
    list(32, "POS-SE-filtered"),
    list(33, "NEG-SE-filtered"),
    list(30, "POS-SE"),
    list(31, "NEG-SE"),
    list(37, "POS-SE-filtered_1"),
    list(38, "POS-SE_1")
)

for (spec in se_test_specs) {
    test_num <- spec[[1]]
    test_id <- spec[[2]]
    prefix <- sprintf("tests/reference_data/test_%03d_%s", test_num, test_id)

    cat(sprintf("Test %d (%s): ", test_num, test_id))

    rds_path <- paste0(prefix, "_input.RDS")
    output_path <- paste0(prefix, "_output.csv")
    params_path <- paste0(prefix, "_params.txt")

    # Skip if RDS not found
    if (!file.exists(rds_path)) {
        cat("SKIPPED (RDS not found)\n")
        validation_results <- rbind(
            validation_results,
            data.frame(
                test_num = test_num,
                test_id = test_id,
                status = "SKIPPED",
                max_abs_diff = NA,
                rmse = NA,
                correlation = NA,
                elapsed_sec = NA,
                notes = "RDS file not found - run create_se_test_data.R first"
            )
        )
        next
    }

    # Load SE and params
    se_input <- readRDS(rds_path)
    params_lines <- readLines(params_path)
    a <- as.numeric(sub("a = ", "", params_lines[1]))
    b <- as.numeric(sub("b = ", "", params_lines[2]))
    k <- as.integer(sub("k = ", "", params_lines[3]))

    # Run impute_block_knn on the SE
    r_error <- NULL
    elapsed <- NA

    t0 <- proc.time()
    se_result <- tryCatch(
        {
            impute_block_knn(se_input, a, b, k)
        },
        error = function(e) {
            r_error <<- e$message
            NULL
        }
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]

    if (!is.null(r_error)) {
        cat(sprintf("R_ERROR: %s\n", substr(r_error, 1, 60)))
        validation_results <- rbind(
            validation_results,
            data.frame(
                test_num = test_num,
                test_id = test_id,
                status = "R_ERROR",
                max_abs_diff = NA,
                rmse = NA,
                correlation = NA,
                elapsed_sec = elapsed,
                notes = r_error
            )
        )
        next
    }

    # Verify return type is SummarizedExperiment
    if (!is(se_result, "SummarizedExperiment")) {
        cat("FAIL - return type is not SummarizedExperiment\n")
        validation_results <- rbind(
            validation_results,
            data.frame(
                test_num = test_num,
                test_id = test_id,
                status = "FAIL",
                max_abs_diff = NA,
                rmse = NA,
                correlation = NA,
                elapsed_sec = elapsed,
                notes = sprintf("Wrong return type: %s", class(se_result)[1])
            )
        )
        next
    }

    # Check no NAs remain in imputed assay
    r_output <- SummarizedExperiment::assay(se_result)
    remaining_na <- sum(is.na(r_output))
    if (remaining_na > 0) {
        cat(sprintf("WARN - %d NAs remain after imputation\n", remaining_na))
    }

    # Compare against MATLAB output if available
    if (file.exists(output_path)) {
        matlab_output <- as.matrix(read.csv(output_path, header = FALSE))
        comparison <- compare_outputs(matlab_output, r_output)
        status <- if (comparison$match) "PASS" else "FAIL"
        cat(sprintf("%s (%.2fs) - %s\n", status, elapsed, comparison$notes))
        validation_results <- rbind(
            validation_results,
            data.frame(
                test_num = test_num,
                test_id = test_id,
                status = status,
                max_abs_diff = comparison$max_abs_diff,
                rmse = comparison$rmse,
                correlation = comparison$correlation,
                elapsed_sec = elapsed,
                notes = comparison$notes
            )
        )
    } else {
        cat(sprintf(
            "PASS (%.2fs) - SE dispatch OK, no MATLAB output to compare\n",
            elapsed
        ))
        validation_results <- rbind(
            validation_results,
            data.frame(
                test_num = test_num,
                test_id = test_id,
                status = "PASS",
                max_abs_diff = NA,
                rmse = NA,
                correlation = NA,
                elapsed_sec = elapsed,
                notes = "SE dispatch verified; MATLAB output not yet available"
            )
        )
    }
}

cat("\n===== VALIDATION SUMMARY =====\n")

# Summary statistics
status_counts <- table(validation_results$status)
print(status_counts)

cat(sprintf("\nTotal tests: %d\n", nrow(validation_results)))
cat(sprintf("Passed: %d\n", sum(validation_results$status == "PASS")))
cat(sprintf("Failed: %d\n", sum(validation_results$status == "FAIL")))
cat(sprintf(
    "Errors: %d\n",
    sum(validation_results$status %in% c("R_ERROR", "MATLAB_ERROR"))
))
cat(sprintf("Skipped: %d\n", sum(validation_results$status == "SKIPPED")))

# Save results
write.csv(validation_results, "validation_results.csv", row.names = FALSE)
cat("\nDetailed results saved to: validation_results.csv\n")

# Print failures if any
failures <- validation_results[validation_results$status == "FAIL", ]
if (nrow(failures) > 0) {
    cat("\n--- FAILED TESTS ---\n")
    print(failures[, c("test_id", "max_abs_diff", "rmse", "notes")])
}

# Print worst numerical differences (for passed tests)
passed <- validation_results[
    validation_results$status == "PASS" & !is.na(validation_results$max_abs_diff),
]
if (nrow(passed) > 0) {
    cat("\n--- WORST NUMERICAL DIFFERENCES (PASSED TESTS) ---\n")
    worst <- passed[order(-passed$max_abs_diff), ][1:min(5, nrow(passed)), ]
    print(worst[, c("test_id", "max_abs_diff", "rmse", "correlation")])
}

cat("\n===== VALIDATION COMPLETE =====\n")
