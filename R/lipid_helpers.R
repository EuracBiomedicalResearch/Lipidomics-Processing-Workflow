#' ===========================================================================
#' Lipidomics Workflow Helper Functions
#' ===========================================================================
#'
#' This file contains reusable functions for the lipidomics data processing

#' and annotation workflow. Source this file at the beginning of your analysis.
#'
#' @author Philippine Louail, Sara Londono-Osorio
#' @date 2026
#' ===========================================================================

# =============================================================================
# VALIDATION / FAILSAFE FUNCTIONS
# =============================================================================

#' Validate the sample metadata Excel file
#'
#' Checks that the metadata file exists, is readable, and contains all
#' required columns with sensible values.
#'
#' @param file_path Path to the sample sequence Excel file
#' @param required_cols Character vector of required column names
#' @param expected_sample_types Character vector of expected sample_type values
#'   (optional; warns if unexpected types found)
#' @param data_dir Directory where mzML files should be (optional; checks files
#'   exist)
#' @return The validated data frame (invisibly). Stops on critical errors,
#'   warns on non-critical issues.
#' @examples
#' validate_metadata("seq_pos.xlsx", data_dir = "data")
validate_metadata <- function(file_path,
                              required_cols = c("file_name", "sample_name",
                                                "sample_type",
                                                "injection_index"),
                              expected_sample_types = NULL,
                              data_dir = NULL) {

  # --- File existence & readability ---
  if (!file.exists(file_path)) {
    stop("Metadata file not found: '", file_path,
         "'\n  Make sure the file is in the current working directory.")
  }

  df <- tryCatch(
    readxl::read_xlsx(file_path, col_names = TRUE) |> as.data.frame(),
    error = function(e) {
      stop("Cannot read metadata file '", file_path, "': ", conditionMessage(e))
    }
  )

  if (nrow(df) == 0) stop("Metadata file '", file_path, "' is empty (0 rows).")

  # --- Required columns ---
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Metadata file is missing required column(s): ",
         paste0("'", missing, "'", collapse = ", "),
         "\n  Found columns: ",
         paste0("'", colnames(df), "'", collapse = ", "))
  }

  # --- Duplicate file names ---
  dup_files <- df$file_name[duplicated(df$file_name)]
  if (length(dup_files) > 0) {
    stop("Duplicate file_name entries in metadata: ",
         paste0("'", unique(dup_files), "'", collapse = ", "))
  }

  # --- Duplicate sample names ---
  dup_samples <- df$sample_name[duplicated(df$sample_name)]
  if (length(dup_samples) > 0) {
    warning("Duplicate sample_name entries: ",
            paste0("'", unique(dup_samples), "'", collapse = ", "),
            "\n  This may cause issues with downstream analyses.")
  }

  # --- NA values in critical columns ---
  for (col in required_cols) {
    n_na <- sum(is.na(df[[col]]) | trimws(df[[col]]) == "")
    if (n_na > 0) {
      stop("Column '", col, "' has ", n_na, " missing/empty value(s).")
    }
  }

  # --- Injection index must be numeric and unique ---
  if ("injection_index" %in% colnames(df)) {
    idx <- suppressWarnings(as.numeric(df$injection_index))
    if (any(is.na(idx))) {
      stop("'injection_index' contains non-numeric values: ",
           paste0("'", df$injection_index[is.na(idx)], "'", collapse = ", "))
    }
    if (any(duplicated(idx))) {
      warning("'injection_index' has duplicate values: ",
              paste(idx[duplicated(idx)], collapse = ", "))
    }
  }

  # --- Sample type check ---
  actual_types <- unique(df$sample_type)
  if (!is.null(expected_sample_types)) {
    unexpected <- setdiff(actual_types, expected_sample_types)
    if (length(unexpected) > 0) {
      warning("Unexpected sample_type values: ",
              paste0("'", unexpected, "'", collapse = ", "),
              "\n  Expected: ",
              paste0("'", expected_sample_types, "'", collapse = ", "))
    }
  }

  # --- Must have at least one QC ---
  if (!"QC" %in% actual_types) {
    warning("No 'QC' samples found in metadata. ",
            "QC-based normalization and RSD filtering will fail.")
  }

  # --- Check mzML files exist ---
  if (!is.null(data_dir)) {
    paths <- file.path(data_dir, df$file_name)
    missing_files <- paths[!file.exists(paths)]
    if (length(missing_files) > 0) {
      stop(length(missing_files), " mzML file(s) listed in metadata not found:\n  ",
           paste(basename(head(missing_files, 5)), collapse = "\n  "),
           if (length(missing_files) > 5)
             paste0("\n  ... and ", length(missing_files) - 5, " more"))
    }
  }

  message("\u2713 Metadata validated: ", nrow(df), " samples, ",
          length(actual_types), " sample types (",
          paste(actual_types, collapse = ", "), ")")
  invisible(df)
}


#' Validate the lipid reference set (internal standard) Excel file
#'
#' Checks the file exists and contains the expected columns with valid values.
#'
#' @param file_path Path to the reference lipid Excel file
#' @param required_cols Required columns (default: short_name, mz, RT)
#' @return The validated data frame (invisibly)
validate_reference_lipids <- function(file_path,
                                      required_cols = c("short_name",
                                                        "mz", "RT")) {
  if (!file.exists(file_path)) {
    stop("Reference lipid file not found: '", file_path, "'")
  }

  df <- tryCatch(
    readxl::read_xlsx(file_path) |> as.data.frame(),
    error = function(e) {
      stop("Cannot read reference lipid file '", file_path, "': ",
           conditionMessage(e))
    }
  )

  if (nrow(df) == 0) stop("Reference lipid file is empty.")

  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Reference lipid file missing column(s): ",
         paste0("'", missing, "'", collapse = ", "),
         "\n  Found: ", paste0("'", colnames(df), "'", collapse = ", "))
  }

  # Numeric checks
  if (!is.numeric(df$mz) || any(is.na(df$mz)) || any(df$mz <= 0)) {
    stop("Column 'mz' must contain positive numeric values with no NAs.")
  }
  if (!is.numeric(df$RT) || any(is.na(df$RT)) || any(df$RT <= 0)) {
    stop("Column 'RT' must contain positive numeric values with no NAs.")
  }

  # Duplicate short_name
  if (any(duplicated(df$short_name))) {
    stop("Duplicate short_name entries: ",
         paste(df$short_name[duplicated(df$short_name)], collapse = ", "))
  }

  message("\u2713 Reference lipids validated: ", nrow(df), " compounds (",
          paste(df$short_name, collapse = ", "), ")")
  invisible(df)
}


#' Validate the lipid database Excel file
#'
#' Checks the database file and specific sheet for expected structure.
#'
#' @param db_path Path to the SRM 1950-derived lipid database Excel file
#' @param sheet Sheet number to validate
#' @param polarity "pos" or "neg"
#' @param rt_col Name of the retention time column
#' @param required_cols Additional required column names
#' @return The raw data frame (invisibly)
validate_lipid_database <- function(db_path,
                                    sheet,
                                    polarity = c("pos", "neg"),
                                    rt_col,
                                    required_cols = c("lipid name",
                                                      "MOLECULAR FORMULA",
                                                      "mz", "Adduct",
                                                      "rank")) {
  polarity <- match.arg(polarity)

  if (!file.exists(db_path)) {
    stop("Lipid database file not found: '", db_path, "'")
  }

  # Check the sheet exists
  available_sheets <- readxl::excel_sheets(db_path)
  if (sheet > length(available_sheets)) {
    stop("Sheet ", sheet, " does not exist in '", db_path,
         "'. Available sheets (", length(available_sheets), "): ",
         paste0("'", available_sheets, "'", collapse = ", "))
  }

  df <- tryCatch(
    readxl::read_xlsx(db_path, sheet = sheet),
    error = function(e) {
      stop("Cannot read sheet ", sheet, " from '", db_path, "': ",
           conditionMessage(e))
    }
  )

  if (nrow(df) == 0) {
    stop("Lipid database sheet ", sheet, " is empty.")
  }

  # Check RT column
  all_required <- c(required_cols, rt_col)
  missing <- setdiff(all_required, colnames(df))
  if (length(missing) > 0) {
    stop("Lipid database (sheet ", sheet, ") missing column(s): ",
         paste0("'", missing, "'", collapse = ", "),
         "\n  Found: ",
         paste0("'", head(colnames(df), 15), "'", collapse = ", "),
         if (ncol(df) > 15) "...")
  }

  # Rank column should contain numeric values including 0 and 1
  ranks <- df$rank
  if (!is.numeric(ranks)) {
    stop("Column 'rank' must be numeric (found ", class(ranks), ").")
  }
  if (!any(ranks == 1, na.rm = TRUE)) {
    warning("No rank=1 entries found in the database. ",
            "Rank 1 matching will return no results.")
  }

  # Check mz is numeric
  mz_vals <- suppressWarnings(as.numeric(df$mz))
  n_na_mz <- sum(is.na(mz_vals) & !is.na(df$mz))
  if (n_na_mz > 0) {
    warning(n_na_mz, " non-numeric value(s) in 'mz' column.")
  }

  n_lipids <- length(unique(df$`lipid name`[!is.na(df$`lipid name`)]))
  message("\u2713 Lipid database validated (sheet ", sheet, "): ",
          nrow(df), " entries, ", n_lipids, " unique lipids")
  invisible(df)
}


#' Validate configuration parameters
#'
#' Checks that the user-defined configuration parameters are within
#' reasonable ranges and consistent with each other.
#'
#' @param polarity "pos" or "neg" ionization mode
#' @param data_source "local", "sqlite", or "metaboLights"
#' @param cores_nb Number of CPU cores
#' @param rt_filter Range for RT filter c(min, max)
#' @param peak_width Peak width range c(min, max)
#' @param ppm PPM tolerance
#' @param sn_threshold Signal-to-noise threshold
#' @param match_ppm PPM for matching
#' @param match_rt_tol RT tolerance for matching
#' @param isopeak_sim Isotope similarity threshold
#' @param rsd_threshold RSD threshold
#' @return TRUE invisibly. Stops or warns on issues.
validate_config <- function(polarity,
                            data_source,
                            cores_nb,
                            rt_filter = c(10, 800),
                            peak_width = c(4, 8),
                            ppm = 10,
                            sn_threshold = 2,
                            match_ppm = 20,
                            match_rt_tol = 20,
                            isopeak_sim = 0.78,
                            rsd_threshold = 0.3) {

  # Polarity
  if (!polarity %in% c("pos", "neg")) {
    stop("POLARITY must be 'pos' or 'neg', got '", polarity, "'")
  }

  # Data source
  if (!data_source %in% c("local", "sqlite", "metaboLights")) {
    stop("DATA_SOURCE must be 'local', 'sqlite', or 'metaboLights', got '",
         data_source, "'")
  }

  # Cores
  max_cores <- parallel::detectCores()
  if (!is.numeric(cores_nb) || cores_nb < 1) {
    stop("CORES_NB must be a positive integer.")
  }
  if (cores_nb > max_cores) {
    warning("CORES_NB (", cores_nb, ") exceeds available cores (", max_cores,
            "). Setting to ", max_cores - 1, " is recommended.")
  }

  # RT filter range
  if (rt_filter[1] >= rt_filter[2]) {
    stop("RT_FILTER_MIN (", rt_filter[1], ") must be less than RT_FILTER_MAX (",
         rt_filter[2], ").")
  }
  if (rt_filter[1] < 0) stop("RT_FILTER_MIN cannot be negative.")

  # Peak width
  if (peak_width[1] >= peak_width[2]) {
    stop("PEAK_WIDTH min (", peak_width[1], ") must be less than max (",
         peak_width[2], ").")
  }
  if (peak_width[1] <= 0) stop("PEAK_WIDTH values must be positive.")

  # PPM
  if (ppm <= 0 || ppm > 100) {
    warning("PPM value (", ppm, ") is outside typical range (1-50).")
  }

  # SN threshold
  if (sn_threshold < 1) {
    warning("SN_THRESHOLD (", sn_threshold,
            ") is below 1 — very permissive peak detection.")
  }

  # Match PPM
  if (match_ppm <= 0 || match_ppm > 100) {
    warning("MATCH_PPM (", match_ppm, ") is outside typical range (5-50).")
  }

  # Match RT tolerance
  if (match_rt_tol <= 0) stop("MATCH_RT_TOL must be positive.")
  if (match_rt_tol > 60) {
    warning("MATCH_RT_TOL (", match_rt_tol,
            "s) is very large. Typical values: 10-30s.")
  }

  # Isotope similarity
  if (isopeak_sim < 0 || isopeak_sim > 1) {
    stop("ISOPEAK_SIM_THRESHOLD must be between 0 and 1.")
  }
  if (isopeak_sim < 0.5) {
    warning("ISOPEAK_SIM_THRESHOLD (", isopeak_sim,
            ") is low — may accept poor isotope matches.")
  }

  # RSD threshold
  if (rsd_threshold <= 0 || rsd_threshold > 1) {
    stop("RSD_THRESHOLD must be between 0 and 1 (e.g., 0.3 = 30%).")
  }

  message("\u2713 Configuration parameters validated")
  invisible(TRUE)
}


#' Validate a preprocessed MsExperiment object
#'
#' Checks that the MsExperiment has expected properties after preprocessing.
#'
#' @param mse An MsExperiment object
#' @param expected_n_samples Expected number of samples (optional)
#' @param check_peaks Logical, check that chromatographic peaks exist
#' @param min_peaks_per_sample Minimum peaks expected per sample (warns below)
#' @return TRUE invisibly. Stops or warns on issues.
validate_mse <- function(mse,
                         expected_n_samples = NULL,
                         check_peaks = TRUE,
                         min_peaks_per_sample = 100) {

  if (!inherits(mse, "MsExperiment")) {
    stop("Object is not an MsExperiment (got ", class(mse)[1], ")")
  }

  n <- length(mse)
  if (n == 0) stop("MsExperiment is empty (0 samples).")

  if (!is.null(expected_n_samples) && n != expected_n_samples) {
    warning("Expected ", expected_n_samples, " samples but MsExperiment has ",
            n, ".")
  }

  # Check spectra
  n_spectra <- length(spectra(mse))
  if (n_spectra == 0) {
    stop("MsExperiment contains 0 spectra.")
  }

  # Check peak detection results
  if (check_peaks && inherits(mse, "XcmsExperiment")) {
    cp <- chromPeaks(mse)
    if (is.null(cp) || nrow(cp) == 0) {
      stop("No chromatographic peaks found. Run findChromPeaks() first.")
    }

    peaks_per_sample <- table(cp[, "sample"])
    low_samples <- names(peaks_per_sample)[peaks_per_sample < min_peaks_per_sample]
    if (length(low_samples) > 0) {
      warning(length(low_samples), " sample(s) have fewer than ",
              min_peaks_per_sample, " peaks. Consider checking data quality.")
    }
  }

  message("\u2713 MsExperiment validated: ", n, " samples, ",
          n_spectra, " spectra")
  invisible(TRUE)
}


#' Validate a SummarizedExperiment result object
#'
#' Checks feature counts, missing value rates, and assay integrity.
#'
#' @param res A SummarizedExperiment object
#' @param assay_name Assay to check (default: first available)
#' @param max_na_pct Maximum acceptable NA percentage (warns above)
#' @param min_features Minimum expected features (warns below)
#' @return TRUE invisibly
validate_result <- function(res,
                            assay_name = NULL,
                            max_na_pct = 50,
                            min_features = 10) {

  if (!inherits(res, "SummarizedExperiment")) {
    stop("Object is not a SummarizedExperiment (got ", class(res)[1], ")")
  }

  if (nrow(res) == 0) stop("SummarizedExperiment has 0 features.")
  if (ncol(res) == 0) stop("SummarizedExperiment has 0 samples.")

  if (nrow(res) < min_features) {
    warning("Only ", nrow(res), " features detected (expected >= ", min_features,
            "). Check peak detection parameters.")
  }

  # Check assay
  if (is.null(assay_name)) assay_name <- assayNames(res)[1]
  if (!assay_name %in% assayNames(res)) {
    stop("Assay '", assay_name, "' not found. Available: ",
         paste(assayNames(res), collapse = ", "))
  }

  mat <- assay(res, assay_name)
  na_pct <- sum(is.na(mat)) / length(mat) * 100

  if (na_pct > max_na_pct) {
    warning("High missing value rate in assay '", assay_name, "': ",
            round(na_pct, 1), "% (threshold: ", max_na_pct, "%).",
            "\n  Consider reviewing peak detection or gap filling.")
  }

  # Check for all-NA features
  all_na_features <- rowSums(!is.na(mat)) == 0
  if (any(all_na_features)) {
    warning(sum(all_na_features),
            " feature(s) have NA across all samples in assay '",
            assay_name, "'.")
  }

  # Check for negative values
  neg_vals <- sum(mat < 0, na.rm = TRUE)
  if (neg_vals > 0) {
    warning(neg_vals, " negative values found in assay '", assay_name, "'.",
            " This is unexpected for abundance data.")
  }

  message("\u2713 Result validated: ", nrow(res), " features x ", ncol(res),
          " samples, ", round(na_pct, 1), "% NA in '", assay_name, "'")
  invisible(TRUE)
}


#' Validate matching results
#'
#' Checks that the matched data frame has expected structure and
#' reasonable quality metrics.
#'
#' @param mtched_data Matched data frame from match_features_to_database()
#' @param min_matches Minimum expected matches (warns below)
#' @param max_ppm_error Maximum acceptable median ppm error (warns above)
#' @return TRUE invisibly
validate_matches <- function(mtched_data,
                             min_matches = 5,
                             max_ppm_error = 15) {

  if (!is.data.frame(mtched_data)) {
    stop("mtched_data must be a data.frame")
  }

  if (nrow(mtched_data) == 0) {
    stop("No matches found. Check:\n",
         "  - PPM and RT tolerance parameters\n",
         "  - RT adjustment model quality\n",
         "  - Database polarity (pos/neg)")
  }

  if (nrow(mtched_data) < min_matches) {
    warning("Very few matches (", nrow(mtched_data), "). Expected at least ",
            min_matches, ".")
  }

  # Check required columns
  expected_cols <- c("feature_id", "mzmed", "rtmed",
                     "target_lipid_name_unique")
  missing_cols <- setdiff(expected_cols, colnames(mtched_data))
  if (length(missing_cols) > 0) {
    stop("Matched data missing expected columns: ",
         paste0("'", missing_cols, "'", collapse = ", "))
  }

  # PPM error distribution
  if ("ppm_error" %in% colnames(mtched_data)) {
    med_ppm <- median(abs(mtched_data$ppm_error), na.rm = TRUE)
    if (med_ppm > max_ppm_error) {
      warning("Median absolute ppm error is ", round(med_ppm, 1),
              " (threshold: ", max_ppm_error,
              "). m/z calibration may need review.")
    }
  }

  # Duplicate feature_id check (informational)
  n_dup <- sum(duplicated(mtched_data$feature_id))
  n_unique <- length(unique(mtched_data$feature_id))
  n_lipids <- length(unique(mtched_data$target_lipid_name_unique))

  message("\u2713 Matches validated: ", nrow(mtched_data), " rows, ",
          n_unique, " unique features, ", n_lipids, " unique lipids",
          if (n_dup > 0) paste0(" (", n_dup, " feature ambiguities)"))
  invisible(TRUE)
}


#' Validate RT adjustment fit quality
#'
#' Checks R-squared and residuals of the RT adjustment model.
#'
#' @param rt_fit List returned by fit_rt_adjustment()
#' @param min_r_squared Minimum acceptable R-squared (warns below)
#' @param max_residual Maximum acceptable residual in seconds (warns above)
#' @return TRUE invisibly
validate_rt_adjustment <- function(rt_fit,
                                   min_r_squared = 0.90,
                                   max_residual = 15) {

  if (!is.list(rt_fit) || !"fit" %in% names(rt_fit)) {
    stop("rt_fit must be a list with a 'fit' element (from fit_rt_adjustment)")
  }

  fit <- rt_fit$fit
  r2 <- if (inherits(fit, "scam")) summary(fit)$r.sq else summary(fit)$r.squared

  if (r2 < min_r_squared) {
    warning("RT adjustment R\u00b2 = ", round(r2, 4),
            " (threshold: ", min_r_squared,
            "). Model fit is poor — review reference lipid EICs.")
  }

  # Check residuals
  resids <- abs(residuals(fit))
  max_res <- max(resids, na.rm = TRUE)
  if (max_res > max_residual) {
    warning("Largest RT adjustment residual is ", round(max_res, 1),
            "s (threshold: ", max_residual,
            "s). Some reference lipids may be poorly detected.")
  }

  # Check for NAs in experimental RT
  n_na <- sum(is.na(rt_fit$exp_rt))
  if (n_na > 0) {
    warning(n_na, " reference lipid(s) were not detected. ",
            "Consider removing them from the reference list.")
  }

  message("\u2713 RT adjustment validated: R\u00b2 = ", round(r2, 4),
          ", max residual = ", round(max_res, 1), "s")
  invisible(TRUE)
}

# =============================================================================
# PACKAGE LOADING
# =============================================================================

#' Load all required packages for the lipidomics workflow
#'
#' @param verbose Logical, whether to print loading messages
#' @return NULL (invisibly)
#' @examples
#' load_lipid_packages()
load_lipid_packages <- function(verbose = TRUE) {
  required_packages <- c(
    # Core data handling
    "knitr", "readxl", "writexl",
    # MS data handling
    "MsExperiment", "MsIO", "alabaster.se", "MsBackendMetaboLights",
    "SummarizedExperiment", "xcms", "Spectra", "MetaboCoreUtils",
    # SQL backend (for SQLite data loading)
    "MsBackendSql", "RSQLite",
    # Statistics
    "limma", "matrixStats", "scam",
    # Visualization
    "pander", "RColorBrewer", "pheatmap", "vioplot",
    "ggplot2", "ggfortify", "gridExtra",
    # Database & annotation
    "BiocFileCache", "dbplyr", "AnnotationHub", "CompoundDb",
    "MetaboAnnotation", "enviPat",
    # Misc
    "ggVennDiagram", "UpSetR"
  )

  for (pkg in required_packages) {
    if (verbose) message("Loading: ", pkg)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  # Verify MsIO version
  msio_ver <- as.character(packageVersion("MsIO"))
  if (msio_ver != "0.0.15") {
    warning("MsIO version ", msio_ver, " is loaded, but 0.0.15 is required. ",
            "Run: install.packages('MsIO', repos = c('https://rformassspectrometry.r-universe.dev', 'https://cloud.r-project.org'))
")
  }

  if (verbose) message("\n✓ All packages loaded successfully!")
  invisible(NULL)
}

# =============================================================================
# FOLDER SETUP
# =============================================================================

#' Set up the folder structure for the analysis
#'
#' Creates all necessary output directories for figures and objects.
#'
#' @param polarity Character, either "pos" or "neg"
#' @param base_path Character, base path for the project (default: current dir)
#' @return A list with all folder paths
#' @examples
#' folders <- setup_folders("pos")
setup_folders <- function(polarity = c("pos", "neg"), base_path = ".") {
  polarity <- match.arg(polarity)

  folders <- list(
    base = base_path,
    data = file.path(base_path, "data"),
    objects = file.path(base_path, "objects"),
    figures = file.path(base_path, "figures"),
    eic_is = file.path(base_path, "figures", "EIC_internal_standards"),
    iso_pattern = file.path(base_path, "figures", "iso_pattern_check"),
    peak_detection = file.path(base_path, "figures", "peak_detection_ref_lipid"),
    ref_lipid = file.path(base_path, "figures", "ref_lipid_image")
  )

  # Create directories
  for (folder in folders) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    }
  }

  message("✓ Folder structure created for ", toupper(polarity), " ionization mode")
  return(folders)
}

# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

#' Load data from a SQLite database
#'
#' Loads MS data from a pre-built SQLite database. This is much faster than
#' reading mzML files directly. Use the R/create_sqlite_database.R script
#' to create a database from mzML files first.
#'
#' @param db_path Path to the SQLite database file
#' @param sample_data Data frame with sample metadata. Must contain a
#'   'file_name' column matching the original mzML file names.
#' @return An MsExperiment object with linked sample data
#'
#' @examples
#' seq_data <- readxl::read_xlsx("seq_pos.xlsx") |> as.data.frame()
#' mse <- load_from_sqlite("data/MICROSAMPLING_pos.sqlite", seq_data)
load_from_sqlite <- function(db_path, sample_data) {
  # Check database exists
  if (!file.exists(db_path)) {
    stop("SQLite database not found: ", db_path,
         "\n\nTo create a database from mzML files, run:",
         "\n  create_sqlite_database(STUDY_ID, POLARITY, SEQ_FILE)")
  }

  # Check sample_data has file_name

  if (!"file_name" %in% colnames(sample_data)) {
    stop("sample_data must contain a 'file_name' column")
  }

  # Load required packages
  if (!requireNamespace("MsBackendSql", quietly = TRUE)) {
    stop("Package 'MsBackendSql' required. Install with:\n",
         "  BiocManager::install('MsBackendSql')")
  }

  message("Loading data from SQLite: ", basename(db_path))

  # Load spectra from database
  dbf <- normalizePath(db_path)
  s <- Spectra(dbname = dbf, source = MsBackendSql::MsBackendOfflineSql(),
               drv = RSQLite::SQLite())

  # Add base file name for linking
  s$base_file <- basename(s$dataOrigin)
  s$base_file <- sub(".*[\\\\/]", "", s$base_file)

  # Create MsExperiment and link sample data
  mse <- MsExperiment(spectra = s, sampleData = sample_data)
  mse <- linkSampleData(mse, with = "sampleData.file_name = spectra.base_file")

  message("✓ Loaded ", length(mse), " samples from SQLite database")
  return(mse)
}


#' Load data from MetaboLights repository
#'
#' Downloads and loads MS data directly from a MetaboLights study.
#' Sample metadata columns are mapped according to `column_mapping`.
#'
#' @param mtbls_id MetaboLights study ID (e.g., "MTBLS10722")
#' @param assay_name Name of the assay file in MetaboLights
#' @param column_mapping Named list mapping standard names (list names) to
#'   MetaboLights column names (list values). Required: sample_name, file_name,
#'   sample_type. You can add any additional columns you need.
#' @param filter_column Optional column name for filtering samples
#' @param filter_value Optional value to filter by (keeps matching rows)
#' @param exclude_sample_types Character vector of sample_type values to exclude
#'   (e.g., c("Blank", "blank"))
#' @param polarity Polarity to add to sample data ("pos" or "neg")
#'
#' @details
#' ## Finding column names for your MetaboLights study:
#'
#' 1. Go to https://www.ebi.ac.uk/metabolights/MTBLS{your_id}
#' 2. Download the assay file (a_*.txt) or sample file (s_*.txt)
#' 3. Open in Excel/R to see available column names
#' 4. Create your column_mapping based on which columns you need
#'
#' ## Example column_mapping:
#' column_mapping <- list(
#'   sample_name = "Sample Name",
#'   file_name = "Raw Spectral Data File",
#'   sample_type = "Factor Value[Sample type]",
#'   injection_index = "Comment[injection_index]",
#'   batch = "Factor Value[Batch]"
#' )
#'
#' The list names become the new column names, the values are the original
#' MetaboLights column names. Add or remove entries as needed for your study.
#'
#' @return An MsExperiment object with mapped sample data columns
#'
#' @examples
#' # Load MTBLS10722 positive mode data
#' mse <- load_from_metaboLights(
#'   mtbls_id = "MTBLS10722",
#'   assay_name = "a_MTBLS10722_LC-MS_positive_reverse-phase_metabolite_profiling.txt",
#'   column_mapping = list(
#'     sample_name = "Sample Name",
#'     file_name = "Raw Spectral Data File",
#'     sample_type = "Factor Value[Sample type]",
#'     injection_index = "Comment[injection_index]"
#'   ),
#'   filter_column = "Factor Value[Cohort]",
#'   filter_value = "Cohort_study",
#'   exclude_sample_types = "Blank",
#'   polarity = "pos"
#' )
load_from_metaboLights <- function(mtbls_id,
                                    assay_name,
                                    column_mapping,
                                    filter_column = NULL,
                                    filter_value = NULL,
                                    exclude_sample_types = NULL,
                                    polarity = NULL) {
  # Validate column_mapping
  required <- c("sample_name", "file_name", "sample_type")
  missing <- setdiff(required, names(column_mapping))
  if (length(missing) > 0) {
    stop("column_mapping must include: ", paste(missing, collapse = ", "),
         "\n\nExample:\n",
         "  column_mapping = list(\n",
         "    sample_name = 'Sample Name',\n",
         "    file_name = 'Raw Spectral Data File',\n",
         "    sample_type = 'Factor Value[Sample type]'\n",
         "  )")
  }

  message("Loading data from MetaboLights: ", mtbls_id)
  message("  Assay: ", assay_name)

  # Create MetaboLights parameter
  param <- MsBackendMetaboLights::MetaboLightsParam(
    mtblsId = mtbls_id,
    assayName = assay_name,
    filePattern = ".mzML"
  )

  # Load data
  mse <- readMsObject(
    MsExperiment(),
    param,
    keepOntology = FALSE,
    keepProtocol = FALSE,
    simplify = TRUE
  )

  # Apply filter if specified
  if (!is.null(filter_column) && !is.null(filter_value)) {
    if (!filter_column %in% colnames(sampleData(mse))) {
      stop("Filter column '", filter_column, "' not found in sample data.\n",
           "Available columns:\n  ",
           paste(colnames(sampleData(mse)), collapse = "\n  "))
    }
    mse <- mse[sampleData(mse)[[filter_column]] == filter_value]
    message("  Filtered to ", length(mse), " samples where ",
            filter_column, " == '", filter_value, "'")
  }

  # Map column names
  sad <- sampleData(mse)

  # Check all mapping columns exist
  all_cols <- unlist(column_mapping)
  missing_cols <- setdiff(all_cols, colnames(sad))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in MetaboLights data: ",
         paste(missing_cols, collapse = ", "),
         "\n\nAvailable columns:\n  ",
         paste(colnames(sad), collapse = "\n  "))
  }

  # Create new data frame with mapped columns
  new_sad <- data.frame(stringsAsFactors = FALSE)
  for (new_name in names(column_mapping)) {
    old_name <- column_mapping[[new_name]]
    new_sad[[new_name]] <- sad[[old_name]]
  }

  # Convert injection_index to numeric if present
  if ("injection_index" %in% names(new_sad)) {
    new_sad$injection_index <- as.numeric(new_sad$injection_index)
  }

  # Add polarity if provided
  if (!is.null(polarity)) {
    new_sad$polarity <- polarity
  }

  sampleData(mse) <- DataFrame(new_sad)

  # Exclude sample types if specified
  if (!is.null(exclude_sample_types)) {
    keep <- !sampleData(mse)$sample_type %in% exclude_sample_types
    mse <- mse[keep]
    message("  Excluded ", sum(!keep), " samples of type: ",
            paste(exclude_sample_types, collapse = ", "))
  }

  # Sort by injection index if available
  if ("injection_index" %in% colnames(sampleData(mse))) {
    mse <- mse[order(sampleData(mse)$injection_index), ]
  }

  message("✓ Loaded ", length(mse), " samples from MetaboLights")
  return(mse)
}

# =============================================================================
# COLOR PALETTE MANAGEMENT
# =============================================================================

#' Create a color palette for a categorical variable
#'
#' Generates a named color vector for the unique values in a column.
#' Can use predefined palettes or auto-generate colors.
#'
#' @param values Character vector of values to create colors for, or
#'   a data frame/MsExperiment from which to extract values
#' @param column Column name to extract values from (if values is df/mse)
#' @param palette RColorBrewer palette name or "auto" for automatic selection
#' @param custom_colors Optional named vector of custom colors to use
#' @return Named character vector of colors
#'
#' @examples
#' # From a vector
#' pal <- create_color_palette(c("QC", "Sample", "Blank"))
#'
#' # From sample data
#' pal <- create_color_palette(mse, column = "sample_type")
#'
#' # With custom colors
#' pal <- create_color_palette(mse, "sample_type",
#'   custom_colors = c("QC" = "grey", "Sample" = "blue"))
create_color_palette <- function(values,
                                  column = NULL,
                                  palette = "Set1",
                                  custom_colors = NULL) {

  # Extract values from MsExperiment or data frame if needed
  if (inherits(values, "MsExperiment")) {
    if (is.null(column)) stop("column must be specified for MsExperiment")
    values <- sampleData(values)[[column]]
  } else if (is.data.frame(values)) {
    if (is.null(column)) stop("column must be specified for data frame")
    values <- values[[column]]
  }

  unique_vals <- unique(as.character(values))
  n_vals <- length(unique_vals)

  # Use custom colors if provided
  if (!is.null(custom_colors)) {
    # Start with custom, fill missing with auto
    colors <- custom_colors[unique_vals]
    missing <- is.na(colors)
    if (any(missing)) {
      n_missing <- sum(missing)
      if (palette == "auto" || n_missing > 9) {
        auto_cols <- grDevices::hcl.colors(n_missing, "Dark 3")
      } else {
        auto_cols <- RColorBrewer::brewer.pal(
          max(3, n_missing), palette
        )[seq_len(n_missing)]
      }
      colors[missing] <- auto_cols
    }
    names(colors) <- unique_vals
    return(colors)
  }

  # Auto-generate colors
  if (palette == "auto" || n_vals > 9) {
    colors <- grDevices::hcl.colors(n_vals, "Dark 3")
  } else {
    colors <- RColorBrewer::brewer.pal(max(3, n_vals), palette)[seq_len(n_vals)]
  }

  names(colors) <- unique_vals
  return(colors)
}

#' Get sample colors from a palette
#'
#' Maps sample values to colors using a named palette.
#'
#' @param values Vector of sample values (e.g., sample types)
#' @param palette Named vector of colors (from create_color_palette)
#' @return Vector of colors matching the input values
#'
#' @examples
#' palette <- create_color_palette(mse, "sample_type")
#' colors <- get_sample_colors(sampleData(mse)$sample_type, palette)
get_sample_colors <- function(values, palette) {
  colors <- palette[as.character(values)]
  if (any(is.na(colors))) {
    warning("Some values not found in palette, using grey")
    colors[is.na(colors)] <- "grey50"
  }
  return(unname(colors))
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

#' Plot Base Peak Chromatogram (BPC)
#'
#' @param mse MsExperiment object
#' @param palette Named color palette (or will auto-generate from color_by)
#' @param color_by Column name in sampleData to color by (default: "sample_type")
#' @param main Plot title
#' @param alpha Transparency (0-1, default: 0.5)
#' @param save_path Optional path to save the plot
#' @return NULL (invisibly)
#'
#' @examples
#' # Using default sample_type coloring
#' plot_bpc(mse)
#'
#' # Color by a different column
#' plot_bpc(mse, color_by = "batch")
#'
#' # With custom palette
#' my_pal <- c("QC" = "black", "Sample" = "blue")
#' plot_bpc(mse, palette = my_pal)
plot_bpc <- function(mse,
                     palette = NULL,
                     color_by = "sample_type",
                     main = "BPC",
                     alpha = 0.5,
                     save_path = NULL) {

  # Get color values from sample data
  color_values <- sampleData(mse)[[color_by]]
  if (is.null(color_values)) {
    stop("Column '", color_by, "' not found in sampleData")
  }

  # Create or use palette
  if (is.null(palette)) {
    palette <- create_color_palette(color_values)
  }

  col_sample <- get_sample_colors(color_values, palette)

  # Add alpha transparency
  col_alpha <- paste0(col_sample, format(as.hexmode(round(alpha * 255)), width = 2))

  bpc <- chromatogram(mse, aggregationFun = "max")

  if (!is.null(save_path)) {
    png(save_path, width = 2000, height = 1500, res = 300)
  }

  plot(bpc, col = col_alpha, main = main, lwd = 1.5)
  grid()
  legend("topright", col = palette,
         legend = names(palette), lty = 1, lwd = 2, cex = 0.6,
         horiz = TRUE, bty = "n")

  if (!is.null(save_path)) dev.off()
  invisible(NULL)
}

#' Plot EIC for reference lipids and save to folder
#'
#' @param eic_object XChromatogram object with EIC data
#' @param fdata Feature data with short_name, rt columns
#' @param palette Named color palette
#' @param color_by Column in eic_object phenoData to color by
#' @param alpha Transparency (0-1)
#' @param output_dir Directory to save plots
#' @param suffix Optional suffix for filenames
#' @return NULL (invisibly)
plot_eic_batch <- function(eic_object,
                           fdata,
                           palette = NULL,
                           color_by = "sample_type",
                           alpha = 0.5,
                           output_dir,
                           suffix = "") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Get colors
  color_values <- pData(eic_object)[[color_by]]
  if (is.null(color_values)) {
    warning("Column '", color_by, "' not found, using default colors")
    col_sample <- rep("black", ncol(eic_object))
  } else {
    if (is.null(palette)) {
      palette <- create_color_palette(color_values)
    }
    col_sample <- get_sample_colors(color_values, palette)
  }

  col_alpha <- paste0(col_sample, format(as.hexmode(round(alpha * 255)), width = 2))

  for (i in seq_len(nrow(eic_object))) {
    eic <- eic_object[i, ]
    filename <- paste0(rownames(fdata)[i], suffix, ".png")
    filepath <- file.path(output_dir, filename)

    png(filepath, width = 800, height = 600, res = 120)
    plot(eic,
         main = rownames(fdata)[i],
         cex.axis = 0.8, cex.main = 0.8,
         col = col_alpha)
    grid()
    abline(v = fdata$rt[i], col = "red", lty = 3)
    if (!is.null(palette)) {
      legend("topright", col = palette, legend = names(palette),
             lty = 1, lwd = 2, cex = 0.6, bty = "n")
    }
    dev.off()
  }

  message("✓ Saved ", nrow(eic_object), " EIC plots to: ", output_dir)
  invisible(NULL)
}

#' Plot isotope pattern comparison (mirror plot)
#'
#' @param mtched_data Matched data frame
#' @param iso_spectra Experimental isotope spectra
#' @param theoretical_spectra Theoretical isotope spectra
#' @param output_dir Directory to save plots
#' @return NULL (invisibly)
plot_isotope_patterns <- function(mtched_data, iso_spectra, theoretical_spectra,
                                  output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_len(nrow(mtched_data))) {
    cmp <- paste0(mtched_data$target_lipid_name_unique[i],
                  "_", mtched_data$target_adduct[i])
    cmp <- gsub("[():/ +]", "", cmp)
    fid <- mtched_data$feature_id[i]

    out_file <- file.path(output_dir, paste0(cmp, "_", fid, ".png"))
    png(filename = out_file, width = 8, height = 8, units = "cm",
        res = 600, pointsize = 4)
    plotSpectraMirror(iso_spectra[i], theoretical_spectra[i], ppm = 20,
                      main = paste(cmp, fid, sep = " - "))
    dev.off()
  }

  message("✓ Saved ", nrow(mtched_data), " isotope pattern plots to: ", output_dir)
  invisible(NULL)
}

#' PCA plot of an abundance matrix colored by sample type
#'
#' Log2-transforms, scales, runs prcomp and returns a ggplot.
#'
#' @param mat Numeric matrix (features x samples) of abundances
#' @param sample_type Character vector of sample types, length = ncol(mat)
#' @param palette Named color palette for sample_type values
#' @param title Plot title
#' @return A ggplot object
plot_pca <- function(mat, sample_type, palette, title) {
  vals <- mat |> log2() |> t() |> scale(center = TRUE, scale = TRUE)
  pca_res <- prcomp(vals, scale = FALSE, center = FALSE)
  vals_st <- cbind(vals, sample_type = sample_type)
  ggplot2::autoplot(pca_res, data = vals_st, colour = "sample_type", scale = 0) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title)
}

#' Plot lipid distribution donut chart
#'
#' @param category_counts Table of lipid category counts
#' @param category_names Named vector mapping abbreviations to full names
#' @param title Plot title
#' @return NULL (invisibly)
plot_lipid_donut <- function(category_counts, category_names = NULL,
                             title = "Lipid Distribution") {
  if (!is.null(category_names)) {
    names(category_counts) <- category_names[names(category_counts)]
  }

  bar_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
  bar_colors <- rep(bar_colors, length.out = length(category_counts))

  par(mar = c(1, 1, 2, 1))
  pie(category_counts, col = bar_colors, labels = category_counts,
      border = "white", main = title, cex = 0.9)
  symbols(0, 0, circles = 0.5, inches = FALSE, add = TRUE, bg = "white")
  legend("topright", legend = names(category_counts), fill = bar_colors,
         border = NA, bty = "n", title = "Lipid Category", cex = 0.8)

  invisible(NULL)
}

#' Prepare subclass count data frame for plotting
#'
#' Extracts bracketed abbreviations from the category/subclass rowData columns,
#' validates that each subclass maps to exactly one category, and returns a
#' count data frame ordered for faceted bar plots.
#'
#' @param se SummarizedExperiment with annotated rowData
#' @param cat_col rowData column for lipid category
#' @param sub_col rowData column for lipid subclass
#' @param cat_levels Display order for lipid categories. Categories absent from
#'   this vector are appended at the end. Defaults to the five classes present
#'   in the current study database; adjust if your database covers additional
#'   LIPID MAPS top-level classes (e.g. PK, PR, SL).
#' @return data.frame with columns category (factor), subclass (factor), Freq
prepare_subclass_counts <- function(
        se,
        cat_col = "target_LIPID.CATEGORY..ABBREV.",
        sub_col = "target_LIPID.SUBCLASS..ABBREV.",
        cat_levels = c("FA", "GL", "GP", "SP", "ST")) {

    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    df <- data.frame(
        category = sub(".*\\[(.*)\\]$", "\\1", rd[[cat_col]]),
        subclass  = sub(".*\\[(.*)\\]$", "\\1", rd[[sub_col]]),
        stringsAsFactors = FALSE
    )

    mismatch <- unique(df[, c("subclass", "category")]) |>
        dplyr::group_by(subclass) |>
        dplyr::filter(dplyr::n() > 1) |>
        as.data.frame()
    if (nrow(mismatch) > 0) {
        stop("Subclass abbreviation maps to multiple categories:\n",
             paste(capture.output(print(mismatch)), collapse = "\n"))
    }

    cat_levels <- c(cat_levels, setdiff(unique(df$category), cat_levels))
    counts_df <- as.data.frame(table(category = df$category,
                                     subclass  = df$subclass))
    counts_df <- counts_df[counts_df$Freq > 0, ]
    counts_df$category <- factor(counts_df$category, levels = cat_levels)
    counts_df <- counts_df[order(counts_df$category, -counts_df$Freq), ]
    counts_df$subclass <- factor(counts_df$subclass,
                                 levels = unique(counts_df$subclass))
    counts_df
}

# =============================================================================
# DATA PREPROCESSING HELPERS
# =============================================================================

#' Prepare reference lipid (internal standard) data
#'
#' Loads the internal standard list and prepares EIC extraction ranges.
#'
#' @param file_path Path to the Excel file with internal standards
#' @param rt_window_left Seconds to subtract from RT for rtmin
#' @param rt_window_right Seconds to add to RT for rtmax
#' @param mz_tolerance Absolute m/z tolerance
#' @return Data frame with prepared internal standard data
prepare_reference_lipids <- function(file_path, rt_window_left = 30,
                                     rt_window_right = 15, mz_tolerance = 0.008) {
  intern_standard <- readxl::read_xlsx(file_path) |> as.data.frame()
  intern_standard$rtmin <- intern_standard$RT - rt_window_left
  intern_standard$rtmax <- intern_standard$RT + rt_window_right
  intern_standard$mzmin <- intern_standard$mz - mz_tolerance
  intern_standard$mzmax <- intern_standard$mz + mz_tolerance
  rownames(intern_standard) <- intern_standard$short_name

  return(intern_standard)
}

#' Validate LRS subclass coverage against the lipid database
#'
#' Checks that every lipid subclass present in the lipid database has at
#' least one entry in the Lipid Reference Set (LRS). Emits a warning listing
#' any subclasses that are not represented in the LRS.
#'
#' The lipid database usually stores subclasses as full labels with an
#' abbreviation in square brackets, e.g. "Sphingomyelins [SM]". Reference lipid
#' spreadsheets do not always include a true subclass column; when
#' `lrs_subclass_col` is `NULL`, subclasses are inferred from common LRS columns
#' such as `short_name` and `lipid_name`.
#'
#' @param intern_standard Data frame returned by `prepare_reference_lipids()`.
#' @param lipid_database Data frame returned by `prepare_lipid_database()`.
#' @param lrs_subclass_col Optional column in `intern_standard` holding the
#'   subclass. If `NULL`, infer subclasses from the LRS.
#' @param db_subclass_col Column in `lipid_database` holding the subclass
#'   (default: "LIPID.SUBCLASS..ABBREV.").
#' @return Character vector of uncovered subclasses (invisibly); empty if
#'   every database subclass is represented in the LRS.
validate_lrs_coverage <- function(intern_standard,
                                  lipid_database,
                                  lrs_subclass_col = NULL,
                                  db_subclass_col = "LIPID.SUBCLASS..ABBREV.") {
  extract_subclass_abbrev <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    bracketed <- sub("^.*\\[([^]]+)\\].*$", "\\1", x)
    has_brackets <- grepl("\\[[^]]+\\]", x)
    x[has_brackets] <- bracketed[has_brackets]
    x <- sub("_[0-9]+$", "", x)
    x <- sub("^d([A-Z][A-Za-z0-9 /-]*)$", "\\1", x)
    x <- sub("^([A-Za-z]+)\\s*\\(.*$", "\\1", x)
    x <- sub("^Chol\\s+Ester.*$", "CE", x, ignore.case = TRUE)
    x <- sub("^C[0-9]+\\s+Sphinganine.*$", "SPB", x, ignore.case = TRUE)
    x <- sub("^Sphinganine.*$", "SPB", x, ignore.case = TRUE)
    x <- toupper(trimws(x))
    x[!is.na(x) & x != ""]
  }

  if (!is.null(lrs_subclass_col) &&
      !lrs_subclass_col %in% colnames(intern_standard)) {
    stop("LRS subclass column '", lrs_subclass_col,
         "' not found in intern_standard. Available: ",
         paste(colnames(intern_standard), collapse = ", "))
  }
  if (!db_subclass_col %in% colnames(lipid_database)) {
    db_subclass_aliases <- c(
      "LIPID.SUBCLASS..ABBREV.",
      "LIPID SUBCLASS [ABBREV]",
      "target_LIPID.SUBCLASS..ABBREV.",
      "target_LIPID SUBCLASS [ABBREV]"
    )
    db_subclass_match <- intersect(db_subclass_aliases, colnames(lipid_database))
    if (length(db_subclass_match) > 0) {
      db_subclass_col <- db_subclass_match[[1]]
    } else {
      stop("Database subclass column '", db_subclass_col,
           "' not found in lipid_database. Available: ",
           paste(colnames(lipid_database), collapse = ", "))
    }
  }

  db_subclasses <- unique(lipid_database[[db_subclass_col]])
  db_subclasses <- db_subclasses[!is.na(db_subclasses) & db_subclasses != ""]
  db_subclass_abbrev <- unique(extract_subclass_abbrev(db_subclasses))

  if (!is.null(lrs_subclass_col)) {
    lrs_subclasses <- unique(intern_standard[[lrs_subclass_col]])
  } else {
    candidate_cols <- intersect(
      c("subclass", "lipid_subclass", "LIPID SUBCLASS [ABBREV]", "short_name",
        "lipid_name"),
      colnames(intern_standard)
    )
    if (length(candidate_cols) == 0) {
      stop("Could not infer LRS subclasses. Available intern_standard columns: ",
           paste(colnames(intern_standard), collapse = ", "))
    }
    lrs_subclasses <- unlist(intern_standard[candidate_cols], use.names = FALSE)
  }
  lrs_subclass_abbrev <- unique(extract_subclass_abbrev(lrs_subclasses))

  missing_abbrev <- setdiff(db_subclass_abbrev, lrs_subclass_abbrev)
  missing <- db_subclasses[extract_subclass_abbrev(db_subclasses) %in% missing_abbrev]
  missing <- unique(missing)

  if (length(missing) > 0) {
    warning("LRS is missing at least one lipid for ", length(missing),
            " subclass(es) present in the lipid database: ",
            paste(missing, collapse = ", "),
            "\n  Add an internal standard covering each uncovered subclass ",
            "to the LRS Excel, or harmonise the subclass vocabularies.")
  } else {
    message("✓ LRS covers all ", length(db_subclasses),
            " subclasses in the lipid database")
  }

  invisible(missing)
}

#' Extract EIC for internal standards
#'
#' @param mse MsExperiment object
#' @param intern_standard Prepared internal standard data frame
#' @param sample_subset Optional logical/integer/character vector to subset
#'   samples before chromatogram extraction. Defaults to NULL (all samples).
#'   Useful to speed up iteration when tuning rt_window_left / rt_window_right
#'   on large datasets (e.g. plot only QCs).
#' @return XChromatograms object with EIC data
extract_is_eics <- function(mse, intern_standard, sample_subset = NULL) {
  if (!is.null(sample_subset)) {
    mse <- mse[sample_subset]
  }
  eic_is <- chromatogram(
    mse,
    rt = as.matrix(intern_standard[, c("rtmin", "rtmax")]),
    mz = as.matrix(intern_standard[, c("mzmin", "mzmax")])
  )

  # Add metadata
  fData(eic_is)$mz <- intern_standard$mz
  fData(eic_is)$rt <- intern_standard$RT
  fData(eic_is)$rtmin <- intern_standard$rtmin
  fData(eic_is)$rtmax <- intern_standard$rtmax
  fData(eic_is)$short_name <- intern_standard$short_name
  rownames(fData(eic_is)) <- intern_standard$short_name

  return(eic_is)
}

# =============================================================================
# ANNOTATION HELPERS
# =============================================================================

#' Convert enviPat isopattern output to Spectra object
#'
#' @param x List output from enviPat::isopattern
#' @return Spectra object
isopattern_to_spectra <- function(x) {
  df <- data.frame(msLevel = 1L, formula = names(x))
  df$mz <- lapply(x, function(z) z[, 1L])
  df$intensity <- lapply(x, function(z) z[, 2L])
  Spectra(df)
}

#' Collapse an enviPat isopattern to nominal-resolution isotope peaks
#'
#' Sums the isotopologue fine-structure sub-peaks into one centroid per nominal
#' isotope (intensity-weighted m/z) - what routine Orbitrap/QTOF resolution
#' records. Assumes singly-charged patterns.
#'
#' @param x isopattern matrix (column 1 = m/z, column 2 = abundance).
#' @return Two-column matrix (m/z, abundance), one row per nominal isotope.
collapse_isopattern_nominal <- function(x) {
  mz <- x[, 1L]
  ab <- x[, 2L]
  o <- order(mz); mz <- mz[o]; ab <- ab[o]
  k <- round(mz - mz[1L])                    # isotope index: M+0, M+1, ...
  mzc <- tapply(seq_along(mz), k, function(i) sum(mz[i] * ab[i]) / sum(ab[i]))
  cbind("m/z" = as.numeric(mzc), abundance = as.numeric(tapply(ab, k, sum)))
}

# =============================================================================
# FEATURE MATCHING
# =============================================================================

#' Match features to lipid database (Rank 1 matching)
#'
#' Performs m/z and RT matching between experimental features and the lipid
#' database. This is the primary annotation step using rank 1 (most common)
#' adducts only.
#'
#' The function:
#' 1. Extracts feature data from the SummarizedExperiment
#' 2. Filters database to rank 1 entries only
#' 3. Matches by m/z (with ppm tolerance) and RT (with time tolerance)
#' 4. Returns matched data with database annotations
#'
#' @param res SummarizedExperiment with feature data in rowData
#' @param lipid_database Prepared lipid database (from prepare_lipid_database)
#' @param ppm PPM tolerance for m/z matching (default: 20)
#' @param rt_tol RT tolerance in seconds (default: 20)
#' @param verbose Print progress messages
#'
#' @return List with:
#'   - mtched_data: Data frame of matches with all annotations
#'   - query: Feature data frame (needed for adduct matching)
#'   - lipid_db_r1: Rank 1 database subset (for coverage calculation)
#'
#' @examples
#' match_result <- match_features_to_database(
#'   res, lipid_database,
#'   ppm = 20, rt_tol = 20
#' )
#' mtched_data <- match_result$mtched_data
match_features_to_database <- function(res,
                                        lipid_database,
                                        ppm = 20,
                                        rt_tol = 20,
                                        verbose = TRUE) {

  if (verbose) message("Matching features to lipid database...")

  # Prepare query (experimental features)
  query <- as.data.frame(SummarizedExperiment::rowData(res))
  query$feature_id <- rownames(query)

  # Prepare target (rank 1 only - primary adducts)
  lipid_db_r1 <- lipid_database[lipid_database$rank == 1, ]
  target <- as.data.frame(lipid_db_r1)

  if (verbose) message("  - Query: ", nrow(query), " features")
  if (verbose) message("  - Target: ", nrow(target), " database entries (",
                       length(unique(lipid_db_r1$lipid_name_unique)), " unique lipids)")

  # Configure matching parameters
  param_match <- MetaboAnnotation::MzRtParam(
    tolerance = 0.001,
    ppm = ppm,
    toleranceRt = rt_tol
  )

  # Perform matching
  mtch <- MetaboAnnotation::matchValues(
    query, target, param = param_match,
    mzColname = c("mzmed", "mz"),
    rtColname = c("rtmed", "rt_adjusted")
  )

  # Filter to keep only matches
  mtch <- mtch[MetaboAnnotation::whichQuery(mtch)]

  # Extract matched data
  mtched_data <- as.data.frame(MetaboAnnotation::matchedData(mtch))
  mtched_data$feature_id <- sub("\\..*$", "", rownames(mtched_data))

  # Calculate coverage
  n_features <- length(unique(mtched_data$feature_id))
  n_lipids <- length(unique(mtched_data$target_lipid_name_unique))
  coverage <- n_lipids / length(unique(lipid_db_r1$lipid_name_unique)) * 100

  if (verbose) {
    message("  - Matched: ", n_features, " features to ", n_lipids, " lipids")
    message("  - Database coverage: ", round(coverage, 1), "%")
  }

  return(list(
    mtched_data = mtched_data,
    query = query,
    lipid_db_r1 = lipid_db_r1
  ))
}

#' Calculate isotope pattern similarity for matches
#'
#' @param mse MsExperiment object
#' @param mtched_data Matched data frame
#' @param polarity "pos" or "neg"
#' @param isopeak_threshold Minimum number of isotope peaks (default: 2)
#' @param similarity_threshold Minimum similarity score (default: 0.78)
#' @param batch_size Maximum number of features per sequential batch (default:
#'   100). All samples and candidate annotations for a feature stay together.
#'   Use Inf to calculate all features in one call.
#' @return Updated mtched_data with isopeak_count and isopeak_sim columns
calculate_isotope_similarity <- function(mse, mtched_data, polarity = "pos",
                                         isopeak_threshold = 2,
                                         similarity_threshold = 0.78,
                                         batch_size = 100L) {
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      is.na(batch_size) || batch_size <= 0 ||
      (is.finite(batch_size) && batch_size != floor(batch_size))) {
    stop("batch_size must be a positive integer or Inf.")
  }
  features <- unique(mtched_data$feature_id)
  if (!length(features)) {
    mtched_data$isopeak_count <- integer()
    mtched_data$isopeak_sim <- numeric()
    return(list(mtched_data = mtched_data, iso_spectra = Spectra(),
                theoretical_spectra = Spectra()))
  }
  if (length(features) <= batch_size) {
    return(.calculate_isotope_similarity_batch(
      mse, mtched_data, polarity, isopeak_threshold, similarity_threshold
    ))
  }

  # Avoid parallel workers retaining copies of each batch's full MS1 spectra.
  previous_param <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(previous_param), add = TRUE)

  batches <- split(features, ceiling(seq_along(features) / batch_size))
  results <- vector("list", length(batches))
  retained_rows <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    rows <- which(mtched_data$feature_id %in% batches[[i]])
    message("Isotope validation batch ", i, "/", length(batches),
            ": ", length(batches[[i]]), " features, ", length(rows), " candidates")
    results[[i]] <- .calculate_isotope_similarity_batch(
      mse, mtched_data[rows, , drop = FALSE], polarity,
      isopeak_threshold, similarity_threshold
    )
    # scalePeaks queues normalization. Spectra requires that processing be
    # applied before objects from separate batches can be concatenated.
    results[[i]]$iso_spectra <- applyProcessing(results[[i]]$iso_spectra)
    results[[i]]$theoretical_spectra <- applyProcessing(results[[i]]$theoretical_spectra)
    retained_rows[[i]] <- rows[match(
      rownames(results[[i]]$mtched_data), rownames(mtched_data)[rows]
    )]
    # Only the small filtered isotope spectra survive into the next batch.
    gc(verbose = FALSE)
  }

  # Feature batches can interleave candidate rows; restore the original order
  # in both the annotation table and its corresponding mirror-plot spectra.
  rows <- unlist(retained_rows, use.names = FALSE)
  original_order <- order(rows)
  matched <- do.call(rbind, lapply(results, `[[`, "mtched_data"))
  rownames(matched) <- rownames(mtched_data)[rows]
  list(
    mtched_data = matched[original_order, , drop = FALSE],
    iso_spectra = concatenateSpectra(lapply(results, `[[`, "iso_spectra"))[original_order],
    theoretical_spectra = concatenateSpectra(lapply(results, `[[`, "theoretical_spectra"))[original_order]
  )
}

# Scientific calculation for one complete feature batch. Keep all spectra for
# each feature together; batching samples would change the combined spectrum.
.calculate_isotope_similarity_batch <- function(mse, mtched_data, polarity,
                                                isopeak_threshold,
                                                similarity_threshold) {
  data(isotopes, package = "enviPat", envir = environment())

  # Set charge based on polarity
  charge <- ifelse(polarity == "neg", -1, 1)

  sp <- featureSpectra(mse, msLevel = 1L, skipFilled = TRUE,
                       features = unique(mtched_data$feature_id),
                       method = "closest_rt")
  # Convert subset to in-memory backend (needed for combineSpectra)
  sp <- setBackend(sp, MsBackendMemory())

  # Combine spectra per feature
  csp <- combineSpectra(sp, f = sp$feature_id, p = sp$feature_id,
                        peaks = "intersect", ppm = 2)

  # Filter for isotope patterns
  iso_spectra <- spectrapply(csp, function(z) {
    idx <- isotopologues(peaksData(z)[[1L]], ppm = 10, seedMz = z$feature_mzmed)
    if (length(idx) == 1L) {
      addProcessing(z, function(x, i = idx[[1L]], ...) x[i, , drop = FALSE]) |>
        applyProcessing()
    } else {
      filterMzValues(z, mz = z$feature_mzmed, ppm = 10, tolerance = 0) |>
        applyProcessing()
    }
  }) |>
    concatenateSpectra() |>
    scalePeaks(by = max)

  # Calculate theoretical patterns (suppress electron configuration messages)
  chem_checked <- suppressMessages(suppressWarnings(
    check_chemform(isotopes, mtched_data$target_adduct_formula)
  ))
  ip <- suppressMessages(suppressWarnings(
    isopattern(isotopes, chem_checked$new_formula,
               threshold = 0.001, charge = charge, rel_to = 0)
  ))
  ip <- lapply(ip, collapse_isopattern_nominal)
  theoretical_spectra <- isopattern_to_spectra(ip)
  # Normalize theoretical spectra to [0,1] to match experimental scalePeaks output
  theoretical_spectra <- scalePeaks(theoretical_spectra, by = max)

  # Calculate similarity
  match_indices <- match(mtched_data$feature_id, iso_spectra$feature_id)
  mtched_data$isopeak_count <- lengths(iso_spectra)[match_indices]
  mtched_data$isopeak_sim <- diag(
    compareSpectra(iso_spectra[match_indices], theoretical_spectra, ppm = 20)
  )

  # Filter
  keep <- mtched_data$isopeak_count >= isopeak_threshold &
    mtched_data$isopeak_sim >= similarity_threshold
  mtched_data <- mtched_data[keep, ]

  # Subset spectra to match filtered mtched_data
  keep_indices <- match_indices[keep]
  iso_spectra_filtered <- iso_spectra[keep_indices]
  theoretical_spectra_filtered <- theoretical_spectra[keep]

  return(list(mtched_data = mtched_data,
              iso_spectra = iso_spectra_filtered,
              theoretical_spectra = theoretical_spectra_filtered))
}

#' Resolve sn-1/sn-2 regioisomers of lysoglycerophospholipids
#'
#' Lysoglycerophospholipids annotated as both regioisomers, e.g. LPC(18:1/0:0)
#' (acyl at sn-1) and LPC(0:0/18:1) (acyl at sn-2), are resolved by elution
#' order on RP-LC: the sn-2 isomer elutes before the sn-1 isomer.
#' When the candidate features of a species have at least two distinct RT,
#' the sn-1 annotation is removed from the earliest feature and the sn-2
#' annotation from the latest one. When both isomers map to a single RT they
#' cannot be distinguished and both annotations are kept, so that they are
#' handled by the feature-ambiguity merge downstream.
#'
#' @param mtched_data Matched data frame
#' @return Filtered mtched_data with resolved isomers
resolve_lyso_sn_isomers <- function(mtched_data) {
  # "(0:0/", "(O-0:0/", "(P-0:0/": acyl chain at sn-2
  sn2_pattern <- "\\((O-|P-)?0:0/"
  # "/0:0)" or "/0:0(d7)": acyl chain at sn-1
  sn1_pattern <- "/0:0([)(])"
  is_sn2 <- function(name) grepl(sn2_pattern, name)
  is_sn1 <- function(name) grepl(sn1_pattern, name)
  get_lipid_base <- function(names) {
    gsub(sn1_pattern, "\\1", gsub(sn2_pattern, "(\\1", names))
  }

  nms <- mtched_data$target_lipid_name_unique
  base <- get_lipid_base(nms)
  rt <- round(mtched_data$rtmed, 2)
  keep_row <- rep(TRUE, nrow(mtched_data))
  n_unresolved <- 0L

  for (species in unique(base[is_sn1(nms) | is_sn2(nms)])) {
    idx <- which(base == species)
    if (!(any(is_sn1(nms[idx])) && any(is_sn2(nms[idx])))) next

    unique_rts <- sort(unique(rt[idx]))
    if (length(unique_rts) < 2) {
      n_unresolved <- n_unresolved + 1L
      next
    }
    rt_early <- unique_rts[1]
    rt_late <- unique_rts[length(unique_rts)]
    keep_row[idx[rt[idx] == rt_early & is_sn1(nms[idx])]] <- FALSE
    keep_row[idx[rt[idx] == rt_late & is_sn2(nms[idx])]] <- FALSE
  }

  final_resolved <- mtched_data[keep_row, ]
  final_resolved$ntch_idx <- seq_len(nrow(final_resolved))

  message("✓ sn-1/sn-2 resolution: ", nrow(mtched_data), " -> ",
          nrow(final_resolved), " rows")
  if (n_unresolved > 0L) {
    message("  ", n_unresolved, " species with both isomers at a single ",
            "RT kept as ambiguous")
  }
  return(final_resolved)
}

# =============================================================================
# NORMALIZATION HELPERS
# =============================================================================


#' Apply volume correction factors
#'
#' @param se SummarizedExperiment object
#' @param sample_pattern_factors Named list with patterns and factors
#' @param assay_name Name of assay to correct
#' @param new_assay_name Name for the corrected assay
#' @return SummarizedExperiment with corrected assay
apply_volume_correction <- function(se, sample_pattern_factors,
                                    assay_name = "raw",
                                    new_assay_name = "raw_corr") {
  # Default factors for this study
  if (missing(sample_pattern_factors)) {
    sample_pattern_factors <- list(
      "W" = 2.5,    # W-DBS: smaller volume
      "P" = 1.08    # Plasma: slight adjustment
    )
  }

  factors <- rep(1, ncol(se))
  for (pattern in names(sample_pattern_factors)) {
    idx <- grepl(pattern, colnames(se))
    factors[idx] <- sample_pattern_factors[[pattern]]
  }

  # Convert to matrix to avoid DelayedArray issues with sweep
  mat <- as.matrix(assay(se, assay_name))
  assay(se, new_assay_name) <- sweep(mat, MARGIN = 2, STATS = factors, FUN = "*")
  return(se)
}

#' Normalize by internal standards
#'
#' @param se SummarizedExperiment with annotated features
#' @param input_assay Name of input assay
#' @param output_assay Name for output assay
#' @param is_col Column name in rowData containing IS mapping
#' @param lipid_col Column name in rowData containing lipid names
#' @return SummarizedExperiment with IS-normalized assay
normalize_by_is <- function(se, input_assay = "corr_filled",
                            output_assay = "ISnorm_filled",
                            is_col = "target_IS_norm",
                            lipid_col = "target_lipid.name") {

  assay(se, output_assay) <- assay(se, input_assay)

  is_mapping <- rowData(se)[[is_col]]
  lipid_names <- rowData(se)[[lipid_col]]

  for (i in seq_len(nrow(se))) {
    is_name <- is_mapping[i]

    # Skip if IS not in dataset
    if (!is_name %in% lipid_names) {
      warning(sprintf("IS '%s' for feature '%s' not found. Skipping.",
                      is_name, rownames(se)[i]))
      next
    }

    # Skip if feature IS itself
    if (lipid_names[i] == is_name) next

    if (!is.na(is_name)) {
      is_idx <- which(lipid_names == is_name)
      if (length(is_idx) == 1) {
        assay(se, output_assay)[i, ] <-
          assay(se, input_assay)[i, ] / assay(se, input_assay)[is_idx, ]
      }
    }
  }

  message("✓ Internal standard normalization complete")
  return(se)
}

# =============================================================================
# QUALITY CONTROL
# =============================================================================

#' Filter features by QC RSD threshold
#'
#' @param se SummarizedExperiment object
#' @param threshold RSD threshold (default: 0.3 = 30%)
#' @param qc_col Column name for sample type
#' @param qc_value Value indicating QC samples
#' @return Filtered SummarizedExperiment
filter_by_qc_rsd <- function(se, threshold = 0.3, qc_col = "sample_type",
                             qc_value = "QC") {
  rsd_filter <- RsdFilter(threshold = threshold,
                          qcIndex = colData(se)[[qc_col]] == qc_value)
  se_filtered <- filterFeatures(se, filter = rsd_filter)

  message("✓ RSD filter: ", nrow(se), " -> ", nrow(se_filtered), " features",
          " (", round(nrow(se_filtered) / nrow(se) * 100, 1), "% retained)")
  return(se_filtered)
}

# =============================================================================
# DATABASE ADJUSTMENT
# =============================================================================

#' Fit RT adjustment model using reference lipids
#'
#' @param eic_is EIC object for internal standards
#' @param intern_standard Reference lipid data frame
#' @param param_group PeakDensityParam for grouping
#' @param method Fitting method: "scam" (monotone increasing P-spline, default)
#'   or "poly" (polynomial lm, legacy)
#' @param poly_degree Polynomial degree; only used when method = "poly"
#' @param output_dir Directory to save diagnostic plots
#' @return List with fit model and experimental RT values
fit_rt_adjustment <- function(eic_is, intern_standard, param_group,
                              method = c("scam", "poly"),
                              poly_degree = 6, output_dir = NULL) {
  method <- match.arg(method)

  eic_is_corr <- groupChromPeaks(eic_is, param = param_group)

  exp_rt <- vector("numeric", length = nrow(eic_is_corr))

  for (i in seq_len(nrow(eic_is_corr))) {
    table <- featureDefinitions(eic_is_corr[i, ])
    if (!nrow(table)) {
      message("Compound not detected: ", fData(eic_is_corr)$short_name[i])
      exp_rt[i] <- NA
      next
    }

    # Save plot if output directory provided
    if (!is.null(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      png(filename = file.path(output_dir, paste0("feature_compound_", i, ".png")))
      plot(eic_is_corr[i, ],
           main = rownames(fData(eic_is_corr))[i],
           cex.axis = 0.8, cex.main = 0.8)
      # Show experimental feature RT (blue solid) and reference RT (red dashed)
      if (nrow(table) > 0) {
        if (nrow(table) > 1) {
          feat_rt <- table$rtmed[which.max(
            rowSums(featureValues(eic_is_corr[i, ]), na.rm = TRUE)
          )]
        } else {
          feat_rt <- table$rtmed
        }
        abline(v = feat_rt, col = "blue", lty = 1, lwd = 1.5)
      }
      abline(v = fData(eic_is_corr)$rt[i], col = "red", lty = 3)
      legend("topright",
             legend = c("Experimental RT", "Reference RT"),
             col = c("blue", "red"), lty = c(1, 3), lwd = c(1.5, 1),
             bty = "n", cex = 0.7)
      dev.off()
    }

    if (nrow(table) > 1) {
      idx_max <- which.max(rowSums(featureValues(eic_is_corr[i, ]), na.rm = TRUE))
      rt <- table[idx_max, ]$rtmed
    } else {
      rt <- table$rtmed
    }
    exp_rt[i] <- rt
  }

  ref_rt <- intern_standard$RT

  df_fit <- data.frame(ref_rt = ref_rt, exp_rt = exp_rt)
  df_fit <- df_fit[!is.na(df_fit$ref_rt) & !is.na(df_fit$exp_rt), ]

  if (method == "scam") {
    fit <- scam::scam(exp_rt ~ s(ref_rt, bs = "mpi"), data = df_fit)
    r2 <- summary(fit)$r.sq
  } else {
    fit <- lm(exp_rt ~ poly(ref_rt, poly_degree, raw = TRUE),
              data = df_fit)
    r2 <- summary(fit)$r.squared
  }

  # Calibration range = span of training ref_rt (input axis). Stored as
  # an attribute so apply_rt_adjustment() can use the correct axis for
  # the in/out-of-range decision instead of fit$fitted.values (output).
  attr(fit, "ref_rt_range") <- range(df_fit$ref_rt)

  message("✓ RT adjustment model fitted (R² = ", round(r2, 4), ")")

  return(list(fit = fit, exp_rt = exp_rt, ref_rt = ref_rt))
}

#' Apply RT adjustment to lipid database
#'
#' Applies an RT adjustment model (fitted by `fit_rt_adjustment()`) to a
#' lipid database. Database entries whose reference RT falls within the
#' calibration range are corrected by direct model prediction. Entries
#' outside the calibration range (below the earliest or above the latest
#' reference standard) are handled according to `extrapolate`.
#'
#' @param lipid_database Data frame with lipid database
#' @param fit Fitted model from `fit_rt_adjustment()` (either an `lm`
#'   polynomial fit or a `scam` monotone-spline fit)
#' @param rt_col Column name containing RT in seconds
#' @param extrapolate How to correct entries whose reference RT lies
#'   outside the calibration range spanned by the reference standards.
#'   One of:
#'   \describe{
#'     \item{`"auto"` (default)}{Model-aware behaviour: use natural model
#'       prediction (`predict()`) for `scam` fits, whose monotonicity
#'       constraint keeps extrapolation physically sensible, and the
#'       constant-offset fallback for polynomial (`lm`) fits, which can
#'       diverge wildly outside the fitted range.}
#'     \item{`TRUE`}{Always use natural model prediction (`predict()`)
#'       outside the calibration range, regardless of model type. Use with
#'       care for high-degree polynomial fits.}
#'     \item{`FALSE`}{Always use the constant-offset fallback: apply a
#'       zero-order ("hold last value") shift equal to the adjustment at
#'       the nearest calibration boundary.}
#'   }
#' @return Data frame with `rt_adjusted` column added
apply_rt_adjustment <- function(lipid_database, fit, rt_col = "rt_sd",
                                extrapolate = "auto") {

  if (!(identical(extrapolate, "auto") ||
        identical(extrapolate, TRUE) ||
        identical(extrapolate, FALSE))) {
    stop('`extrapolate` must be "auto", TRUE, or FALSE')
  }
  use_predict <- if (identical(extrapolate, "auto")) {
    inherits(fit, "scam")
  } else {
    isTRUE(extrapolate)
  }

  rt_sd <- lipid_database[[rt_col]]

  # Calibration range = training ref_rt span (input axis), attached by
  # fit_rt_adjustment(). Must be present.
  rng <- attr(fit, "ref_rt_range")
  if (is.null(rng))
    stop("fit is missing 'ref_rt_range' attribute; refit with fit_rt_adjustment()")
  rt_min <- rng[1]
  rt_max <- rng[2]

  in_range_idx <- which(rt_sd >= rt_min & rt_sd <= rt_max)

  corrected_rt <- rep(NA_real_, nrow(lipid_database))
  if (length(in_range_idx))
    corrected_rt[in_range_idx] <- predict(
      fit, newdata = data.frame(ref_rt = rt_sd[in_range_idx])
    )

  below <- which(rt_sd < rt_min)
  if (length(below)) {
    if (use_predict) {
      corrected_rt[below] <- predict(
        fit, newdata = data.frame(ref_rt = rt_sd[below])
      )
    } else {
      offset <- rt_min - predict(
        fit, newdata = data.frame(ref_rt = rt_min)
      )
      corrected_rt[below] <- rt_sd[below] - offset
    }
  }

  above <- which(rt_sd > rt_max)
  if (length(above)) {
    if (use_predict) {
      corrected_rt[above] <- predict(
        fit, newdata = data.frame(ref_rt = rt_sd[above])
      )
    } else {
      offset <- rt_max - predict(
        fit, newdata = data.frame(ref_rt = rt_max)
      )
      corrected_rt[above] <- rt_sd[above] - offset
    }
  }

  lipid_database$rt_adjusted <- corrected_rt

  message("✓ RT adjustment applied. Extrapolated ",
          nrow(lipid_database) - length(in_range_idx),
          " compounds outside fitted range")

  return(lipid_database)
}

# =============================================================================
# LIPID DATABASE PREPARATION
# =============================================================================

#' Prepare lipid database for matching
#'
#' This function performs all necessary preprocessing steps on the lipid database:
#'
#' 1. **Rank filtering**: Removes rank=0 entries (typically unconfirmed identifications)
#' 2. **RT conversion**: Converts retention time from minutes to seconds
#' 3. **Unique naming**: Creates unique lipid identifiers by appending RT rank
#'    (handles lipids with same name but different RT)
#' 4. **Adduct processing**: Standardizes adduct notation and calculates
#'    theoretical m/z using MetaboCoreUtils::adductFormula()
#' 5. **Deuterium handling**: Fixes formula notation for deuterated compounds
#'    (D -> [2H] for calculation, then back)
#' 6. **NA removal**: Removes empty rows from Excel import
#' 7. **RT adjustment**: Applies SCAM RT adjustment based on reference lipids
#'
#' @param db_path Path to the lipid database Excel file
#' @param sheet Sheet number (4 for POS, 5 for NEG)
#' @param polarity "pos" or "neg"
#' @param rt_col Column name for retention time
#' @param rt_fit Fitted RT adjustment model (from fit_rt_adjustment)
#' @param extrapolate Extrapolation policy for database entries outside
#'   the calibration range; forwarded to `apply_rt_adjustment()`. One of
#'   `"auto"` (default; model-aware), `TRUE`, or `FALSE`. See
#'   `?apply_rt_adjustment` for details.
#' @param verbose Print progress messages
#'
#' @return Prepared lipid database data frame with columns:
#'   - rt_sd: RT in seconds
#'   - lipid_name_unique: Unique identifier for each lipid/RT combination
#'   - adduct: Cleaned adduct name
#'   - adduct_formula: Calculated adduct formula
#'   - rt_adjusted: Corrected RT based on experimental reference lipids
#'   - mz: Numeric m/z value
#'
#' @examples
#' lipid_db <- prepare_lipid_database(
#'   db_path = "LipidDatabase_R.xlsx",
#'   sheet = 4L,
#'   polarity = "pos",
#'   rt_col = "RT ESI(+) (min)",
#'   rt_fit = rt_fit
#' )
prepare_lipid_database <- function(db_path,
                                    sheet,
                                    polarity,
                                    rt_col,
                                    rt_fit,
                                    extrapolate = "auto",
                                    verbose = TRUE) {

  if (verbose) message("Loading SRM 1950-derived lipid database from sheet ", sheet, "...")

  # Step 1: Load and filter by rank
  lipid_database <- readxl::read_xlsx(db_path, sheet = sheet)
  n_initial <- nrow(lipid_database)
  lipid_database <- lipid_database[lipid_database$rank != 0, ]
  if (verbose) message("  - Removed ", n_initial - nrow(lipid_database),
                       " rank=0 entries")

  # Step 2: Convert RT to seconds and sort
  lipid_database$rt_sd <- as.numeric(lipid_database[[rt_col]]) * 60
  lipid_database <- lipid_database[order(lipid_database$rt_sd), ]

  # Step 3: Create unique lipid names (same lipid can have multiple RTs)
  get_rt_rank <- function(x) match(x, sort(unique(x)))
  rt_ranks <- ave(
    lipid_database[[rt_col]],
    lipid_database$`lipid name`,
    FUN = get_rt_rank
  )
  lipid_database$lipid_name_unique <- paste(
    lipid_database$`lipid name`, rt_ranks, sep = "_"
  )

  # Step 4: Process adduct information
  lipid_database$adduct <- gsub("Adduct", "", lipid_database$Adduct)
  lookup_ads <- get_adduct_lookup(polarity)

  # Step 5: Fix deuterium formula notation for adductFormula()
  # D notation must be converted to [2H] for MetaboCoreUtils
  lipid_database$fix_form <- gsub("D([0-9]+)", "[2H\\1]",
                                  lipid_database$`MOLECULAR FORMULA`)

  # Step 6: Remove NA rows (empty Excel rows that get imported)
  n_before_na <- nrow(lipid_database)
  lipid_database <- lipid_database[!is.na(lipid_database$`lipid name`), ]
  if (verbose && n_before_na != nrow(lipid_database)) {
    message("  - Removed ", n_before_na - nrow(lipid_database), " empty rows")
  }

  # Step 7: Calculate adduct formulas using MetaboCoreUtils
  if (verbose) message("  - Calculating adduct formulas...")
  custom_ads <- lookup_ads[match(lipid_database$adduct, lookup_ads$name), ]
  rownames(custom_ads) <- NULL
  sladduct <- split(custom_ads, seq(nrow(custom_ads)))

  form <- mapply(
    lipid_database$fix_form,
    adduct = sladduct,
    FUN = MetaboCoreUtils::adductFormula,
    standardize = TRUE
  )
  # Convert back from [2H] to D notation
  form <- gsub("\\[2H([0-9]+)\\]", "D\\1", form)
  lipid_database$adduct_formula <- unname(gsub("^\\[|\\].*$", "", form))

  # Step 8: Apply RT adjustment from reference lipids
  if (verbose) message("  - Applying RT adjustment...")
  lipid_database <- apply_rt_adjustment(lipid_database, rt_fit$fit, "rt_sd",
                                        extrapolate = extrapolate)
  lipid_database$mz <- as.numeric(lipid_database$mz)

  if (verbose) message("✓ Database prepared: ", nrow(lipid_database),
                       " entries, ",
                       length(unique(lipid_database$lipid_name_unique)),
                       " unique lipids")

  return(lipid_database)
}

# =============================================================================
# ADDUCT LOOKUP TABLES
# =============================================================================

#' Get adduct lookup table for positive or negative mode
#'
#' Returns a data frame with adduct parameters needed for m/z calculation.
#' These are used by MetaboCoreUtils::adductFormula().
#'
#' @param polarity "pos" or "neg"
#' @return Data frame with columns: name, mass_multi, mass_add, formula_add,
#'         formula_sub, charge, positive
get_adduct_lookup <- function(polarity = c("pos", "neg")) {
  polarity <- match.arg(polarity)

  if (polarity == "pos") {
    lookup_ads <- data.frame(
      name = c("[M+H]+", "[M+Na]+", "[M+C2H7N2]+", "[M+K]+", "[M+NH4]+"),
      mass_multi = 1.0,
      mass_add = c(1.007276, 22.989220, 59.060380, 38.963160, 18.033830),
      formula_add = c("H", "Na", "C2H7N2", "K", "NH4"),
      formula_sub = "C0",
      charge = 1,
      positive = TRUE,
      stringsAsFactors = FALSE
    )
  } else {
    lookup_ads <- data.frame(
      name = c("[M-H]-", "[M+CH3COO]-", "[M+HCOO]-", "[M+Cl]-",
               "[M-H+(CH3COONa)]-", "[M+HCOO+(CH3COONa)]-",
               "[M+CH3COO+(CH3COONa)]-", "[M+CH3COO+(CH3COONa)2]-",
               "[M+CH3COO+(CH3COONa)3]-"),
      mass_multi = 1.0,
      mass_add = c(-1.007276, 59.01385, 44.99820, 34.96940, 81.00551,
                   127.0110, 141.0267, 223.0395, 305.0524),
      formula_add = c("C0", "CH3COO", "HCOO", "Cl", "C2H3O2Na",
                      "C3H4O4Na", "C4H6O4Na", "C6H9O6Na2", "C8H12O8Na3"),
      formula_sub = c("H", "C0", "C0", "C0", "H", "C0", "C0", "C0", "C0"),
      charge = -1,
      positive = FALSE,
      stringsAsFactors = FALSE
    )
  }

  return(lookup_ads)
}

# =============================================================================
# ADDUCT MATCHING
# =============================================================================

#' Match secondary adducts for confirmed lipid matches
#'
#' After primary (rank 1) matching, this function searches for additional
#' adducts (rank > 1) that support the identification. Finding multiple
#' adducts for the same lipid increases confidence in the annotation.
#'
#' The function:
#' 1. Creates a link table from confirmed rank 1 matches
#' 2. Merges with rank > 1 database entries for the same lipids
#' 3. Searches for these secondary adducts in the feature list
#' 4. Calculates adduct detection ratio (found/expected) per match
#'
#' @param mtched_data Data frame of rank 1 matches (from matchValues)
#' @param lipid_database Full lipid database with all ranks
#' @param query Feature data frame with mzmed and rtmed columns
#' @param ppm PPM tolerance for m/z matching
#' @param rt_tol RT tolerance in seconds (typically tight, e.g., 5s)
#' @param verbose Print progress messages
#'
#' @return Updated mtched_data with adduct_count and adduct_ratio columns
#'
#' @examples
#' mtched_data <- match_adducts(mtched_data, lipid_database, query, ppm = 20)
match_adducts <- function(mtched_data,
                          lipid_database,
                          query,
                          ppm = 20,
                          rt_tol = 5,
                          verbose = TRUE) {

  if (verbose) message("Matching secondary adducts (rank > 1)...")

  # Add index for tracking
  mtched_data$ntch_idx <- seq_len(nrow(mtched_data))

  # Create link table from confirmed matches
  link_table <- mtched_data[, c("feature_id", "target_lipid_name_unique",
                                "rtmed", "ntch_idx")]
  colnames(link_table) <- c("ref_feature_id", "lipid_name_unique",
                            "rt_experimental", "ntch_idx")

  # Get potential adducts (rank > 1) and merge with confirmed lipids
  db_adducts <- lipid_database[lipid_database$rank != 1, ]
  target_r2 <- merge(as.data.frame(db_adducts), link_table,
                     by = "lipid_name_unique")
  target_r2$mz <- as.numeric(target_r2$mz)

  if (nrow(target_r2) == 0) {
    if (verbose) message("  - No secondary adducts to search for")
    mtched_data$adduct_ratio <- 0
    mtched_data$adduct_count <- 0
    return(mtched_data)
  }

  # Match adducts using tight RT tolerance (same compound)
  param_adduct <- MetaboAnnotation::MzRtParam(
    tolerance = 0.001, ppm = ppm, toleranceRt = rt_tol
  )
  mtch_adducts <- MetaboAnnotation::matchValues(
    query, target_r2, param = param_adduct,
    mzColname = c("mzmed", "mz"),
    rtColname = c("rtmed", "rt_experimental")
  )
  mtch_adducts <- mtch_adducts[MetaboAnnotation::whichQuery(mtch_adducts)]

  # Extract adduct hits
  adduct_hits <- as.data.frame(MetaboAnnotation::matchedData(mtch_adducts))
  adduct_hits$feature_id <- sub("\\..*$", "", rownames(adduct_hits))

  # Calculate adduct ratios
  mtched_data$adduct_ratio <- 0
  mtched_data$adduct_count <- 0

  if (nrow(adduct_hits) > 0) {
    counts_per_mtch <- table(adduct_hits$target_ntch_idx)
    total_count <- table(target_r2$ntch_idx)
    names(total_count) <- as.character(names(total_count))

    ratios_per_match <- table(adduct_hits$target_ntch_idx) /
      total_count[names(table(adduct_hits$target_ntch_idx))]

    mtched_data$adduct_count[as.numeric(names(counts_per_mtch))] <-
      as.numeric(counts_per_mtch)
    mtched_data$adduct_ratio[as.numeric(names(ratios_per_match))] <-
      as.numeric(ratios_per_match)
  }

  if (verbose) {
    n_with_adducts <- sum(mtched_data$adduct_ratio > 0)
    message("  - ", n_with_adducts, " matches have supporting adducts (",
            round(n_with_adducts / nrow(mtched_data) * 100, 1), "%)")
  }

  return(mtched_data)
}

# =============================================================================
# AMBIGUITY HANDLING
# =============================================================================

#' Export ambiguity tables for manual curation
#'
#' Creates two Excel files for manual review of ambiguous matches:
#'
#' 1. **Lipid ambiguities**: One lipid matches multiple features
#'    - May indicate isomers or incorrect matches
#'    - User should keep the most likely match based on RT, isotope score, etc.
#'
#' 2. **Feature ambiguities**: One feature matches multiple lipids
#'    - Common for isobaric species
#'    - If not resolved, annotations will be merged (e.g., "LPC 16:0; LPC O-16:1")
#'
#' @param mtched_data Data frame of matches
#' @param output_dir Directory for output files (default: current dir)
#' @param lipid_file Filename for lipid ambiguity table
#' @param feature_file Filename for feature ambiguity table
#' @param verbose Print progress messages
#'
#' @return List with lipid_ambiguities and feature_ambiguities data frames
#'
#' @examples
#' amb <- export_ambiguity_tables(mtched_data)
export_ambiguity_tables <- function(mtched_data,
                                     output_dir = ".",
                                     lipid_file = "lipid_ambiguity_resolution.xlsx",
                                     feature_file = "feature_ambiguity_resolution.xlsx",
                                     verbose = TRUE) {

  # Columns to include in ambiguity tables
  cols_amb <- c("feature_id", "target_mz", "mzmed", "target_rt_adjusted",
                "rtmed", "target_lipid_name_unique", "target_adduct",
                "ppm_error", "score_rt", "target_rank", "score",
                "isopeak_count", "isopeak_sim", "adduct_ratio", "ntch_idx")
  cols_use <- cols_amb[cols_amb %in% colnames(mtched_data)]

  # Lipid ambiguities: one lipid -> multiple features
  key <- paste(mtched_data$target_lipid_name_unique,
               mtched_data$target_adduct, sep = "_")
  is_duplicated <- key %in% key[duplicated(key)]
  amblip <- mtched_data[is_duplicated, cols_use]
  # Auto-resolve (Type 1): keep the feature closest to the expected RT
  # (smallest |score_rt|). Drops redundant features of the same lipid, never a
  # lipid; keep_row is editable.
  amblip$keep_row <- FALSE
  rt_dev <- if ("score_rt" %in% names(amblip)) abs(amblip$score_rt) else
    rep(0, nrow(amblip))
  rt_dev[is.na(rt_dev)] <- Inf
  tie <- if ("isopeak_sim" %in% names(amblip))
    ifelse(is.na(amblip$isopeak_sim), 0, amblip$isopeak_sim) else
    rep(0, nrow(amblip))
  alip_key <- paste(amblip$target_lipid_name_unique, amblip$target_adduct,
                    sep = "_")
  for (kk in unique(alip_key)) {
    grp <- which(alip_key == kk)
    # closest RT wins; break exact RT ties by higher isotope similarity
    amblip$keep_row[grp[order(rt_dev[grp], -tie[grp])][1]] <- TRUE
  }

  # Feature ambiguities: one feature -> multiple lipids (left unresolved;
  # keep_row TRUE, merged downstream to keep all candidate lipids).
  fids_ambiguous <- mtched_data$feature_id[duplicated(mtched_data$feature_id)]
  table_amb <- mtched_data[mtched_data$feature_id %in% fids_ambiguous, cols_use]
  table_amb <- table_amb[order(table_amb$feature_id), ]
  table_amb$keep_row <- TRUE

  # Export files
  if (nrow(amblip) > 0) {
    writexl::write_xlsx(amblip, path = file.path(output_dir, lipid_file))
    if (verbose) {
      message("Exported ", nrow(amblip), " lipid ambiguities to ", lipid_file,
              " (keep_row pre-filled: closest-RT feature per lipid).")
      message("  Review recommended: edit keep_row where you disagree.")
    }
  } else {
    if (verbose) message("No lipid ambiguities to export")
  }

  if (nrow(table_amb) > 0) {
    writexl::write_xlsx(table_amb, path = file.path(output_dir, feature_file))
    if (verbose) message("Exported ", nrow(table_amb),
                         " feature ambiguities to ", feature_file)
  } else {
    if (verbose) message("No feature ambiguities to export")
  }

  return(list(
    lipid_ambiguities = amblip,
    feature_ambiguities = table_amb
  ))
}
