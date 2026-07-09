#' =============================================================================
#' Create SQLite Database from mzML Files
#' =============================================================================
#'
#' Converts mzML mass spectrometry files into a SQLite database for faster
#' subsequent data loading. Called from Preprocessing_pos.qmd or
#' Preprocessing_neg.qmd when DATA_SOURCE = "sqlite".
#'
#' @param study_id  Character. Study identifier (e.g. "pilot", "exercise").
#' @param polarity  Character. "pos" or "neg".
#' @param seq_file  Character. Path to the sequence Excel file.
#' @param data_dir  Character. Directory containing the .mzML files.
#'                  Defaults to `data` in the working directory.
#' @param output_db Character or NULL. Path for the output .sqlite file.
#'                  Defaults to `<data_dir>/<study_id>_<polarity>.sqlite`.
#' @return Invisible path to the created database file.
#'
#' @examples
#' create_sqlite_database("pilot", "pos", "seq_pos_pilot.xlsx")
#' =============================================================================
create_sqlite_database <- function(study_id,
                                   polarity = c("pos", "neg"),
                                   seq_file,
                                   data_dir = "data",
                                   output_db = NULL) {
  polarity <- match.arg(polarity)

  if (is.null(output_db)) {
    output_db <- file.path(data_dir, paste0(study_id, "_", polarity, ".sqlite"))
  }

  # ---- Load Required Packages ----
  if (!requireNamespace("MsBackendSql", quietly = TRUE)) {
    stop("Package 'MsBackendSql' required. Install with:\n",
         "  BiocManager::install('MsBackendSql')")
  }
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Package 'RSQLite' required. Install with:\n",
         "  install.packages('RSQLite')")
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' required. Install with:\n",
         "  install.packages('readxl')")
  }

  library(MsBackendSql)
  library(RSQLite)
  library(readxl)

  # ---- Read Sequence File ----
  message("Reading sequence file: ", seq_file)

  if (!file.exists(seq_file)) {
    stop("Sequence file not found: ", seq_file)
  }

  seq_data <- read_xlsx(seq_file, col_names = TRUE) |> as.data.frame()

  if (!"file_name" %in% colnames(seq_data)) {
    stop("Sequence file must contain a 'file_name' column")
  }

  message("  Found ", nrow(seq_data), " samples in sequence file")

  # ---- Build File Paths ----
  spectra_files <- file.path(data_dir, seq_data$file_name)

  missing_files <- spectra_files[!file.exists(spectra_files)]
  if (length(missing_files) > 0) {
    stop("Missing ", length(missing_files), " mzML files:\n",
         paste("  -", head(missing_files, 5), collapse = "\n"),
         if (length(missing_files) > 5)
           paste0("\n  ... and ", length(missing_files) - 5, " more"))
  }

  message("  All ", length(spectra_files), " mzML files found")

  # ---- Remove Existing Database ----
  if (file.exists(output_db)) {
    message("Removing existing database: ", output_db)
    unlink(output_db)
  }

  # ---- Create SQLite Database ----
  message("\nCreating SQLite database: ", output_db)
  message("Processing ", length(spectra_files), " files...")
  message("This may take several minutes...\n")

  start_time <- Sys.time()
  con <- dbConnect(SQLite(), output_db)
  createMsBackendSqlDatabase(dbcon = con, spectra_files)
  dbDisconnect(con)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
  file_size_mb <- round(file.info(output_db)$size / 1024^2, 1)

  message("\n✓ Database created successfully!")
  message("  Path: ", output_db)
  message("  Size: ", file_size_mb, " MB")
  message("  Time: ", elapsed, " minutes")

  invisible(output_db)
}
