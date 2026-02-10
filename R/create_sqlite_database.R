#' =============================================================================
#' Create SQLite Database from mzML Files
#' =============================================================================
#'
#' This script creates a SQLite database from mzML mass spectrometry files.
#' Using a SQLite database significantly speeds up subsequent data loading.
#'
#' USAGE:
#'   1. Set the parameters in the "Configuration" section below
#'   2. Source this script: source("R/create_sqlite_database.R")
#'   
#' Or run interactively by executing each section.
#' =============================================================================

# ---- Configuration ----
# Modify these parameters for your dataset

# Polarity: "pos" or "neg"
POLARITY <- "pos"

# Directory containing mzML files
DATA_DIR <- paste0(toupper(POLARITY), "_data")

# Sequence file (Excel) with file names
SEQ_FILE <- paste0("seq_", POLARITY, ".xlsx")

# Output database path (default: inside data directory)
OUTPUT_DB <- file.path(DATA_DIR, paste0("pilot_", POLARITY, ".sqlite"))

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
message("Reading sequence file: ", SEQ_FILE)

if (!file.exists(SEQ_FILE)) {
  stop("Sequence file not found: ", SEQ_FILE)
}

seq_data <- read_xlsx(SEQ_FILE, col_names = TRUE) |> as.data.frame()

# Check for file_name column
if (!"file_name" %in% colnames(seq_data)) {
  stop("Sequence file must contain a 'file_name' column")
}

message("  Found ", nrow(seq_data), " samples in sequence file")

# ---- Build File Paths ----
spectra_files <- file.path(DATA_DIR, seq_data$file_name)

# Check all files exist
missing_files <- spectra_files[!file.exists(spectra_files)]
if (length(missing_files) > 0) {
  stop("Missing ", length(missing_files), " mzML files:\n",
       paste("  -", head(missing_files, 5), collapse = "\n"),
       if (length(missing_files) > 5) paste0("\n  ... and ", length(missing_files) - 5, " more"))
}

message("  All ", length(spectra_files), " mzML files found")

# ---- Remove Existing Database ----
if (file.exists(OUTPUT_DB)) {
  message("Removing existing database: ", OUTPUT_DB)
  unlink(OUTPUT_DB)
}

# ---- Create SQLite Database ----
message("\nCreating SQLite database: ", OUTPUT_DB)
message("Processing ", length(spectra_files), " files...")
message("This may take several minutes...\n")

# Track timing
start_time <- Sys.time()

# Create database connection
con <- dbConnect(SQLite(), OUTPUT_DB)

# Create the database
createMsBackendSqlDatabase(dbcon = con, spectra_files)

# Close connection
dbDisconnect(con)

# Report timing
elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
file_size_mb <- round(file.info(OUTPUT_DB)$size / 1024^2, 1)
