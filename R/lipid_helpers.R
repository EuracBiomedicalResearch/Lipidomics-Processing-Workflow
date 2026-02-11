#' ===========================================================================
#' Lipidomics Workflow Helper Functions
#' ===========================================================================
#'
#' This file contains reusable functions for the lipidomics data processing

#' and annotation workflow. Source this file at the beginning of your analysis.
#'
#' @author Philippine Louail, Sara Londono
#' @date 2026
#' ===========================================================================

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
    "limma", "matrixStats",
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
    data = file.path(base_path, paste0(toupper(polarity), "_data")),
    objects = file.path(base_path, "objects"),
    figures = file.path(base_path, "figures"),
    eic_istd = file.path(base_path, "figures", "EIC_internal_standards"),
    iso_pattern = file.path(base_path, paste0("figures/iso_pattern_check",
                                              ifelse(polarity == "neg", "_neg", ""))),
    peak_detection = file.path(base_path, paste0(polarity, "_peak_detection_ref_lipid")),
    ref_lipid = file.path(base_path, paste0("figures/ref_lipid_image",
                                            ifelse(polarity == "neg", "_neg", ""))),
    istd_matched = file.path(base_path, paste0("ISTD_mtched_data",
                                               ifelse(polarity == "neg", "_neg", "")))
  )

  # Create directories
  for (folder in folders) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    }
  }

  message("✓ Folder structure created for ", toupper(polarity), " mode")
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
#' mse <- load_from_sqlite("POS_data/pilot_pos.sqlite", seq_data)
load_from_sqlite <- function(db_path, sample_data) {
  # Check database exists
  if (!file.exists(db_path)) {
    stop("SQLite database not found: ", db_path,
         "\n\nTo create a database from mzML files, run:",
         "\n  source('R/create_sqlite_database.R')")
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

    out_file <- file.path(output_dir,
                          paste0(cmp, "_", iso_spectra$feature_id[i], ".png"))
    png(filename = out_file, width = 8, height = 8, units = "cm",
        res = 600, pointsize = 4)
    plotSpectraMirror(iso_spectra[i], theoretical_spectra[i], ppm = 20,
                      main = paste(cmp, iso_spectra$feature_id[i], sep = " - "))
    dev.off()
  }

  message("✓ Saved ", nrow(mtched_data), " isotope pattern plots to: ", output_dir)
  invisible(NULL)
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

#' Extract EICs for internal standards
#'
#' @param mse MsExperiment object
#' @param intern_standard Prepared internal standard data frame
#' @return XChromatograms object with EIC data
extract_is_eics <- function(mse, intern_standard) {
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
#' @return Updated mtched_data with isopeak_count and isopeak_sim columns
calculate_isotope_similarity <- function(mse, mtched_data, polarity = "pos",
                                         isopeak_threshold = 2,
                                         similarity_threshold = 0.78) {
  data(isotopes, package = "enviPat", envir = environment())

  # Set charge based on polarity
  charge <- ifelse(polarity == "neg", -1, 1)

  # Extract experimental spectra
  spectra(mse) <- setBackend(spectra(mse), MsBackendMemory())
  sp <- featureSpectra(mse, msLevel = 1L, skipFilled = TRUE,
                       features = unique(mtched_data$feature_id),
                       method = "closest_rt")

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

  # Calculate theoretical patterns
  chem_checked <- check_chemform(isotopes, mtched_data$target_adduct_formula)
  ip <- isopattern(isotopes, chem_checked$new_formula,
                   threshold = 0.001, charge = charge, rel_to = 0)
  theoretical_spectra <- isopattern_to_spectra(ip)

  # Calculate similarity
  match_indices <- match(mtched_data$feature_id, iso_spectra$feature_id)
  mtched_data$isopeak_count <- lengths(iso_spectra)[match_indices]
  mtched_data$isopeak_sim <- diag(
    compareSpectra(iso_spectra[match_indices], theoretical_spectra, ppm = 20)
  )

  # Filter
  mtched_data <- subset(mtched_data,
                        isopeak_count >= isopeak_threshold &
                          isopeak_sim >= similarity_threshold)

  return(list(mtched_data = mtched_data,
              iso_spectra = iso_spectra,
              theoretical_spectra = theoretical_spectra))
}

#' Resolve SM1/SM2 isomer ambiguity
#'
#' @param mtched_data Matched data frame
#' @return Filtered mtched_data with resolved isomers
resolve_sm_isomers <- function(mtched_data) {
  # Helper functions
  get_lipid_base <- function(names) gsub("0:0/|/0:0", "", names)
  is_sm2 <- function(name) grepl("0:0/", name, fixed = TRUE)
  is_sm1 <- function(name) grepl("/0:0", name, fixed = TRUE)

  df_clean <- mtched_data
  df_clean$lipid_base_group <- get_lipid_base(df_clean$target_lipid_name_unique)
  df_clean$keep_row <- TRUE

  unique_species <- unique(df_clean$lipid_base_group)

  for (species in unique_species) {
    idx <- which(df_clean$lipid_base_group == species)
    sub_df <- df_clean[idx, ]
    has_sm1 <- any(is_sm1(sub_df$target_lipid_name_unique))
    has_sm2 <- any(is_sm2(sub_df$target_lipid_name_unique))

    if (has_sm1 && has_sm2) {
      unique_rts <- sort(unique(round(sub_df$rtmed, 2)))

      if (length(unique_rts) >= 2) {
        rt_early <- unique_rts[1]
        rt_late <- unique_rts[length(unique_rts)]

        for (i in idx) {
          row_rt <- round(df_clean$rtmed[i], 2)
          row_name <- df_clean$target_lipid_name_unique[i]
          if (row_rt == rt_early && is_sm1(row_name)) {
            df_clean$keep_row[i] <- FALSE
          }
          if (row_rt == rt_late && is_sm2(row_name)) {
            df_clean$keep_row[i] <- FALSE
          }
        }
      } else {
        for (i in idx) {
          if (is_sm2(df_clean$target_lipid_name_unique[i])) {
            df_clean$keep_row[i] <- FALSE
          }
        }
      }
    }
  }

  final_resolved <- df_clean[df_clean$keep_row, ]
  final_resolved$lipid_base_group <- NULL
  final_resolved$keep_row <- NULL
  final_resolved$ntch_idx <- seq_len(nrow(final_resolved))

  message("✓ SM1/SM2 resolution: ", nrow(mtched_data), " -> ", nrow(final_resolved), " rows")
  return(final_resolved)
}

# =============================================================================
# NORMALIZATION HELPERS
# =============================================================================

#' Impute missing values using uniform distribution
#'
#' @param z Numeric vector with potential NA values
#' @return Vector with NAs replaced by random values
na_unidis <- function(z) {
  na <- is.na(z)
  if (any(na)) {
    min_val <- min(z, na.rm = TRUE)
    z[na] <- runif(sum(na), min = min_val / 2, max = min_val)
  }
  z
}

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
# DATABASE CORRECTION
# =============================================================================

#' Fit RT correction model using reference lipids
#'
#' @param eic_is EIC object for internal standards
#' @param intern_standard Reference lipid data frame
#' @param param_group PeakDensityParam for grouping
#' @param poly_degree Polynomial degree for fitting
#' @param output_dir Directory to save diagnostic plots
#' @return List with fit model and experimental RT values
fit_rt_correction <- function(eic_is, intern_standard, param_group,
                              poly_degree = 6, output_dir = NULL) {

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
      abline(v = fData(eic_is_corr)$rt[i], col = "red", lty = 3)
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

  # Fit polynomial model
  fit <- lm(exp_rt ~ poly(ref_rt, poly_degree, raw = TRUE))

  message("✓ RT correction model fitted (R² = ",
          round(summary(fit)$r.squared, 4), ")")

  return(list(fit = fit, exp_rt = exp_rt, ref_rt = ref_rt))
}

#' Apply RT correction to lipid database
#'
#' @param lipid_database Data frame with lipid database
#' @param fit Linear model from fit_rt_correction
#' @param rt_col Column name containing RT in seconds
#' @return Data frame with rt_adjusted column added
apply_rt_correction <- function(lipid_database, fit, rt_col = "rt_sd") {

  rt_sd <- lipid_database[[rt_col]]

  # Find values in fitted range
  in_range_idx <- which(rt_sd >= min(fit$fitted.values) &
                          rt_sd <= max(fit$fitted.values))

  # Apply correction
  corrected_rt <- rep(NA, length = nrow(lipid_database))
  corrected_rt[in_range_idx] <- predict(
    fit, newdata = data.frame(ref_rt = rt_sd[in_range_idx])
  )

  # Extrapolate outside range
  idx <- which(rt_sd < min(fit$fitted.values))
  lidx <- length(idx)
  if (lidx) {
    corrected_rt[idx] <- rt_sd[idx] -
      (rt_sd[lidx + 1L] - corrected_rt[lidx + 1L])
  }

  idx <- which(rt_sd > max(fit$fitted.values))
  if (length(idx)) {
    corrected_rt[idx] <- rt_sd[idx] -
      (rt_sd[idx[1L] - 1L] - corrected_rt[idx[1L] - 1L])
  }

  lipid_database$rt_adjusted <- corrected_rt

  message("✓ RT correction applied. Extrapolated ",
          nrow(lipid_database) - length(in_range_idx), " compounds outside fitted range")

  return(lipid_database)
}

# =============================================================================
# LIPID DATABASE PREPARATION
# =============================================================================

#' Prepare lipid database for matching
#'
#' This function performs all necessary preprocessing steps on the lipid database:
#'
#' 1. **Filter by rank**: Removes rank=0 entries (typically unconfirmed identifications)
#' 2. **RT conversion**: Converts retention time from minutes to seconds
#' 3. **Unique naming**: Creates unique lipid identifiers by appending RT rank
#'    (handles lipids with same name but different RT)
#' 4. **Adduct processing**: Standardizes adduct notation and calculates
#'    theoretical m/z using MetaboCoreUtils::adductFormula()
#' 5. **Deuterium handling**: Fixes formula notation for deuterated compounds
#'    (D -> [2H] for calculation, then back)
#' 6. **NA removal**: Removes empty rows from Excel import
#' 7. **RT correction**: Applies polynomial RT correction based on reference lipids
#'
#' @param db_path Path to the lipid database Excel file
#' @param sheet Sheet number (4 for POS, 5 for NEG)
#' @param polarity "pos" or "neg"
#' @param rt_col Column name for retention time
#' @param rt_fit Fitted RT correction model (from fit_rt_correction)
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
                                    verbose = TRUE) {

  if (verbose) message("Loading lipid database from sheet ", sheet, "...")

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

  # Step 8: Apply RT correction from reference lipids
  if (verbose) message("  - Applying RT correction...")
  lipid_database <- apply_rt_correction(lipid_database, rt_fit$fit, "rt_sd")
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
  amblip$keep_row <- TRUE

  # Feature ambiguities: one feature -> multiple lipids
  fids_ambiguous <- mtched_data$feature_id[duplicated(mtched_data$feature_id)]
  table_amb <- mtched_data[mtched_data$feature_id %in% fids_ambiguous, cols_use]
  table_amb <- table_amb[order(table_amb$feature_id), ]
  table_amb$keep_row <- TRUE

  # Export files
  if (nrow(amblip) > 0) {
    writexl::write_xlsx(amblip, path = file.path(output_dir, lipid_file))
    if (verbose) message("Exported ", nrow(amblip),
                         " lipid ambiguities to ", lipid_file)
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


