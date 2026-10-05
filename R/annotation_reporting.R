# Reporting snapshots only: no matching, filtering or abundance transformations.
annotation_report_frame <- function(x, polarity) {
  if (methods::is(x, "SummarizedExperiment")) {
    d <- as.data.frame(SummarizedExperiment::rowData(x))
    d$feature_id <- rownames(x)
  } else {
    d <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!"feature_id" %in% names(d)) stop("Report input has no feature_id.")
  get_col <- function(candidates, default = NA_character_) {
    found <- intersect(candidates, names(d))
    if (length(found)) d[[found[1]]] else rep(default, nrow(d))
  }
  mode <- get_col("ionization_mode", polarity)
  annotation <- get_col(c("target_lipid_name_unique", "target_lipid.name"))
  lipid <- get_col("target_lipid.name")
  standard <- get_col("target_IS_norm")
  origin <- toupper(as.character(get_col(c("origin", "target_ORIGIN", "ORIGIN"))))
  is_standard <- origin %in% c("IS", "ISTD") |
    (!is.na(lipid) & !is.na(standard) & lipid == standard)
  out <- data.frame(feature_id = as.character(d$feature_id),
                    polarity = as.character(mode), annotation = as.character(annotation),
                    .is_standard = is_standard, stringsAsFactors = FALSE)
  columns <- list(adduct = c("target_Adduct", "target_adduct"),
                  mz = "mzmed", rt_seconds = "rtmed", ppm_error = "ppm_error",
                  rt_error_seconds = "score_rt", isotope_peaks = "isopeak_count",
                  isotope_similarity = "isopeak_sim", adduct_support = "adduct_ratio",
                  qc_cv = "qc_cv")
  for (name in names(columns)) {
    found <- intersect(columns[[name]], names(d))
    if (length(found)) out[[name]] <- d[[found[1]]]
  }
  out
}

annotation_report_pairs <- function(d) {
  empty <- data.frame(feature_key = character(), pair_key = character(),
                      source_row = integer(), lipid = character())
  if (!nrow(d)) return(empty)
  parts <- strsplit(ifelse(is.na(d$annotation), "", d$annotation), ";", fixed = TRUE)
  rows <- rep(seq_len(nrow(d)), lengths(parts))
  lipids <- trimws(unlist(parts, use.names = FALSE))
  keep <- !is.na(lipids) & nzchar(lipids)
  rows <- rows[keep]; lipids <- lipids[keep]
  if (!length(rows)) return(empty)
  keys <- paste(d$polarity[rows], d$feature_id[rows], sep = "::")
  out <- data.frame(feature_key = keys, pair_key = paste(keys, lipids, sep = "||"),
                    source_row = rows, lipid = lipids, stringsAsFactors = FALSE)
  out[!duplicated(out$pair_key), , drop = FALSE]
}

new_annotation_report <- function(study_id, polarity) {
  list(version = 1L, study_id = study_id, polarity = polarity, phases = list())
}

capture_annotation_phase <- function(report, phase, x, annotations = TRUE) {
  if (phase %in% names(report$phases)) stop("Duplicate reporting phase: ", phase)
  d <- annotation_report_frame(x, report$polarity)
  if (!annotations) d$annotation <- NA_character_
  report$phases[[phase]] <- list(data = d, annotations = annotations)
  report
}

save_annotation_report <- function(report, path) {
  if (!"manual_curation" %in% names(report$phases)) stop("Missing manual curation.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(report, path)
}

insert_rank1_mz_phase <- function(report, matches) {
  if ("rank1_mz" %in% names(report$phases)) stop("rank1_mz snapshot already exists.")
  if (!all(c("preprocessed_features", "rank1_mz_rt") %in% names(report$phases)))
    stop("Missing preprocessing or m/z-RT snapshot.")
  snapshot <- annotation_report_frame(matches, report$polarity)
  mz_pairs <- annotation_report_pairs(snapshot)$pair_key
  rt_pairs <- annotation_report_pairs(report$phases$rank1_mz_rt$data)$pair_key
  if (!all(rt_pairs %in% mz_pairs))
    stop("m/z-only candidates do not contain the saved m/z-RT candidates; check inputs and ppm.")
  if (!all(snapshot$feature_id %in% report$phases$preprocessed_features$data$feature_id))
    stop("m/z-only candidates include features absent from preprocessing.")
  report$phases <- c(report$phases["preprocessed_features"],
                     list(rank1_mz = list(data = snapshot, annotations = TRUE)),
                     report$phases[setdiff(names(report$phases), "preprocessed_features")])
  report
}

# Reporting-only refresh for older snapshots. It recalculates only
# mass candidates from saved preprocessing features and the study database;
# all original RT, isotope, curation, abundance and merge snapshots survive.
add_rank1_mz_reporting <- function(study_dir, study_id) {
  if (!requireNamespace("alabaster.se", quietly = TRUE))
    stop("Reading saved preprocessing features requires alabaster.se.")
  helpers <- new.env(parent = globalenv())
  sys.source(file.path(study_dir, "../../R/lipid_helpers.R"), envir = helpers)
  for (suffix in c("pos", "neg")) {
    polarity <- if (suffix == "pos") "positive" else "negative"
    directory <- file.path(study_dir, polarity)
    report_path <- file.path(directory, "objects",
                             paste0(study_id, "_annotation_report_", suffix, ".rds"))
    report <- readRDS(report_path)
    if ("rank1_mz" %in% names(report$phases)) next
    # Read the actual workflow's mass tolerance, rather than assuming a default.
    lines <- readLines(file.path(directory, paste0("Annotation_", suffix, ".qmd")))
    ppm_line <- lines[grepl("^MATCH_PPM\\s*<-", lines)]
    if (length(ppm_line) != 1L) stop("Cannot determine MATCH_PPM for ", directory)
    ppm <- eval(parse(text = ppm_line)[[1L]][[3L]], envir = baseenv())
    if (!is.numeric(ppm) || length(ppm) != 1L || !is.finite(ppm) || ppm < 0)
      stop("Invalid MATCH_PPM for ", directory)
    res <- alabaster.base::readObject(file.path(directory, "objects",
      paste0(study_id, "_preprocessed_res_", suffix)))
    if (!identical(rownames(res), report$phases$preprocessed_features$data$feature_id))
      stop("Saved preprocessing features do not match the annotation snapshot.")
    database <- helpers$prepare_lipid_database(
      file.path(study_dir, "LipidDatabase_R.xlsx"),
      sheet = if (suffix == "pos") 4L else 5L, polarity = suffix,
      rt_col = if (suffix == "pos") "RT ESI(+) (min)" else "RT ESI(\u2013) (min)"
    )
    matches <- helpers$match_features_mz_only(res, database, ppm = ppm)
    report <- insert_rank1_mz_phase(report, matches)
    save_annotation_report(report, report_path)
    message(study_id, " ", polarity, ": added rank1_mz reporting snapshot")
  }
  invisible(NULL)
}

annotation_report_summary <- function(report) {
  previous <- NULL
  rows <- lapply(names(report$phases), function(phase) {
    snapshot <- report$phases[[phase]]
    d <- snapshot$data
    features <- unique(paste(d$polarity, d$feature_id, sep = "::"))
    pairs <- annotation_report_pairs(d)
    out <- data.frame(polarity = report$polarity, phase = phase,
                      features = length(features),
                      annotation_pairs = if (snapshot$annotations) nrow(pairs) else NA_integer_,
                      removed_features = if (is.null(previous)) NA_integer_ else
                        length(setdiff(previous, features)),
                      ambiguous_features = if (snapshot$annotations)
                        sum(table(pairs$feature_key) > 1L) else NA_integer_)
    previous <<- features
    out
  })
  do.call(rbind, rows)
}

annotation_reference_comparison <- function(report) {
  truth <- report$phases$manual_curation$data
  truth <- truth[!truth$.is_standard, , drop = FALSE]
  truth_pairs <- annotation_report_pairs(truth)
  truth_features <- unique(truth_pairs$feature_key)
  phases <- intersect(c("rank1_mz", "rank1_mz_rt", "isotope_filter", "adduct_scored",
                        "ambiguity_auto", "manual_curation"), names(report$phases))
  ratio <- function(n, d) if (d) n / d else NA_real_
  do.call(rbind, lapply(phases, function(phase) {
    d <- report$phases[[phase]]$data
    # Also exclude alternative candidates on curated standard features.
    all_truth <- report$phases$manual_curation$data
    standard_keys <- paste(all_truth$polarity[all_truth$.is_standard],
                           all_truth$feature_id[all_truth$.is_standard], sep = "::")
    keys <- paste(d$polarity, d$feature_id, sep = "::")
    d <- d[!d$.is_standard & !keys %in% standard_keys, , drop = FALSE]
    pairs <- annotation_report_pairs(d)
    on_curated <- pairs$feature_key %in% truth_features
    matched <- sum(pairs$pair_key %in% truth_pairs$pair_key)
    alternative <- sum(on_curated & !pairs$pair_key %in% truth_pairs$pair_key)
    missed <- sum(!truth_pairs$pair_key %in% pairs$pair_key)
    data.frame(polarity = report$polarity, phase = phase,
               curated_matches = matched, alternative_assignments = alternative,
               uncurated_assignments = sum(!on_curated), missed_curated_pairs = missed,
               precision = ratio(matched, matched + alternative),
               recall = ratio(matched, nrow(truth_pairs)),
               f1 = ratio(2 * matched, 2 * matched + alternative + missed))
  }))
}

annotation_metrics_tables <- function(positive, negative, merged) {
  reports <- list(positive, negative, merged)
  stopifnot(positive$polarity == "positive", negative$polarity == "negative",
            merged$polarity == "merged",
            identical(positive$study_id, negative$study_id),
            identical(positive$study_id, merged$study_id))
  tables <- list(Summary = do.call(rbind, lapply(reports, annotation_report_summary)),
                 Curated_reference_comparison = rbind(
                   annotation_reference_comparison(positive),
                   annotation_reference_comparison(negative)))
  tables$Summary <- tables$Summary[, c("polarity", "phase", "features", "removed_features")]
  tables$Curated_reference_comparison$alternative_assignments <- NULL
  detail_phases <- c("rank1_mz", "rank1_mz_rt", "isotope_filter", "adduct_scored",
                     "ambiguity_auto", "manual_curation", "qc_rsd_filtered",
                     "within_mode_adduct_resolution", "resolved_annotation_preference",
                     "duplicate_resolution", "internal_standard_removal")
  for (phase in detail_phases) {
    frames <- lapply(reports, function(report) report$phases[[phase]]$data)
    frames <- Filter(Negate(is.null), frames)
    if (!length(frames)) next
    cols <- unique(unlist(lapply(frames, names)))
    frames <- lapply(frames, function(d) {
      for (col in setdiff(cols, names(d))) d[[col]] <- rep(NA, nrow(d))
      d[, cols, drop = FALSE]
    })
    d <- do.call(rbind, frames)
    d$.is_standard <- NULL
    # Final details are one row per retained feature, matching the results XLSX.
    sheet <- if (phase == "internal_standard_removal") "final_annotations" else phase
    tables[[sheet]] <- d
  }
  tables
}

annotation_metrics_descriptions <- c(
  polarity = "Positive, negative, or merged. Feature IDs are qualified by polarity.",
  phase = "Actual workflow stage, in execution order; no alternative-method experiments.",
  features = "Distinct detected features remaining, including standards until their removal phase.",
  removed_features = "Features from the preceding stage absent here. Blank for the first row of each polarity.",
  curated_matches = "Candidate feature-lipid pairs matching the manually curated reference.",
  uncurated_assignments = "Assignments on features outside the curated reference; excluded from precision.",
  missed_curated_pairs = "Curated reference pairs not recovered by this stage.",
  precision = "Matches / (matches + alternative assignments), on curated features only.",
  recall = "Matches / all curated reference pairs.",
  f1 = "Harmonic mean of precision and recall.",
  feature_id = "Detected feature identifier; use polarity with this ID to identify a feature.",
  annotation = "Lipid annotation name, including RT-specific suffixes; unresolved alternatives separated by semicolons.",
  adduct = "Assigned adduct, or semicolon-separated alternatives.",
  mz = "Measured feature m/z.", rt_seconds = "Measured feature retention time, in seconds.",
  ppm_error = "Mass matching error, in ppm.", rt_error_seconds = "RT matching score/error, in seconds.",
  isotope_peaks = "Observed isotope peak count.", isotope_similarity = "Experimental/theoretical isotope similarity.",
  adduct_support = "Secondary-adduct support ratio; evidence, not a mandatory filter.",
  qc_cv = "QC RSD from IS-normalized abundances; not the raw-assay RSD used for filtering."
)

annotation_metrics_payload <- function(tables) {
  sheets <- lapply(names(tables), function(name) {
    d <- tables[[name]]
    descriptions <- unname(annotation_metrics_descriptions[names(d)])
    descriptions[is.na(descriptions)] <- ""
    list(name = name, columns = names(d), descriptions = descriptions,
         note = if (name == "Curated_reference_comparison")
           paste("Agreement with manually curated feature-lipid assignments on curated features,",
                 "not independent chemical ground truth or database accuracy.",
                 "Internal standards excluded. Uncurated assignments are not treated as wrong.",
                 "Manual curation agrees perfectly by construction.") else NULL,
         rows = unname(lapply(seq_len(nrow(d)), function(i) unname(as.list(d[i, ])))))
  })
  list(sheets = sheets)
}

# openxlsx can leave relationships to absent, unused drawing parts in otherwise
# plain workbooks. Remove only unreferenced drawing links so strict XLSX readers
# can open the report. Never discard a drawing actually referenced by a sheet.
clean_annotation_workbook_relationships <- function(path) {
  path <- normalizePath(path)
  directory <- tempfile("annotation-xlsx-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = directory)
  rels <- list.files(directory, pattern = "\\.rels$", recursive = TRUE, full.names = TRUE)
  changed <- FALSE
  for (rel in rels) {
    xml <- paste(readLines(rel, warn = FALSE), collapse = "\n")
    tags <- regmatches(xml, gregexpr("<Relationship\\b[^>]*?/>", xml, perl = TRUE))[[1]]
    for (tag in tags) {
      if (!grepl('/(drawing|vmlDrawing)"', tag)) next
      target <- sub('.* Target="([^"]+)".*', '\\1', tag)
      owner_dir <- dirname(dirname(rel))
      if (file.exists(file.path(owner_dir, target))) next
      id <- sub('.* Id="([^"]+)".*', '\\1', tag)
      owner <- file.path(owner_dir, sub("\\.rels$", "", basename(rel)))
      source <- paste(readLines(owner, warn = FALSE), collapse = "\n")
      if (grepl(paste0('r:id="', id, '"'), source, fixed = TRUE))
        stop("Export references a missing drawing part: ", target)
      xml <- sub(tag, "", xml, fixed = TRUE)
      changed <- TRUE
    }
    writeLines(xml, rel, useBytes = TRUE)
  }
  if (changed) zip::zipr(path, list.files(directory, recursive = TRUE, all.files = TRUE,
                                         no.. = TRUE), root = directory, mode = "mirror")
  invisible(path)
}

write_annotation_metrics_workbook <- function(positive, negative, merged, path) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Metrics export requires openxlsx. Install it with install.packages('openxlsx').")
  tables <- annotation_metrics_tables(positive, negative, merged)
  payload <- annotation_metrics_payload(tables)
  wb <- openxlsx::createWorkbook()
  header <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                                  fontColour = "#FFFFFF", fgFill = "#284B63", wrapText = TRUE)
  description <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
                                       fontColour = "#555555", wrapText = TRUE, valign = "top")
  for (sheet in payload$sheets) {
    name <- sheet$name
    d <- tables[[name]]
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    start <- if (is.null(sheet$note)) 1L else 3L
    if (!is.null(sheet$note)) {
      openxlsx::writeData(wb, name, sheet$note, startRow = 1, colNames = FALSE)
      openxlsx::mergeCells(wb, name, cols = seq_len(ncol(d)), rows = 1)
      openxlsx::addStyle(wb, name, description, rows = 1, cols = 1)
      openxlsx::setRowHeights(wb, name, 1, 60)
    }
    openxlsx::writeData(wb, name, t(sheet$columns), startRow = start, colNames = FALSE)
    openxlsx::writeData(wb, name, t(sheet$descriptions), startRow = start + 1L, colNames = FALSE)
    openxlsx::writeData(wb, name, d, startRow = start + 2L, colNames = FALSE, keepNA = FALSE)
    if (nrow(d)) openxlsx::addStyle(wb, name,
      openxlsx::createStyle(fontName = "Arial", fontSize = 10),
      rows = seq.int(start + 2L, start + 1L + nrow(d)),
      cols = seq_len(ncol(d)), gridExpand = TRUE)
    openxlsx::addStyle(wb, name, header, rows = start, cols = seq_len(ncol(d)), gridExpand = TRUE)
    openxlsx::addStyle(wb, name, description, rows = start + 1L,
                       cols = seq_len(ncol(d)), gridExpand = TRUE)
    widths <- if (name == "Summary") c(14, 30, 20, 24) else
      if (name == "Curated_reference_comparison") c(14, 30, rep(22, ncol(d) - 2L)) else
        c(16, 14, 44, rep(20, ncol(d) - 3L))
    openxlsx::setColWidths(wb, name, seq_len(ncol(d)), widths)
    openxlsx::setRowHeights(wb, name, start, 30)
    openxlsx::setRowHeights(wb, name, start + 1L, 90)
    openxlsx::freezePane(wb, name, firstActiveRow = start + 2L, firstActiveCol = 3)
    if (nrow(d)) {
      rates <- which(names(d) %in% c("precision", "recall", "f1"))
      if (length(rates)) openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = "0.0%"),
        rows = seq.int(start + 2L, start + 1L + nrow(d)), cols = rates, gridExpand = TRUE)
    }
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  clean_annotation_workbook_relationships(path)
  invisible(path)
}

regenerate_annotation_metrics <- function(study_dir, study_id) {
  add_rank1_mz_reporting(study_dir, study_id)
  read_report <- function(polarity, suffix) {
    path <- file.path(study_dir, polarity, "objects",
                      paste0(study_id, "_annotation_report_", suffix, ".rds"))
    if (!file.exists(path)) stop("Missing workflow snapshot: ", path,
                                ". Rerun the corresponding annotation workflow.")
    report <- readRDS(path)
    required <- c("preprocessed_features", "rank1_mz", "rank1_mz_rt", "isotope_filter", "adduct_scored",
                  "ambiguity_auto", "manual_curation", "normalization", "qc_rsd_filtered",
                  "within_mode_adduct_resolution", "qc_samples_removed")
    if (!identical(report$version, 1L) || !identical(report$study_id, study_id) ||
        !identical(report$polarity, polarity) || !identical(names(report$phases), required))
      stop("Incomplete or incompatible annotation snapshot: ", path)
    report
  }
  positive <- read_report("positive", "pos")
  negative <- read_report("negative", "neg")
  merged_path <- file.path(study_dir, "objects", paste0(study_id, "_annotation_report_merged.rds"))
  if (!file.exists(merged_path)) stop("Missing merge snapshot; rerun POS_NEG_merge.qmd.")
  merged <- readRDS(merged_path)
  if (!identical(merged$version, 1L) ||
      !identical(names(merged$phases), c("merge_input", "resolved_annotation_preference",
                                       "duplicate_resolution", "internal_standard_removal")))
    stop("Incomplete or incompatible merge snapshot: ", merged_path)
  expected <- rbind(positive$phases$qc_samples_removed$data,
                    negative$phases$qc_samples_removed$data)
  cols <- c("feature_id", "polarity", "annotation")
  if (!identical(expected[, cols], merged$phases$merge_input$data[, cols]))
    stop("Annotation snapshots do not match merge inputs; rerun POS_NEG_merge.qmd.")
  final <- merged$phases$internal_standard_removal$data
  if (any(final$.is_standard)) stop("Internal standards remain in final annotations.")
  write_annotation_metrics_workbook(positive, negative, merged,
    file.path(study_dir, "objects", paste0(study_id, "_annotation_metrics.xlsx")))
}
