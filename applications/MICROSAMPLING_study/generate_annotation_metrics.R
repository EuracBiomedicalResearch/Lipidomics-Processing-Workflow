#!/usr/bin/env Rscript
# Refresh metrics from workflow snapshots, including the mass-only reporting phase.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this file using Rscript.")
study_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(study_dir, "../../R/annotation_reporting.R"))
regenerate_annotation_metrics(study_dir, "MICROSAMPLING")
