suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lipid_helpers.R")
source("R/annotation_reporting.R")
chunk_code <- function(path, label) {
  lines <- readLines(path)
  start <- which(lines == paste0("#| label: ", label))
  stopifnot(length(start) == 1L)
  tail <- lines[seq.int(start + 1L, length(lines))]
  tail[seq_len(which(tail == "```")[1L] - 1L)]
}
fixture <- function(ids, lipid, adduct, rt, study, qc) {
  x <- SummarizedExperiment(assays = list(ISnorm_filled_imputed = cbind(study, qc)),
    colData = S4Vectors::DataFrame(sample_name = c("S1", "QC1"), sample_type = c("Plasma", "QC")),
    rowData = S4Vectors::DataFrame(target_lipid_name_unique = lipid, target_Adduct = adduct,
      target_lipid.name = lipid, target_IS_norm = rep("standard", length(ids)), rtmed = rt))
  rownames(x) <- ids
  x
}
annotated_res <- fixture(c("P1", "P2", "P3", "P4"),
  c("A", "A", "B; C", "D"), c("H", "Na", "H", "H"), c(10, 12, 20, 30),
  c(10, 5, 10, 10), c(1, 1000, 1, 1))
POLARITY <- "pos"; ADDUCT_RT_THRESHOLD <- 4
annotation_report <- new_annotation_report("TEST", "positive")
eval(parse(text = chunk_code("applications/MICROSAMPLING_study/positive/Annotation_pos.qmd",
                             "resolve-within-mode-adducts")))
stopifnot(identical(rownames(annotated_res), c("P1", "P3", "P4")),
          annotation_report_summary(annotation_report)$features == 3L)
res_pos <- annotated_res[, 1L, drop = FALSE]
res_neg <- fixture(c("N1", "N2"), c("B", "D"), c("H", "H"), c(999, 1000),
                   c(10, 100), c(1, 1))[, 1L, drop = FALSE]
merge_report <- new_annotation_report("TEST", "merged")
for (label in c("prefer-resolved-sn-position", "find-duplicates", "remove-duplicates")) {
  eval(parse(text = chunk_code("applications/MICROSAMPLING_study/POS_NEG_merge.qmd", label)))
}
stopifnot(identical(rownames(res_pos_filtered), "P1"),
          identical(rownames(res_neg_filtered), c("N1", "N2")),
          nrow(merge_report$phases$resolved_annotation_preference$data) == 4L,
          comparison_df$rt_diff == 970, comparison_df$keep_mode == "negative")
cat("PASS: within-mode dedup excludes QC, resolved annotations replace ambiguous rows, and cross-mode dedup ignores RT eligibility\n")
