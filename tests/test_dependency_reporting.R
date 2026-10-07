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
# Lipid ambiguities now span adducts and select closest RT, with isotope ties.
candidates <- data.frame(feature_id = c("P1", "P2", "P3", "P4", "P5"),
  target_lipid_name_unique = c("A", "A", "B", "B", "C"),
  target_adduct = c("H", "Na", "H", "Na", "H"),
  score_rt = c(2, 1, -1, 1, 0), isopeak_sim = c(1, .9, .8, .95, 1))
folder <- tempfile(); dir.create(folder)
ambiguities <- export_ambiguity_tables(candidates, output_dir = folder, verbose = FALSE)
unlink(folder, recursive = TRUE)
stopifnot(identical(ambiguities$lipid_ambiguities$feature_id, c("P1", "P2", "P3", "P4")),
          identical(ambiguities$lipid_ambiguities$keep_row, c(FALSE, TRUE, FALSE, TRUE)))
# Only genuine sn-1/sn-2 pairs are replaced across modes; other ambiguities survive.
pair <- "LPC(0:0/18:1)_1; LPC(18:1/0:0)_1"
stopifnot(length(sn_pair_options(pair)) == 2L,
          length(sn_pair_options("A; B")) == 0L,
          length(sn_pair_options("LPC(0:0/18:1)_1; LPC(18:2/0:0)_1")) == 0L,
          length(sn_pair_options(NA_character_)) == 0L)
res_pos <- fixture(c("P1", "P3", "P4", "P5"),
  c("A", "B; C", "D", pair), rep("H", 4), c(10, 20, 30, 40),
  rep(10, 4), rep(1, 4))[, 1L, drop = FALSE]
res_neg <- fixture(c("N1", "N2", "N3"), c("B", "D", "LPC(18:1/0:0)_1"),
  rep("H", 3), c(999, 1000, 1100), c(10, 100, 10), rep(1, 3))[, 1L, drop = FALSE]
merge_report <- new_annotation_report("TEST", "merged")
for (label in c("prefer-resolved-sn-position", "find-duplicates", "remove-duplicates")) {
  eval(parse(text = chunk_code("applications/MICROSAMPLING_study/POS_NEG_merge.qmd", label)))
}
stopifnot(identical(rownames(res_pos_filtered), c("P1", "P3")),
          identical(rownames(res_neg_filtered), c("N1", "N2", "N3")),
          nrow(merge_report$phases$resolved_annotation_preference$data) == 6L,
          comparison_df$rt_diff == 970, comparison_df$keep_mode == "negative")
cat("PASS: ambiguity resolution spans adducts, only sn-isomer pairs are replaced, and cross-mode dedup ignores RT eligibility\n")
