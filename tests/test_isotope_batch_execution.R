# Verify the execution settings and row/spectrum alignment at the batching boundary.
suppressPackageStartupMessages(library(Spectra))
suppressPackageStartupMessages(library(BiocParallel))
source("R/lipid_helpers.R")

input <- data.frame(feature_id = c("a", "b", "a", "c", "b"), value = 1:5)
.calculate_isotope_similarity_batch <- function(mse, mtched_data, polarity,
                                                isopeak_threshold,
                                                similarity_threshold, BPPARAM) {
  stopifnot(inherits(BPPARAM, "SerialParam"))
  mtched_data <- mtched_data[mtched_data$value != 2L, , drop = FALSE]
  mtched_data$isopeak_count <- rep(2L, nrow(mtched_data))
  mtched_data$isopeak_sim <- rep(1, nrow(mtched_data))
  d <- data.frame(msLevel = rep(1L, nrow(mtched_data)),
                  feature_id = mtched_data$feature_id, value = mtched_data$value)
  d$mz <- lapply(mtched_data$value, function(x) c(100, 101))
  d$intensity <- lapply(mtched_data$value, function(x) c(x, x / 2))
  sp <- scalePeaks(Spectra(d), by = max)
  list(mtched_data = mtched_data, iso_spectra = sp, theoretical_spectra = sp)
}

register(MulticoreParam(2L))
previous <- bpparam()
single <- calculate_isotope_similarity(NULL, input, batch_size = Inf)
stopifnot(identical(bpparam(), previous))
batched <- calculate_isotope_similarity(NULL, input, batch_size = 1L)
stopifnot(identical(rownames(single$mtched_data), rownames(batched$mtched_data)),
          all(vapply(names(single$mtched_data), function(name)
            identical(single$mtched_data[[name]], batched$mtched_data[[name]]), logical(1))),
          identical(batched$mtched_data$value, c(1L, 3L, 4L, 5L)),
          identical(rownames(batched$mtched_data), c("1", "3", "4", "5")),
          identical(batched$iso_spectra$value, batched$mtched_data$value),
          identical(batched$theoretical_spectra$value, batched$mtched_data$value),
          identical(bpparam(), previous),
          all(vapply(peaksData(batched$iso_spectra), function(x)
            identical(unname(x[, 2L]), c(1, 0.5)), logical(1))))

with_empty_batch <- input[1:4, ]
with_empty_batch$feature_id <- c("a", "b", "c", "a")
empty_batch <- calculate_isotope_similarity(NULL, with_empty_batch, batch_size = 1L)
stopifnot(identical(empty_batch$mtched_data$value, c(1L, 3L, 4L)),
          identical(empty_batch$iso_spectra$value, c(1L, 3L, 4L)))
empty <- calculate_isotope_similarity(NULL, input[FALSE, ], batch_size = 1L)
stopifnot(nrow(empty$mtched_data) == 0L, length(empty$iso_spectra) == 0L)

.calculate_isotope_similarity_batch <- function(...) stop("intentional failure")
stopifnot(inherits(try(calculate_isotope_similarity(NULL, input, batch_size = 1L),
                       silent = TRUE), "try-error"),
          identical(bpparam(), previous))
cat("PASS: serial isotope batches, row/spectrum alignment, empty batches, and registry restoration.\n")
