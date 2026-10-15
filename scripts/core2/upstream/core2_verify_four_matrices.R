# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Fresh-session verification of the four pre-QC matrices; no biological analysis.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task_dir <- normalizePath(args[[1]])
.libPaths(c(file.path(task_dir, "work", "R_library"), .libPaths()))
stopifnot(requireNamespace("GenomicRanges", quietly = TRUE),
          requireNamespace("Signac", quietly = TRUE),
          requireNamespace("Matrix", quietly = TRUE))
project_dir <- Sys.getenv("P3_ROOT", getwd())
samples <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
tags <- c("control1", "control2", "stress1", "stress2")
expected_nuclei <- c(3431L, 6631L, 7827L, 6585L)
common <- readRDS(file.path(task_dir, "outputs", "core2_common_peaks.rds"))
expected_rows <- Signac::GRangesToString(common)
stopifnot(length(common) == 27677L, GenomicRanges::isDisjoint(common))

summary <- lapply(seq_along(samples), function(i) {
  output_dir <- file.path(task_dir, "outputs", paste0(tags[[i]], "_upstream"))
  stopifnot(file.exists(file.path(output_dir, "upstream_completion.json")))
  source_barcode_path <- list.files(
    file.path(project_dir, "data", "rna", samples[[i]]),
    pattern = "_barcodes\\.tsv\\.gz$", full.names = TRUE
  )
  stopifnot(length(source_barcode_path) == 1L)
  original_cells <- readLines(gzfile(source_barcode_path))
  saved_cells <- readLines(file.path(output_dir, "matrix_barcodes.tsv"))
  counts <- readRDS(file.path(output_dir, paste0(tags[[i]], "_common_peak_counts.rds")))
  peaks <- readRDS(file.path(output_dir, "common_peaks.rds"))
  stopifnot(identical(peaks, common), identical(saved_cells, original_cells),
            identical(rownames(counts), expected_rows),
            identical(colnames(counts), original_cells),
            nrow(counts) == 27677L, ncol(counts) == expected_nuclei[[i]],
            !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
            all(is.finite(counts@x)), all(counts@x >= 0),
            all(counts@x == floor(counts@x)))
  checks <- utils::read.csv(file.path(output_dir, "independent_recount_checks.csv"))
  stopifnot(nrow(checks) >= 40L, all(checks$passed),
            all(checks$nuclei_checked == expected_nuclei[[i]]))
  manifest <- utils::read.csv(file.path(output_dir, "count_matrix_manifest.csv"))
  stopifnot(nrow(manifest) == 1L, identical(manifest$sample, samples[[i]]),
            identical(manifest$nucleus_QC_applied, FALSE))
  data.frame(sample = samples[[i]], run = manifest$run,
             peaks = nrow(counts), nuclei_before_QC = ncol(counts),
             nonzero_entries = length(counts@x),
             matrix_fragment_overlap_sum = sum(counts@x),
             zero_count_nuclei_retained = sum(Matrix::colSums(counts) == 0),
             zero_count_peak_rows_retained = sum(Matrix::rowSums(counts) == 0),
             identical_common_peak_coordinates_and_order = TRUE,
             original_barcode_order_preserved = TRUE,
             independent_peak_rows_checked = nrow(checks),
             nucleus_QC_applied = FALSE)
})
summary <- do.call(rbind, summary)
utils::write.csv(summary, file.path(task_dir, "outputs", "core2_four_sample_matrix_checks.csv"),
                 row.names = FALSE)
print(summary, row.names = FALSE)
capture.output(utils::sessionInfo(),
               file = file.path(task_dir, "outputs", "core2_four_sample_audit_sessionInfo.txt"))
