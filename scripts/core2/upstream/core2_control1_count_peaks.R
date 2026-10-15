# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Recount common peaks from Control 1 fragments; no nucleus QC or normalization.
# Usage: Rscript core2_control1_count_peaks.R /absolute/path/to/task
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
task_dir <- normalizePath(args[[1]])
.libPaths(c(file.path(task_dir, "work", "R_library"), .libPaths()))
stopifnot(requireNamespace("GenomicRanges", quietly = TRUE),
          requireNamespace("Signac", quietly = TRUE),
          requireNamespace("Matrix", quietly = TRUE),
          requireNamespace("Rsamtools", quietly = TRUE))

output_dir <- file.path(task_dir, "outputs", "control1_upstream")
fragment_path <- file.path(output_dir, "control1.fragments.tsv.gz")
common_peaks <- readRDS(file.path(task_dir, "outputs", "core2_common_peaks.rds"))
cells <- readLines(file.path(output_dir, "matrix_barcodes.tsv"))
stopifnot(length(common_peaks) == 27677L, length(cells) == 3431L,
          !anyDuplicated(cells), GenomicRanges::isDisjoint(common_peaks))

# Validate the counting unit and boundary behavior with a tiny known fixture.
# Multiplicity 7 must contribute one unique fragment, not seven PCR copies.
fixture_dir <- file.path(task_dir, "work", "control1_upstream", "count_fixture")
dir.create(fixture_dir, showWarnings = FALSE)
fixture <- file.path(fixture_dir, "fragments.tsv")
writeLines(c(
  paste("1", 100, 120, cells[[1]], 7, sep = "\t"),
  paste("1", 110, 130, cells[[1]], 3, sep = "\t"),
  paste("1", 140, 150, cells[[2]], 1, sep = "\t")
), fixture)
fixture_gz <- Rsamtools::bgzip(fixture, overwrite = TRUE)
Rsamtools::indexTabix(fixture_gz, format = "bed")
fixture_features <- GenomicRanges::GRanges(
  seqnames = rep("1", 4),
  ranges = IRanges::IRanges(start = c(101L, 121L, 131L, 141L),
                           end = c(110L, 130L, 140L, 150L))
)
fixture_object <- Signac::CreateFragmentObject(
  path = fixture_gz, validate.fragments = FALSE, verbose = FALSE
)
fixture_counts <- Signac::FeatureMatrix(
  fragments = fixture_object, features = fixture_features,
  cells = cells[1:2], keep_all_features = TRUE, verbose = FALSE
)
expected_fixture <- matrix(c(1, 1, 0, 0, 0, 0, 0, 1), nrow = 4L)
stopifnot(identical(unname(as.matrix(fixture_counts)), expected_fixture))
fixture_chunked <- Signac::FeatureMatrix(
  fragments = fixture_object, features = fixture_features,
  cells = cells[1:2], keep_all_features = FALSE, process_n = 1L,
  verbose = FALSE
)
stopifnot(identical(unname(as.matrix(fixture_chunked)), expected_fixture))
# Record the behavior of an unused option; production uses the tested FALSE.
fixture_true_chunked <- Signac::FeatureMatrix(
  fragments = fixture_object, features = fixture_features,
  cells = cells[1:2], keep_all_features = TRUE, process_n = 1L,
  verbose = FALSE
)
capture.output(list(
  expected = expected_fixture,
  production_setting_FALSE = unname(as.matrix(fixture_chunked)),
  unused_setting_TRUE = unname(as.matrix(fixture_true_chunked)),
  production_setting_passed = TRUE,
  unused_TRUE_setting_matches_expected = identical(
    unname(as.matrix(fixture_true_chunked)), expected_fixture
  )
), file = file.path(output_dir, "logs", "counting_fixture_diagnostic.txt"))
if (length(args) == 2L && identical(args[[2]], "--fixture-only")) {
  message("Counting unit, half-open boundaries, zero rows and chunking fixture passed")
  quit(status = 0L)
}

# Require real coverage of all five nuclear contigs before retaining zero rows.
tabix <- Rsamtools::TabixFile(fragment_path)
indexed_contigs <- Rsamtools::seqnamesTabix(tabix)
stopifnot(all(as.character(1:5) %in% indexed_contigs))
fragment <- Signac::CreateFragmentObject(
  path = fragment_path, validate.fragments = FALSE, verbose = TRUE
)
counts <- Signac::FeatureMatrix(
  fragments = fragment, features = common_peaks,
  # All five chromosomes are proven present above, so FALSE still keeps every
  # input feature. In Signac 1.17.1, TRUE can fill every chunk with zero rows for
  # other chunks, creating duplicate row names and selecting erroneous zeros.
  cells = cells, keep_all_features = FALSE, process_n = 500L,
  verbose = TRUE
)
expected_rows <- Signac::GRangesToString(common_peaks)
stopifnot(nrow(counts) == length(common_peaks), ncol(counts) == length(cells),
          identical(rownames(counts), expected_rows),
          setequal(colnames(counts), cells),
          !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)))
counts <- counts[, cells, drop = FALSE]
stopifnot(identical(colnames(counts), cells),
          all(is.finite(counts@x)), all(counts@x >= 0),
          all(counts@x == floor(counts@x)))

# Independent recount from raw BED fields for selected real peak rows.
# Explicit 0-based half-open overlap avoids silently accepting an off-by-one.
chr <- as.character(GenomicRanges::seqnames(common_peaks))
set.seed(20261010)
check_rows <- sort(unique(c(
  sample.int(length(common_peaks), 30L),
  which(Matrix::rowSums(counts) > 0)[seq_len(10L)],
  unlist(lapply(split(seq_along(chr), chr), function(x) c(x[[1]], tail(x, 1L))))
)))
stopifnot(!anyNA(check_rows))
open(tabix)
checks <- lapply(check_rows, function(i) {
  region <- common_peaks[i]
  lines <- Rsamtools::scanTabix(tabix, param = region)[[1]]
  expected <- integer(length(cells))
  if (length(lines)) {
    fields <- strsplit(lines, "\t", fixed = TRUE)
    keep <- vapply(fields, function(x) {
      identical(x[[1]], chr[[i]]) &&
        as.integer(x[[2]]) < GenomicRanges::end(region) &&
        as.integer(x[[3]]) > GenomicRanges::start(region) - 1L
    }, logical(1))
    ids <- vapply(fields[keep], `[[`, character(1), 4L)
    column <- match(ids, cells)
    expected <- tabulate(column[!is.na(column)], nbins = length(cells))
  }
  observed <- as.numeric(counts[i, ])
  stopifnot(identical(as.numeric(expected), observed))
  data.frame(peak = expected_rows[[i]], nuclei_checked = length(cells),
             fragment_count_sum = sum(expected), passed = TRUE)
})
close(tabix)
utils::write.csv(do.call(rbind, checks),
                 file.path(output_dir, "independent_recount_checks.csv"), row.names = FALSE)

# Keep all supplied matrix barcodes and all common intervals, including zeros.
saveRDS(counts, file.path(output_dir, "control1_common_peak_counts.rds"))
saveRDS(common_peaks, file.path(output_dir, "common_peaks.rds"))
utils::write.table(
  data.frame(peak_id = expected_rows, chrom = chr,
             start_0based = GenomicRanges::start(common_peaks) - 1L,
             end_0based_exclusive = GenomicRanges::end(common_peaks)),
  file.path(output_dir, "common_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
manifest <- data.frame(
  sample = "control_rep1", run = "SRR24921616", peaks = nrow(counts),
  matrix_nuclei_retained = ncol(counts), nonzero_entries = length(counts@x),
  zero_count_nuclei_retained = sum(Matrix::colSums(counts) == 0),
  zero_count_peak_rows_retained = sum(Matrix::rowSums(counts) == 0),
  count_unit = "unique deduplicated fragments overlapping each peak; PCR multiplicity ignored",
  nucleus_QC_applied = FALSE, fixture_test_passed = TRUE,
  independent_peak_rows_checked = length(check_rows),
  Signac_version = as.character(utils::packageVersion("Signac")),
  Matrix_version = as.character(utils::packageVersion("Matrix"))
)
utils::write.csv(manifest, file.path(output_dir, "count_matrix_manifest.csv"), row.names = FALSE)
capture.output(utils::sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
print(manifest, row.names = FALSE)
