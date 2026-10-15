# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Core 2 step 1: define common intervals; this does not generate count data.
# Cell Ranger ARC feature coordinates use zero-based, half-open BED format.
project_dir <- Sys.getenv("P3_ROOT", getwd())
samples <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
output_dir <- file.path(Sys.getenv("P3_CORE2_WORK"), "outputs")
stopifnot(nzchar(Sys.getenv("P3_CORE2_WORK")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

peak_sets <- setNames(lapply(samples, function(sample) {
  feature_file <- list.files(
    file.path(project_dir, "data", "rna", sample),
    pattern = "_features\\.tsv\\.gz$", full.names = TRUE
  )
  stopifnot(length(feature_file) == 1L)
  features <- utils::read.delim(feature_file, header = FALSE,
                              stringsAsFactors = FALSE)
  peaks <- features[features[[3]] == "Peaks", , drop = FALSE]
  stopifnot(nrow(peaks) > 0L, !anyDuplicated(peaks[[1]]),
            all(peaks[[5]] >= 0L), all(peaks[[6]] > peaks[[5]]),
            all(peaks[[4]] %in% as.character(1:5)),
            identical(peaks[[1]], paste0(peaks[[4]], ":", peaks[[5]], "-", peaks[[6]])))
  peaks[[1]]
}), samples)

peak_ranges <- lapply(peak_sets, function(ids) {
  coordinates <- utils::strcapture(
    "^([^:]+):([0-9]+)-([0-9]+)$", ids,
    proto = list(chr = character(), start = integer(), end = integer())
  )
  stopifnot(!anyNA(coordinates))
  GenomicRanges::makeGRangesFromDataFrame(
    coordinates, starts.in.df.are.0based = TRUE
  )
})
common_peaks_unfiltered <- GenomicRanges::reduce(
  do.call(c, unname(peak_ranges))
)
region_width <- IRanges::width(common_peaks_unfiltered)
common_peaks <- common_peaks_unfiltered[
  region_width > 20L & region_width < 10000L
]
stopifnot(GenomicRanges::isDisjoint(common_peaks))

saveRDS(common_peaks, file.path(output_dir, "core2_common_peaks.rds"))
bed <- as.data.frame(common_peaks)[, c("seqnames", "start", "end")]
bed$start <- bed$start - 1L
utils::write.table(bed, file.path(output_dir, "core2_common_peaks.bed"),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
manifest <- data.frame(
  sample = samples, original_peaks = lengths(peak_sets),
  common_intervals_before_width_filter = length(common_peaks_unfiltered),
  common_intervals_after_width_filter = length(common_peaks)
)
utils::write.csv(manifest, file.path(output_dir, "core2_common_peak_manifest.csv"),
                 row.names = FALSE)
print(manifest, row.names = FALSE)
print(utils::head(common_peaks))
print(summary(IRanges::width(common_peaks)))
