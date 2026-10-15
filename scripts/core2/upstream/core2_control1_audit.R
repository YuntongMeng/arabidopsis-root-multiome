# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Focused Control 1 audit. User-run only; no filtering, alignment or overwriting QC.
# A: compare published ATAC-only counts to rebuilt counts by identical barcode,
#    using all four samples as a baseline (units/peak sets differ).
# B: independently rescan Control 1 with the exact fragment-QC function used for
#    the other three samples, and compare to its previously imported QC columns.

run_core2_control1_audit <- function(
    task_dir = Sys.getenv("P3_CORE2_WORK", ""),
    project_dir = Sys.getenv("P3_ROOT", getwd())) {
  .libPaths(c(file.path(task_dir, "work", "R_library"), .libPaths()))
  stopifnot(requireNamespace("Matrix", quietly = TRUE),
            requireNamespace("GenomicRanges", quietly = TRUE),
            requireNamespace("IRanges", quietly = TRUE))
  qc_dir <- file.path(task_dir, "outputs", "core2_qc")
  output_dir <- file.path(qc_dir, "control1_audit")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  qc <- utils::read.csv(file.path(qc_dir, "all_nuclei_qc.csv"), check.names = FALSE)
  sample_ids <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
  expected_n <- c(3431L, 6631L, 7827L, 6585L)
  tables <- list()
  summaries <- list()
  unique_file <- function(path, pattern) {
    candidates <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(candidates) != 1L) stop("Expected one file matching ", pattern, " in ", path)
    candidates[[1]]
  }
  read_published <- function(sample_id) {
    folder <- file.path(project_dir, "data", "rna", sample_id)
    feature_path <- unique_file(folder, "features\\.tsv\\.gz$")
    barcode_path <- unique_file(folder, "barcodes\\.tsv\\.gz$")
    matrix_path <- unique_file(folder, "matrix\\.mtx\\.gz$")
    input <- gzfile(matrix_path, open = "rt")
    on.exit(close(input), add = TRUE)
    features <- utils::read.table(
      feature_path, sep = "\t", header = FALSE, quote = "", comment.char = "",
      stringsAsFactors = FALSE
    )
    barcode_input <- gzfile(barcode_path, open = "rt")
    on.exit(close(barcode_input), add = TRUE)
    barcodes <- readLines(barcode_input, warn = FALSE)
    mtx <- Matrix::readMM(input)
    stopifnot(nrow(mtx) == nrow(features), ncol(mtx) == length(barcodes),
              !anyDuplicated(barcodes),
              all(features[[3]] %in% c("Peaks", "Gene Expression")))
    peak_rows <- which(features[[3]] == "Peaks")
    stopifnot(length(peak_rows) > 0L)
    # Crucial: gene-expression rows never contribute to this ATAC comparison.
    list(barcode = barcodes, counts = as.numeric(Matrix::colSums(
      mtx[peak_rows, , drop = FALSE]
    )), peaks = length(peak_rows))
  }

  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    message("Audit published versus rebuilt ATAC counts: ", sample_id)
    published <- read_published(sample_id)
    d <- qc[qc$sample == sample_id, , drop = FALSE]
    stopifnot(nrow(d) == expected_n[i], length(published$barcode) == expected_n[i],
              !anyDuplicated(d$barcode), setequal(d$barcode, published$barcode))
    index <- match(d$barcode, published$barcode)
    stopifnot(!anyNA(index))
    d$published_peak_counts <- published$counts[index]
    stopifnot(all(is.finite(d$published_peak_counts)))
    tables[[sample_id]] <- d
    summaries[[sample_id]] <- data.frame(
      sample = sample_id, nuclei = nrow(d), published_peaks = published$peaks,
      median_published_peak_counts = stats::median(d$published_peak_counts),
      median_rebuilt_peak_counts = stats::median(d$nCount_ATAC),
      spearman_published_vs_rebuilt_peaks = stats::cor(
        d$published_peak_counts, d$nCount_ATAC, method = "spearman"
      ),
      spearman_published_vs_nuclear_fragments = stats::cor(
        d$published_peak_counts, d$nuclear_fragments, method = "spearman"
      ),
      all_published_barcodes_matched = TRUE
    )
    utils::write.csv(d, file.path(output_dir, paste0(sample_id, "_paired_counts.csv")),
                     row.names = FALSE)
    rm(published, d)
    invisible(gc())
  }
  concordance <- do.call(rbind, summaries)
  rownames(concordance) <- NULL
  utils::write.csv(concordance, file.path(output_dir, "published_vs_rebuilt_summary.csv"),
                   row.names = FALSE)
  # Descriptive strata only. No nuclei are removed and no cutoff is finalized.
  c1 <- tables[["control_rep1"]]
  group <- ifelse(c1$nuclear_fragments < 500, "nuclear_fragments_below_500",
                  "nuclear_fragments_at_least_500")
  strata <- do.call(rbind, lapply(split(c1, group), function(d) {
    data.frame(
      group = if (d$nuclear_fragments[1] < 500) "below_500" else "at_least_500",
      nuclei = nrow(d), median_published_peak_counts = stats::median(d$published_peak_counts),
      median_rebuilt_peak_counts = stats::median(d$nCount_ATAC),
      median_nuclear_fragments = stats::median(d$nuclear_fragments),
      median_organelle_pct = stats::median(d$percent_organelle)
    )
  }))
  rownames(strata) <- NULL
  utils::write.csv(strata, file.path(output_dir, "control1_descriptive_strata.csv"),
                   row.names = FALSE)
  local({
    grDevices::pdf(file.path(output_dir, "published_vs_rebuilt_plots.pdf"),
                   width = 10, height = 9)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    xlim <- range(unlist(lapply(tables, function(d) log10(d$published_peak_counts + 1))))
    ylim <- range(unlist(lapply(tables, function(d) log10(d$nCount_ATAC + 1))))
    for (sample_id in sample_ids) {
      d <- tables[[sample_id]]
      graphics::plot(
        log10(d$published_peak_counts + 1), log10(d$nCount_ATAC + 1),
        pch = 16, cex = .35,
        col = grDevices::adjustcolor("steelblue", alpha.f = .35),
        xlim = xlim, ylim = ylim, main = sample_id,
        xlab = "log10(published ATAC peak counts + 1)",
        ylab = "log10(rebuilt common-peak counts + 1)"
      )
    }
  })

  # Load helpers without triggering run_core2_batch_qc(). Restore the switch.
  helpers <- new.env(parent = .GlobalEnv)
  local({
    old <- Sys.getenv("CORE2_QC_DEFINE_ONLY", unset = NA_character_)
    on.exit({
      if (is.na(old)) Sys.unsetenv("CORE2_QC_DEFINE_ONLY")
      else Sys.setenv(CORE2_QC_DEFINE_ONLY = old)
    }, add = TRUE)
    Sys.setenv(CORE2_QC_DEFINE_ONLY = "1")
    sys.source(file.path(Sys.getenv("P3_ROOT", getwd()), "scripts", "core2", "upstream", "core2_batch_qc.R"), envir = helpers)
  })
  message("Rescan Control 1 with the SAME independent-fragment QC function as the other samples")
  peaks <- readRDS(file.path(task_dir, "outputs", "core2_common_peaks.rds"))
  fresh <- helpers$core2_qc_scan_fragments(
    path = file.path(task_dir, "outputs", "control1_upstream", "control1.fragments.tsv.gz"),
    peaks = peaks, cells = c1$barcode
  )
  stopifnot(identical(as.vector(fresh$barcode), as.vector(c1$barcode)))
  metrics <- c("total_fragments", "nuclear_fragments", "fragments_in_common_peaks",
               "FRiP_nuclear", "percent_organelle")
  comparison <- do.call(rbind, lapply(metrics, function(metric) {
    a <- fresh[[metric]]
    b <- c1[[metric]]
    same_missing <- is.na(a) & is.na(b)
    finite <- is.finite(a) & is.finite(b)
    tolerance <- if (metric %in% c("FRiP_nuclear", "percent_organelle")) 1e-8 else 0
    agree <- same_missing | (finite & abs(a - b) <= tolerance)
    agree[is.na(agree)] <- FALSE
    data.frame(
      metric = metric, nuclei_compared = length(a), nuclei_disagree = sum(!agree),
      max_absolute_difference = if (any(finite)) max(abs(a[finite] - b[finite])) else NA_real_,
      within_tolerance = all(agree)
    )
  }))
  utils::write.csv(fresh, file.path(output_dir, "control1_fresh_fragment_qc.csv"),
                   row.names = FALSE)
  utils::write.csv(comparison, file.path(output_dir, "control1_old_vs_fresh_qc.csv"),
                   row.names = FALSE)
  utils::capture.output(utils::sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
  cat("\nPublished versus rebuilt counts (different counting units; no equality requirement):\n")
  print(concordance, row.names = FALSE)
  cat("\nControl 1 descriptive groups, without filtering:\n")
  print(strata, row.names = FALSE)
  cat("\nControl 1 imported QC versus fresh same-function QC:\n")
  print(comparison, row.names = FALSE)
  message("AUDIT COMPLETE: ", output_dir)
  if (!all(comparison$within_tolerance)) {
    warning("Control 1 imported QC differs from fresh calculation. Do not use the old values for filtering.")
  }
  invisible(list(concordance = concordance, strata = strata, comparison = comparison))
}

core2_control1_audit_result <- run_core2_control1_audit()
