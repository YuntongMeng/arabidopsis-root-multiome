# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Project 3 / Core 2: calculate comparable QC metrics for four ATAC samples.
# Source defines helpers; call run_core2_batch_qc(task_dir = Sys.getenv("P3_CORE2_WORK")) explicitly.
# No FASTQ processing, peak calling, normalization, cell filtering, or doublet calls.
# Existing original fragments and count matrices are never overwritten.

core2_qc_log <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
}

core2_qc_save <- function(x, path) {
  tmp <- tempfile(pattern = ".qc_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish checkpoint: ", path)
  invisible(path)
}

core2_qc_file_signature <- function(path, hash = FALSE) {
  path <- normalizePath(path, mustWork = TRUE)
  info <- file.info(path)
  list(path = path, size = unname(info$size),
       modified = as.numeric(info$mtime),
       md5 = if (hash) unname(tools::md5sum(path)) else NULL)
}

core2_qc_ranges_key <- function(x, strand = FALSE) {
  key <- paste0(as.character(GenomicRanges::seqnames(x)), ":",
                GenomicRanges::start(x), "-", GenomicRanges::end(x))
  if (strand) key <- paste0(key, ":", as.character(GenomicRanges::strand(x)))
  key
}

# One pass through the ORIGINAL fragments file. Each deduplicated row counts once;
# column 5 (PCR support) is deliberately not a weight. BED start is converted +1.
# If requested, also write a nuclear-only copy, retaining all original barcodes.
core2_qc_scan_fragments <- function(path, peaks, cells, nuclear_plain = NULL,
                                    chunk_size = 100000L, verbose = TRUE) {
  stopifnot(!anyDuplicated(cells), length(cells) > 0L)
  input <- gzfile(path, open = "rt")
  on.exit(close(input), add = TRUE)
  if (!is.null(nuclear_plain)) {
    output <- file(nuclear_plain, open = "wt")
    on.exit(close(output), add = TRUE)
  }
  total <- nuclear <- hits <- organelle <- numeric(length(cells))
  peak_index <- GenomicRanges::GNCList(peaks)
  scanned <- 0
  next_report <- 5000000

  repeat {
    lines <- readLines(input, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    # Our upstream files have five columns and no header/comment records.
    if (any(lengths(regmatches(lines, gregexpr("\t", lines, fixed = TRUE))) != 4L)) {
      stop("Expected exactly five tab-separated columns in: ", path)
    }
    fragments <- utils::read.table(
      text = lines, sep = "\t", header = FALSE, quote = "", comment.char = "",
      colClasses = c("character", "integer", "integer", "character", "NULL")
    )
    if (anyNA(fragments) || any(fragments[[2]] < 0L) ||
        any(fragments[[3]] <= fragments[[2]])) {
      stop("Invalid fragment coordinates/barcodes in: ", path)
    }
    is_nuclear <- fragments[[1]] %in% as.character(1:5)
    is_organelle <- fragments[[1]] %in% c("Mt", "Pt")
    if (any(!(is_nuclear | is_organelle))) {
      stop("Unexpected contig: require only 1-5, Mt and Pt in: ", path)
    }
    if (!is.null(nuclear_plain)) writeLines(lines[is_nuclear], output)

    cell_index <- match(fragments[[4]], cells)
    keep <- !is.na(cell_index)
    total <- total + tabulate(cell_index[keep], nbins = length(cells))
    organelle <- organelle + tabulate(
      cell_index[keep & is_organelle], nbins = length(cells)
    )
    keep_nuclear <- keep & is_nuclear
    nuclear <- nuclear + tabulate(
      cell_index[keep_nuclear], nbins = length(cells)
    )
    if (any(keep_nuclear)) {
      fragment_ranges <- GenomicRanges::GRanges(
        seqnames = fragments[[1]][keep_nuclear],
        ranges = IRanges::IRanges(
          start = fragments[[2]][keep_nuclear] + 1L,
          end = fragments[[3]][keep_nuclear]
        )
      )
      in_peak <- IRanges::overlapsAny(
        query = fragment_ranges, subject = peak_index, ignore.strand = TRUE
      )
      hits <- hits + tabulate(
        cell_index[keep_nuclear][in_peak], nbins = length(cells)
      )
    }
    scanned <- scanned + length(lines)
    if (verbose && scanned >= next_report) {
      core2_qc_log("Scanned ", format(scanned, big.mark = ","), " fragment rows")
      next_report <- scanned + 5000000
    }
  }
  stopifnot(all(total == nuclear + organelle), all(hits <= nuclear))
  frip <- percent_organelle <- rep(NA_real_, length(cells))
  frip[nuclear > 0] <- hits[nuclear > 0] / nuclear[nuclear > 0]
  percent_organelle[total > 0] <-
    100 * organelle[total > 0] / total[total > 0]
  data.frame(
    barcode = cells, total_fragments = total, nuclear_fragments = nuclear,
    organelle_fragments = organelle, fragments_in_common_peaks = hits,
    FRiP_nuclear = frip, percent_organelle = percent_organelle,
    stringsAsFactors = FALSE
  )
}

core2_qc_fragment_stage <- function(original, nuclear_path, peaks, cells,
                                    signature, cache_path) {
  if (file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    stopifnot(identical(cache$signature, signature),
              identical(as.vector(cache$metrics$barcode), as.vector(cells)),
              file.exists(nuclear_path), file.exists(paste0(nuclear_path, ".tbi")),
              identical(cache$nuclear_signature,
                        core2_qc_file_signature(nuclear_path)))
    core2_qc_log("Resume: independent-fragment QC already complete")
    return(cache$metrics)
  }
  # Never overwrite a previously created nuclear copy.
  if (file.exists(nuclear_path) && !file.exists(paste0(nuclear_path, ".tbi"))) {
    Rsamtools::indexTabix(nuclear_path, format = "bed")
  }
  make_copy <- !file.exists(nuclear_path)
  plain <- if (make_copy) tempfile(
    pattern = ".nuclear_", tmpdir = dirname(nuclear_path), fileext = ".tsv"
  ) else NULL
  tmp_gz <- if (make_copy) paste0(plain, ".gz") else NULL
  on.exit(unlink(c(plain, tmp_gz, if (!is.null(tmp_gz)) paste0(tmp_gz, ".tbi"))),
          add = TRUE)
  core2_qc_log("Count total/nuclear/organelle fragments and exact FRiP in one pass")
  metrics <- core2_qc_scan_fragments(original, peaks, cells, nuclear_plain = plain)
  if (make_copy) {
    core2_qc_log("Compress and index the nuclear-only fragment copy")
    Rsamtools::bgzip(plain, dest = tmp_gz, overwrite = FALSE)
    Rsamtools::indexTabix(tmp_gz, format = "bed")
    stopifnot(setequal(Rsamtools::seqnamesTabix(tmp_gz), as.character(1:5)))
    if (!file.rename(tmp_gz, nuclear_path)) stop("Cannot publish nuclear fragments")
    if (!file.rename(paste0(tmp_gz, ".tbi"), paste0(nuclear_path, ".tbi"))) {
      stop("Nuclear fragments saved, but index rename failed; rerun to index")
    }
  }
  stopifnot(setequal(Rsamtools::seqnamesTabix(nuclear_path), as.character(1:5)))
  # Separate Signac tally independently checks the scanner and any existing copy.
  core2_qc_log("Cross-check nuclear fragment totals with Signac::CountFragments")
  tally <- Signac::CountFragments(nuclear_path, cells = cells, max_lines = NULL)
  index <- match(cells, tally$CB)
  observed <- numeric(length(cells))
  observed[!is.na(index)] <- tally$frequency_count[index[!is.na(index)]]
  stopifnot(all(observed == metrics$nuclear_fragments))
  core2_qc_save(list(signature = signature, metrics = metrics,
                     nuclear_signature = core2_qc_file_signature(nuclear_path)),
                cache_path)
  metrics
}

core2_qc_validate_object <- function(object, counts, peaks) {
  stopifnot(inherits(object, "Seurat"),
            inherits(object[["ATAC"]], "ChromatinAssay"),
            identical(as.vector(colnames(object)), as.vector(colnames(counts))),
            identical(as.vector(rownames(object[["ATAC"]])), as.vector(rownames(counts))),
            identical(core2_qc_ranges_key(GenomicRanges::granges(object[["ATAC"]])),
                      core2_qc_ranges_key(peaks)))
  stored <- SeuratObject::GetAssayData(object = object[["ATAC"]], layer = "counts")
  # fastmatch can attach an internal .match.hash attribute to character labels.
  # Compare the labels above and the integer matrix values here, ignoring that cache.
  stopifnot(identical(dim(stored), dim(counts)), Matrix::nnzero(stored - counts) == 0,
            all(object$nCount_ATAC == Matrix::colSums(counts)),
            all(object$nFeature_ATAC == Matrix::colSums(counts > 0)))
  invisible(TRUE)
}

core2_qc_table <- function(object, sample_id, condition) {
  fields <- c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment",
              "total_fragments", "nuclear_fragments", "fragments_in_common_peaks",
              "FRiP_nuclear", "percent_organelle")
  md <- object[[]]
  stopifnot(all(fields %in% colnames(md)))
  answer <- data.frame(
    sample = sample_id, condition = condition, barcode = rownames(md),
    nucleus_id = paste(sample_id, rownames(md), sep = ":"),
    md[, fields, drop = FALSE], row.names = NULL, check.names = FALSE
  )
  stopifnot(all(answer$total_fragments >= answer$nuclear_fragments),
            all(answer$fragments_in_common_peaks <= answer$nuclear_fragments),
            all(answer$fragments_in_common_peaks <= answer$nCount_ATAC))
  expected_frip <- rep(NA_real_, nrow(answer))
  has_nuclear <- answer$nuclear_fragments > 0
  expected_frip[has_nuclear] <- answer$fragments_in_common_peaks[has_nuclear] /
    answer$nuclear_fragments[has_nuclear]
  expected_organelle <- rep(NA_real_, nrow(answer))
  has_total <- answer$total_fragments > 0
  expected_organelle[has_total] <- 100 * (
    answer$total_fragments[has_total] - answer$nuclear_fragments[has_total]
  ) / answer$total_fragments[has_total]
  stopifnot(isTRUE(all.equal(answer$FRiP_nuclear, expected_frip)),
            isTRUE(all.equal(answer$percent_organelle, expected_organelle)))
  answer$qc_metrics_finite <- has_nuclear & has_total &
    is.finite(answer$nucleosome_signal) & is.finite(answer$TSS.enrichment) &
    is.finite(answer$FRiP_nuclear) & is.finite(answer$percent_organelle)
  answer
}

core2_qc_make_plots <- function(tables, path) {
  grDevices::pdf(path, width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  scatter <- function(x, y, title, xlab, ylab) {
    keep <- is.finite(x) & is.finite(y)
    if (any(keep)) {
      graphics::plot(x[keep], y[keep], pch = 16, cex = 0.3,
                     col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
                     main = title, xlab = xlab, ylab = ylab)
    } else {
      graphics::plot.new()
      graphics::title(main = title, sub = "No finite values")
    }
  }
  for (sample_id in names(tables)) {
    d <- tables[[sample_id]]
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    graphics::hist(log10(d$nCount_ATAC + 1), breaks = 50,
                   main = paste(sample_id, "peak counts"),
                   xlab = "log10(peak overlap counts + 1)", col = "grey80")
    graphics::hist(log10(d$nFeature_ATAC + 1), breaks = 50,
                   main = "Detected peaks", xlab = "log10(detected peaks + 1)",
                   col = "grey80")
    x <- log10(d$nuclear_fragments + 1)
    scatter(x, d$FRiP_nuclear, "Nuclear FRiP", "log10(nuclear fragments + 1)", "FRiP")
    scatter(x, d$TSS.enrichment, "TSS enrichment", "log10(nuclear fragments + 1)",
            "TSS enrichment")
    scatter(x, d$nucleosome_signal, "Nucleosome signal",
            "log10(nuclear fragments + 1)", "Nucleosome signal")
    scatter(x, d$percent_organelle, "Original organelle fraction",
            "log10(nuclear fragments + 1)", "Mt + Pt (%)")
  }
  all_qc <- do.call(rbind, tables)
  all_qc$sample <- factor(all_qc$sample, levels = names(tables))
  graphics::par(mfrow = c(2, 3), mar = c(7, 4, 3, 1))
  metrics <- c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal",
               "TSS.enrichment", "FRiP_nuclear", "percent_organelle")
  for (metric in metrics) {
    y <- all_qc[[metric]]
    if (metric %in% c("nCount_ATAC", "nFeature_ATAC")) y <- log10(y + 1)
    keep <- is.finite(y)
    if (any(keep)) {
      graphics::boxplot(split(y[keep], all_qc$sample[keep]), outline = FALSE,
                        las = 2, col = "lightblue", main = metric, ylab = "", cex.axis = 0.8)
    } else {
      graphics::plot.new()
      graphics::title(main = metric, sub = "No finite values")
    }
  }
}

run_core2_batch_qc <- function(
    task_dir = Sys.getenv("P3_CORE2_WORK", ""),
    reuse_control1 = NULL) {
  revision <- "core2-qc-v1"
  task_dir <- normalizePath(task_dir, mustWork = TRUE)
  .libPaths(c(file.path(task_dir, "work", "R_library"), .libPaths()))
  required <- c("Signac", "SeuratObject", "Matrix", "GenomicRanges", "IRanges",
                "Rsamtools", "rtracklayer")
  absent <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install missing packages first: ", paste(absent, collapse = ", "))
  # Detect API/namespace problems before reading or processing real inputs.
  exports <- list(
    Signac = c("CreateChromatinAssay", "CreateFragmentObject", "Fragments",
               "GetFragmentData", "CountFragments", "NucleosomeSignal", "TSSEnrichment"),
    GenomicRanges = c("GRanges", "GNCList", "resize", "granges", "seqnames",
                      "start", "end", "strand"),
    IRanges = c("IRanges", "overlapsAny"),
    Rsamtools = c("bgzip", "indexTabix", "seqnamesTabix"),
    SeuratObject = c("CreateSeuratObject", "GetAssayData", "AddMetaData")
  )
  for (package in names(exports)) {
    missing <- setdiff(exports[[package]], getNamespaceExports(package))
    if (length(missing)) stop("Unavailable API in ", package, ": ", paste(missing, collapse = ", "))
  }
  samples <- data.frame(
    sample = c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2"),
    tag = c("control1", "control2", "stress1", "stress2"),
    condition = c("control", "control", "stress", "stress"),
    expected_nuclei = c(3431L, 6631L, 7827L, 6585L)
  )
  output_dir <- file.path(task_dir, "outputs", "core2_qc")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  common_path <- file.path(task_dir, "outputs", "core2_common_peaks.rds")
  gtf_path <- file.path(task_dir, "work", "annotation",
                        "Arabidopsis_thaliana.TAIR10.62.gtf.gz")
  stopifnot(file.exists(common_path), file.exists(gtf_path))
  # Check all four upstream inputs before starting the batch.
  for (i in seq_len(nrow(samples))) {
    upstream <- file.path(task_dir, "outputs", paste0(samples$tag[i], "_upstream"))
    original <- file.path(upstream, paste0(samples$tag[i], ".fragments.tsv.gz"))
    stopifnot(file.exists(original), file.exists(paste0(original, ".tbi")),
              file.exists(file.path(upstream, paste0(samples$tag[i], "_common_peak_counts.rds"))),
              file.exists(file.path(upstream, "matrix_barcodes.tsv")))
    contigs <- Rsamtools::seqnamesTabix(original)
    stopifnot(all(as.character(1:5) %in% contigs),
              all(contigs %in% c(as.character(1:5), "Mt", "Pt")))
  }
  common_peaks <- readRDS(common_path)
  stopifnot(inherits(common_peaks, "GRanges"), length(common_peaks) == 27677L,
            all(as.character(GenomicRanges::seqnames(common_peaks)) %in% as.character(1:5)))
  peak_rows <- paste0(as.character(GenomicRanges::seqnames(common_peaks)), "-",
                      GenomicRanges::start(common_peaks), "-", GenomicRanges::end(common_peaks))
  stopifnot(!anyDuplicated(peak_rows))

  core2_qc_log("Load TAIR10 annotation and define one shared protein-coding gene TSS set")
  annotation <- rtracklayer::import(con = gtf_path, format = "gtf")
  genes <- annotation[which(
    annotation$type == "gene" & annotation$gene_biotype == "protein_coding" &
      as.character(GenomicRanges::seqnames(annotation)) %in% as.character(1:5)
  )]
  stopifnot(length(genes) > 0L,
            all(as.character(GenomicRanges::strand(genes)) %in% c("+", "-")))
  tss <- unique(GenomicRanges::resize(genes, width = 1L, fix = "start"))
  # Confirm the TSS set used in the interactive Control 1 session, if available.
  if (exists("tss_positions", envir = .GlobalEnv, inherits = FALSE)) {
    previous_tss <- get("tss_positions", envir = .GlobalEnv)
    stopifnot(identical(core2_qc_ranges_key(previous_tss, strand = TRUE),
                        core2_qc_ranges_key(tss, strand = TRUE)))
  }
  rm(annotation, genes)
  invisible(gc())
  common_signature <- core2_qc_file_signature(common_path, hash = TRUE)
  gtf_signature <- core2_qc_file_signature(gtf_path, hash = TRUE)
  versions <- vapply(required, function(p) as.character(utils::packageVersion(p)), character(1))
  definitions <- list(
    count_unit = "one deduplicated fragment row; PCR support ignored",
    FRiP = "unique nuclear fragments overlapping >=1 common peak / all nuclear fragments",
    organelle = "100 * (all fragments - nuclear fragments) / all fragments",
    nuclear_contigs = as.character(1:5), organelle_contigs = c("Mt", "Pt"),
    TSS = "TAIR10.62 nuclear protein-coding gene starts; strand-aware; deduplicated",
    TSS_count = length(tss), TSS_fast = FALSE, TSS_region_extension = 1000L,
    nucleosome_signal = "Signac::NucleosomeSignal(n=NULL) on nuclear fragments",
    filtering_applied = FALSE, doublet_detection_applied = FALSE
  )
  core2_qc_save(tss, file.path(output_dir, "shared_tss_positions.rds"))
  core2_qc_save(list(revision = revision, common_peaks = common_signature,
                     annotation = gtf_signature, versions = versions, definitions = definitions),
                file.path(output_dir, "qc_definitions.rds"))
  utils::capture.output(utils::sessionInfo(),
                        file = file.path(output_dir, "R_sessionInfo.txt"))

  tables <- list()
  object_paths <- character(nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sample_id <- samples$sample[i]
    tag <- samples$tag[i]
    core2_qc_log("===== ", sample_id, " (", i, "/4) =====")
    upstream <- file.path(task_dir, "outputs", paste0(tag, "_upstream"))
    counts_path <- file.path(upstream, paste0(tag, "_common_peak_counts.rds"))
    barcode_path <- file.path(upstream, "matrix_barcodes.tsv")
    original <- file.path(upstream, paste0(tag, ".fragments.tsv.gz"))
    nuclear_path <- file.path(upstream, paste0(tag, ".nuclear.fragments.tsv.gz"))
    sample_dir <- file.path(output_dir, sample_id)
    dir.create(sample_dir, showWarnings = FALSE)
    object_path <- file.path(sample_dir, "atac_qc_unfiltered.rds")
    object_paths[i] <- object_path
    counts <- readRDS(counts_path)
    cells <- readLines(barcode_path, warn = FALSE)
    stopifnot(inherits(counts, "sparseMatrix"),
              nrow(counts) == length(common_peaks),
              ncol(counts) == samples$expected_nuclei[i],
              identical(as.vector(rownames(counts)), as.vector(peak_rows)),
              identical(as.vector(colnames(counts)), as.vector(cells)), !anyDuplicated(cells),
              all(is.finite(counts@x)), all(counts@x >= 0),
              all(counts@x == floor(counts@x)))
    signature <- list(
      revision = revision, sample = sample_id,
      counts = core2_qc_file_signature(counts_path, hash = TRUE),
      barcodes = core2_qc_file_signature(barcode_path, hash = TRUE),
      original = core2_qc_file_signature(original),
      common_peaks = common_signature, annotation = gtf_signature,
      versions = versions, definitions = definitions
    )
    object <- NULL
    if (file.exists(object_path)) {
      object <- readRDS(object_path)
      provenance <- attr(object, "core2_qc_checkpoint")
      if (is.null(provenance) || !identical(provenance$signature, signature)) {
        stop("Checkpoint inputs/settings changed for ", sample_id,
             ". Existing results were not overwritten: ", object_path)
      }
      stopifnot(file.exists(nuclear_path), file.exists(paste0(nuclear_path, ".tbi")),
                identical(provenance$nuclear_signature,
                          core2_qc_file_signature(nuclear_path)))
      core2_qc_log("Resume saved object: ", provenance$stage)
    } else if (sample_id == "control_rep1") {
      # Import the previously taught, unfiltered Control 1 QC calculation.
      candidate <- reuse_control1
      legacy <- file.path(task_dir, "outputs", "control1_atac_qc_unfiltered.rds")
      if (is.null(candidate) && file.exists(legacy)) candidate <- readRDS(legacy)
      required_fields <- c("nucleosome_signal", "TSS.enrichment", "total_fragments",
                           "nuclear_fragments", "fragments_in_common_peaks",
                           "FRiP_nuclear", "percent_organelle")
      if (!is.null(candidate) && all(required_fields %in% colnames(candidate[[]]))) {
        core2_qc_validate_object(candidate, counts, common_peaks)
        frags <- Signac::Fragments(candidate[["ATAC"]])
        stopifnot(length(frags) == 1L, file.exists(nuclear_path),
                  file.exists(paste0(nuclear_path, ".tbi")),
                  setequal(Rsamtools::seqnamesTabix(nuclear_path), as.character(1:5)),
                  identical(normalizePath(Signac::GetFragmentData(frags[[1]], slot = "path")),
                            normalizePath(nuclear_path)))
        invisible(core2_qc_table(candidate, sample_id, samples$condition[i]))
        object <- candidate
        attr(object, "core2_qc_checkpoint") <- list(
          signature = signature, nuclear_signature = core2_qc_file_signature(nuclear_path),
          stage = "complete", origin = "reused prior interactive Control 1 QC"
        )
        core2_qc_save(object, object_path)
        core2_qc_log("Reused and validated the completed, unfiltered Control 1 object")
      }
      rm(candidate)
    }
    if (is.null(object)) {
      metrics <- core2_qc_fragment_stage(
        original, nuclear_path, common_peaks, cells, signature,
        file.path(sample_dir, "fragment_qc_checkpoint.rds")
      )
      core2_qc_log("Create an ATAC-only ChromatinAssay; retain every supplied nucleus/peak")
      frag <- Signac::CreateFragmentObject(
        path = nuclear_path, cells = stats::setNames(cells, cells)
      )
      assay <- Signac::CreateChromatinAssay(
        counts = counts, ranges = common_peaks, fragments = frag,
        min.cells = 0, min.features = 0
      )
      object <- SeuratObject::CreateSeuratObject(
        counts = assay, assay = "ATAC", project = sample_id,
        min.cells = 0, min.features = 0
      )
      rownames(metrics) <- metrics$barcode
      object <- SeuratObject::AddMetaData(
        object, metadata = metrics[, setdiff(colnames(metrics), "barcode"), drop = FALSE]
      )
      attr(object, "core2_qc_checkpoint") <- list(
        signature = signature, nuclear_signature = core2_qc_file_signature(nuclear_path),
        stage = "fragment metrics complete", origin = "batch calculation"
      )
      core2_qc_save(object, object_path)
      rm(metrics, assay, frag)
    }
    core2_qc_validate_object(object, counts, common_peaks)
    if (!("nucleosome_signal" %in% colnames(object[[]]))) {
      core2_qc_log("Compute nucleosome signal from the entire nuclear fragment file")
      provenance <- attr(object, "core2_qc_checkpoint")
      object <- Signac::NucleosomeSignal(object, assay = "ATAC", n = NULL)
      provenance$stage <- "nucleosome signal complete"
      attr(object, "core2_qc_checkpoint") <- provenance
      core2_qc_save(object, object_path)
    }
    if (!("TSS.enrichment" %in% colnames(object[[]]))) {
      core2_qc_log("Compute TSS enrichment using the same shared TSS positions")
      provenance <- attr(object, "core2_qc_checkpoint")
      object <- Signac::TSSEnrichment(
        object, assay = "ATAC", tss.positions = tss, n = NULL,
        fast = FALSE, region_extension = 1000L
      )
      provenance$stage <- "TSS enrichment complete"
      attr(object, "core2_qc_checkpoint") <- provenance
      core2_qc_save(object, object_path)
    }
    object$sample_id <- sample_id
    object$condition <- samples$condition[i]
    tables[[sample_id]] <- core2_qc_table(object, sample_id, samples$condition[i])
    provenance <- attr(object, "core2_qc_checkpoint")
    provenance$stage <- "complete"
    attr(object, "core2_qc_checkpoint") <- provenance
    core2_qc_save(object, object_path)
    utils::write.csv(tables[[sample_id]], file.path(sample_dir, "nucleus_qc.csv"),
                     row.names = FALSE)
    core2_qc_log("Saved ", ncol(object), " nuclei; no QC filtering applied")
    rm(object, counts)
    invisible(gc())
  }
  all_qc <- do.call(rbind, tables)
  rownames(all_qc) <- NULL
  stopifnot(!anyDuplicated(all_qc$nucleus_id),
            nrow(all_qc) == sum(samples$expected_nuclei))
  utils::write.csv(all_qc, file.path(output_dir, "all_nuclei_qc.csv"), row.names = FALSE)
  metrics <- c("nCount_ATAC", "nFeature_ATAC", "nuclear_fragments",
               "nucleosome_signal", "TSS.enrichment", "FRiP_nuclear", "percent_organelle")
  summary_rows <- list()
  for (sample_id in names(tables)) {
    d <- tables[[sample_id]]
    for (metric in metrics) {
      values <- d[[metric]]
      finite <- values[is.finite(values)]
      q <- if (length(finite)) unname(stats::quantile(finite, c(0, .25, .5, .75, 1)))
      else rep(NA_real_, 5)
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        sample = sample_id, metric = metric, nuclei = nrow(d),
        nonfinite = sum(!is.finite(values)), min = q[1], q25 = q[2], median = q[3],
        mean = if (length(finite)) mean(finite) else NA_real_, q75 = q[4], max = q[5]
      )
    }
  }
  summary_table <- do.call(rbind, summary_rows)
  utils::write.csv(summary_table, file.path(output_dir, "sample_qc_summary.csv"),
                   row.names = FALSE)
  manifest <- data.frame(samples, object_path = object_paths,
                          metrics_finite = vapply(tables, function(x) sum(x$qc_metrics_finite), integer(1)),
                          filtering_applied = FALSE, doublet_detection_applied = FALSE)
  utils::write.csv(manifest, file.path(output_dir, "sample_manifest.csv"), row.names = FALSE)
  core2_qc_make_plots(tables, file.path(output_dir, "four_sample_qc_plots.pdf"))
  completion <- list(status = "QC metrics complete; filtering and doublet detection pending",
                     completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     samples = manifest, definitions = definitions)
  core2_qc_save(completion, file.path(output_dir, "qc_metrics_completion.rds"))
  core2_qc_log("DONE: four samples, ", nrow(all_qc), " unfiltered nuclei. Outputs: ", output_dir)
  print(manifest, row.names = FALSE)
  invisible(list(output_dir = output_dir, manifest = manifest, summary = summary_table))
}

# Repository edition: source defines functions only.
# Explicit call: run_core2_batch_qc(task_dir = Sys.getenv("P3_CORE2_WORK"))
