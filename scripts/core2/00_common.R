# Source defines helpers only. No data are loaded or analyses run here.
core2_ids <- function(x) {
  x <- as.character(x)
  attributes(x) <- NULL
  x
}

core2_read_csv <- function(path, ...) {
  if (!file.exists(path)) path <- paste0(path, ".gz")
  if (!file.exists(path)) stop("Missing CSV or CSV.gz: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, ...)
}

core2_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
}

core2_config <- function() {
  root <- normalizePath(Sys.getenv("P3_ROOT", getwd()), mustWork = TRUE)
  work <- Sys.getenv("P3_CORE2_WORK", "")
  if (!nzchar(work)) stop("Set P3_CORE2_WORK to the local upstream workspace; see docs/core2/running.md.")
  work <- normalizePath(work, mustWork = TRUE)
  library <- Sys.getenv("P3_R_LIB", file.path(work, "work", "R_library"))
  if (dir.exists(library)) .libPaths(c(library, .libPaths()))
  analysis <- Sys.getenv("P3_CORE2_ANALYSIS", file.path(work, "outputs", "core2_rebuild"))
  dir.create(analysis, recursive = TRUE, showWarnings = FALSE)
  list(root = root, work = work, analysis = analysis,
       gtf = file.path(work, "work", "annotation", "Arabidopsis_thaliana.TAIR10.62.gtf.gz"),
       samples = c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2"),
       populations = c("Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage", "Pericycle", "Root_cap"))
}

core2_save <- function(object, config, name) {
  saveRDS(object, file.path(config$analysis, name))
  writeLines(capture.output(utils::sessionInfo()),
             file.path(config$analysis, paste0(name, ".sessionInfo.txt")))
}

# Every counted column must have exactly one metadata row. Never silently intersect.
core2_align_metadata <- function(counts, metadata) {
  cells <- core2_ids(colnames(counts))
  keys <- core2_ids(rownames(metadata))
  stopifnot(length(cells) > 0L, !anyDuplicated(cells), !anyDuplicated(keys),
            setequal(cells, keys))
  metadata <- metadata[match(cells, keys), , drop = FALSE]
  stopifnot(identical(cells, core2_ids(rownames(metadata))))
  metadata
}

# RNA calls belong to this exact nucleus, identified by capture + raw barcode.
core2_selection <- function(metadata, sample, calls) {
  needed <- c("nuclear_fragments", "TSS.enrichment", "nucleosome_signal")
  stopifnot(all(needed %in% names(metadata)),
            all(c("cell", "sample", "qc_pass", "scDblFinder_class", "scDblFinder_score") %in% names(calls)),
            !anyDuplicated(calls$cell))
  cells <- paste0(sample, "_", core2_ids(rownames(metadata)))
  k <- match(cells, calls$cell)
  stopifnot(!anyNA(k), all(calls$sample[k] == sample))
  basic <- is.finite(metadata$nuclear_fragments) &
    is.finite(metadata$TSS.enrichment) & is.finite(metadata$nucleosome_signal) &
    metadata$nuclear_fragments >= 1000 & metadata$TSS.enrichment > 1 &
    metadata$nucleosome_signal < 2
  rna_pass <- !is.na(calls$qc_pass[k]) & calls$qc_pass[k]
  classes <- calls$scDblFinder_class[k]
  paired <- basic & rna_pass & !is.na(classes) & classes == "singlet"
  data.frame(cell = cells, barcode = core2_ids(rownames(metadata)), sample = sample,
             pass_basic_atac = basic, rna_qc_pass = rna_pass,
             rna_doublet_class = classes, rna_doublet_score = calls$scDblFinder_score[k],
             pass_paired_qc = paired, stringsAsFactors = FALSE)
}

core2_pseudobulk_counts <- function(counts, metadata) {
  core2_require("Matrix")
  metadata <- core2_align_metadata(counts, metadata)
  stopifnot(all(c("sample", "condition", "rna_population") %in% names(metadata)),
            !anyNA(metadata[, c("sample", "condition", "rna_population")]))
  group <- paste(metadata$sample, metadata$rna_population, sep = "__")
  levels <- unique(group) # Only observed combinations; no artificial zero samples.
  membership <- Matrix::sparseMatrix(i = seq_along(group), j = match(group, levels),
                                      x = 1, dims = c(length(group), length(levels)),
                                      dimnames = list(core2_ids(colnames(counts)), levels))
  pb <- counts %*% membership
  rows <- lapply(levels, function(id) {
    d <- metadata[group == id, , drop = FALSE]
    stopifnot(length(unique(d$sample)) == 1L, length(unique(d$condition)) == 1L,
              length(unique(d$rna_population)) == 1L)
    data.frame(sample = d$sample[1], condition = d$condition[1],
               rna_population = d$rna_population[1], n_nuclei = nrow(d),
               total_peak_counts = sum(pb[, id]),
               median_nuclear_fragments = stats::median(d$nuclear_fragments),
               median_organelle_percent = stats::median(d$percent_organelle), row.names = id)
  })
  md <- do.call(rbind, rows)
  stopifnot(identical(core2_ids(colnames(pb)), core2_ids(rownames(md))),
            sum(pb) == sum(counts))
  list(counts = pb, metadata = md)
}
