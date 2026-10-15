# Reconstructed from the reviewed interactive procedure. Source defines only.
core2_pair_and_merge <- function(config) {
  core2_require(c("Signac", "SeuratObject"))
  calls <- core2_read_csv(file.path(config$root, "results", "core1", "tables", "qc_doublet_calls_all.csv"))
  objects <- list()
  summary <- list()
  for (sample in config$samples) {
    object <- readRDS(file.path(config$work, "outputs", "core2_qc", sample, "atac_qc_unfiltered.rds"))
    selection <- core2_selection(object[[]], sample, calls)
    object$sample <- sample
    object$condition <- if (startsWith(sample, "control")) "control" else "stress"
    object$barcode <- selection$barcode
    object$rna_qc_pass <- selection$rna_qc_pass
    object$rna_doublet_class <- selection$rna_doublet_class
    object$rna_doublet_score <- selection$rna_doublet_score
    sample_dir <- file.path(config$analysis, "selection", sample)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(selection, file.path(sample_dir, "selection.csv"), row.names = FALSE)
    objects[[sample]] <- base::subset(object, cells = colnames(object)[selection$pass_paired_qc])
    summary[[sample]] <- data.frame(sample = sample, initial = ncol(object),
      basic_atac = sum(selection$pass_basic_atac), paired_singlet = ncol(objects[[sample]]))
  }
  reference <- core2_ids(rownames(objects[[1]][["ATAC"]]))
  stopifnot(all(vapply(objects, function(o) identical(reference, core2_ids(rownames(o[["ATAC"]]))), logical(1))))
  merged <- base::merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects),
                         merge.data = FALSE, merge.dr = FALSE)
  stopifnot(!anyDuplicated(core2_ids(colnames(merged))),
    all(core2_ids(colnames(merged)) == paste0(merged$sample, "_", merged$barcode)))
  utils::write.csv(do.call(rbind, summary), file.path(config$analysis, "paired_selection_summary.csv"), row.names = FALSE)
  core2_save(merged, config, "atac_merged_paired_qc.rds")
  invisible(merged)
}
