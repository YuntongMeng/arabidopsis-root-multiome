# Raw peak counts, not TF-IDF, LSI coordinates or gene activity, enter DESeq2.
core2_pseudobulk <- function(config) {
  core2_require(c("Signac", "SeuratObject", "Matrix"))
  object <- readRDS(file.path(config$analysis, "atac_annotated_checkpoint.rds"))
  counts <- SeuratObject::LayerData(object, assay = "ATAC", layer = "counts")
  input <- core2_pseudobulk_counts(counts, object[[]])
  core2_save(input, config, "atac_pseudobulk_unfiltered.rds")
  utils::write.csv(input$metadata, file.path(config$analysis, "atac_pseudobulk_metadata.csv"))
  invisible(input)
}
