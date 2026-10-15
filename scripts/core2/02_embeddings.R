# The reviewed working choice is LSI 2:30, not a claim of optimal dimensions.
core2_embeddings <- function(config) {
  core2_require(c("Signac", "Seurat", "SeuratObject"))
  object <- readRDS(file.path(config$analysis, "atac_merged_paired_qc.rds"))
  SeuratObject::DefaultAssay(object) <- "ATAC"
  object <- Signac::RunTFIDF(object, method = 1, scale.factor = 10000)
  object <- Signac::FindTopFeatures(object, min.cutoff = 0)
  set.seed(42)
  object <- Signac::RunSVD(object, n = 50, reduction.name = "lsi", reduction.key = "LSI_")
  embeddings <- SeuratObject::Embeddings(object, "lsi")
  stopifnot(ncol(embeddings) >= 30L)
  depth <- data.frame(dimension = colnames(embeddings),
    peak_count_correlation = apply(embeddings, 2, stats::cor, y = log10(object$nCount_ATAC + 1)),
    nuclear_fragment_correlation = apply(embeddings, 2, stats::cor, y = log10(object$nuclear_fragments + 1)))
  utils::write.csv(depth, file.path(config$analysis, "LSI_depth_correlations.csv"), row.names = FALSE)
  object <- Seurat::RunUMAP(object, reduction = "lsi", dims = 2:30,
    umap.method = "uwot", metric = "cosine", n.neighbors = 30L, min.dist = 0.3,
    seed.use = 42, reduction.name = "umap.atac", reduction.key = "ATACUMAP_")
  object <- Seurat::FindNeighbors(object, reduction = "lsi", dims = 2:30, k.param = 20,
                                   graph.name = c("atac_nn", "atac_snn"))
  object <- Seurat::FindClusters(object, graph.name = "atac_snn", cluster.name = "atac_cluster",
                                   resolution = 0.5, algorithm = 1, random.seed = 42)
  core2_save(object, config, "atac_clustered_r0.5.rds")
  invisible(object)
}
