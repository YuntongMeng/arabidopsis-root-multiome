# Exact-nucleus RNA matching; no anchors, WNN, motif analysis or target links.
core2_annotation_and_activity <- function(config) {
  core2_require(c("Signac", "Seurat", "SeuratObject", "GenomicRanges", "rtracklayer"))
  object <- readRDS(file.path(config$analysis, "atac_clustered_r0.5.rds"))
  rna <- core2_read_csv(file.path(config$root, "results", "core1", "tables", "annotated_cell_metadata.csv"))
  review <- core2_read_csv(file.path(config$root, "data", "core1_corrected_reference", "cluster_annotation_review.csv"))
  review <- review[review$branch == "pca", , drop = FALSE]
  stopifnot(!anyDuplicated(rna$cell), !anyDuplicated(review$cluster))
  k <- match(core2_ids(colnames(object)), rna$cell)
  stopifnot(!anyNA(k), all(rna$sample[k] == object$sample), all(rna$condition[k] == object$condition))
  j <- match(rna$pca_cluster[k], review$cluster)
  stopifnot(!anyNA(j), all(rna$pca_population[k] == review$cell_population[j]))
  object$rna_population <- rna$pca_population[k]
  object$rna_annotation_confidence <- review$confidence[j]

  composition <- 100 * prop.table(table(object$atac_cluster, object$rna_population), margin = 1)
  utils::write.csv(composition, file.path(config$analysis, "atac_cluster_rna_composition_percent.csv"))
  # A reviewed description is bound to cell-level cluster assignments, not ID numbers alone.
  old <- core2_read_csv(file.path(config$root, "data", "core2_reference", "reviewed_cell_clusters.csv"))
  old_k <- match(core2_ids(colnames(object)), old$cell)
  reviewed <- !anyNA(old_k) && nrow(old) == ncol(object) &&
    identical(as.character(object$atac_cluster), as.character(old$atac_cluster[old_k]))
  descriptions <- core2_read_csv(file.path(config$root, "data", "core2_reference", "atac_cluster_annotation_review.csv"))
  if (reviewed) {
    z <- match(as.character(object$atac_cluster), as.character(descriptions$atac_cluster))
    stopifnot(!anyNA(z))
    object$atac_cluster_description <- descriptions$description[z]
  } else {
    object$atac_cluster_description <- "Unreviewed ATAC cluster"
    warning("Cluster assignments differ from the reviewed run; old descriptions were not reused.")
  }
  writeLines(paste("Reviewed cell-level cluster assignment match:", reviewed),
             file.path(config$analysis, "annotation_review_status.txt"))

  annotation <- rtracklayer::import(config$gtf, format = "gtf")
  annotation <- annotation[which(annotation$type %in% c("transcript", "exon") &
    annotation$gene_biotype == "protein_coding" &
    as.character(GenomicRanges::seqnames(annotation)) %in% as.character(1:5))]
  annotation$tx_id <- annotation$transcript_id
  missing <- is.na(annotation$gene_name) | annotation$gene_name == ""
  annotation$gene_name[missing] <- annotation$gene_id[missing]
  Signac::Annotation(object[["ATAC"]]) <- annotation
  activity <- Signac::GeneActivity(object, assay = "ATAC", features = NULL,
    biotypes = "protein_coding", extend.upstream = 500, extend.downstream = 0, gene.id = TRUE)
  core2_save(activity, config, "atac_gene_activity_counts.rds")
  object[["GENE_ACTIVITY"]] <- SeuratObject::CreateAssayObject(counts = activity)
  object <- Seurat::NormalizeData(object, assay = "GENE_ACTIVITY", normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObject::DefaultAssay(object) <- "ATAC"
  core2_save(object, config, "atac_annotated_checkpoint.rds")
  invisible(object)
}

core2_marker_and_coverage_plots <- function(config) {
  core2_require(c("Signac", "Seurat", "SeuratObject", "ggplot2", "patchwork"))
  object <- readRDS(file.path(config$analysis, "atac_annotated_checkpoint.rds"))
  out <- file.path(config$analysis, "marker_checks")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  markers <- c("AT1G79840", "AT5G14750", "AT1G62500", "AT1G62510", "AT3G54220", "AT5G57620",
               "AT2G31310", "AT4G30450", "AT1G79580", "AT1G33280", "AT5G49270", "AT1G66470")
  labels <- c("GL2", "WER", "CO2", "AT1G62510", "SCR", "MYB36", "LBD14", "XPP", "SMB", "BRN1", "COBL9", "RHD6")
  stopifnot(all(markers %in% rownames(object[["GENE_ACTIVITY"]])))
  dot <- Seurat::DotPlot(object, assay = "GENE_ACTIVITY", features = markers,
    group.by = "atac_cluster", scale = FALSE) +
    ggplot2::scale_x_discrete(labels = stats::setNames(labels, markers)) +
    ggplot2::labs(title = "Marker gene-region ATAC signal", x = "Gene", y = "ATAC cluster",
                   color = "Mean normalized\nATAC signal", size = "Nuclei with\nATAC signal (%)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(out, "marker_gene_regions.png"), dot, width = 11, height = 6, dpi = 180)
  annotation <- Signac::Annotation(object[["ATAC"]])
  scr <- annotation[which(annotation$type == "transcript" & annotation$gene_id == "AT3G54220")]
  stopifnot(length(scr) > 0L)
  region <- base::range(scr, ignore.strand = TRUE)
  p <- Signac::CoveragePlot(object, assay = "ATAC", region = region,
    group.by = "atac_cluster", extend.upstream = 2000, extend.downstream = 2000,
    annotation = TRUE, peaks = TRUE, links = FALSE) +
    patchwork::plot_annotation(title = "SCR region: ATAC clusters")
  ggplot2::ggsave(file.path(out, "SCR_cluster_coverage.pdf"), p, width = 12, height = 9)
  invisible(out)
}
