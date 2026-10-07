# Project 3 Core 1: evidence-guided cell population annotation.
# Main method: reference/canonical signature scoring, interpreted manually.
# Dataset-specific strict markers + canonical DotPlot validate the labels.
# This was not an automated reference classifier or an external atlas transfer.
# Final labels/confidence are preserved in data/core1_reference/cluster_annotation.csv.
# Source 00_core1_common.R first; source alone does not execute this function.
run_core1_02 <- function(root = getwd(), output = "results/core1_reproduction",
                         checkpoint = "results/reproducibility/p3_core1_checkpoint_2026-09-30.rds",
                         recompute_evidence = FALSE) {
  core1_require(c("Seurat", "SeuratObject", "dplyr", "tidyr", "AnnotationDbi", "org.At.tair.db", "ggplot2"))
  p <- core1_paths(root, output)
  p3 <- readRDS(file.path(p$root, checkpoint))
  SeuratObject::DefaultAssay(p3) <- "RNA"
  frozen <- file.path(p$root, "results/core1_frozen")
  inputs <- dget(file.path(p$root, "data/core1_reference/annotation_inputs.R"))
  annotation <- utils::read.csv(file.path(p$root, "data/core1_reference/cluster_annotation.csv"),
                                colClasses = "character")

  # Numeric cluster IDs can change on a rerun. Never silently attach historical
  # labels to a different partition. Compare cells AND cluster IDs with evidence.
  historical_md <- utils::read.csv(file.path(frozen, "cell_metadata.csv"),
                                   stringsAsFactors = FALSE)
  stopifnot(ncol(p3) == 24277L, !anyDuplicated(historical_md$cell),
            setequal(colnames(p3), historical_md$cell))
  historical_md <- historical_md[match(colnames(p3), historical_md$cell), , drop = FALSE]
  if (!identical(as.character(p3$seurat_clusters), as.character(historical_md$seurat_clusters))) {
    stop("Cluster identities differ from the reviewed partition. Review/remap labels before annotation.")
  }

  # Seurat v5 initially warned that data layers were not joined and marker
  # testing yielded no usable results. JoinLayers repaired the assay structure;
  # it was not a biological filter or batch correction.
  layers_before <- SeuratObject::Layers(p3[["RNA"]])
  p3 <- SeuratObject::JoinLayers(p3, assay = "RNA")
  writeLines(c("Before:", layers_before, "After:", SeuratObject::Layers(p3[["RNA"]])),
             file.path(p$out, "diagnostics/02_joinlayers.txt"))
  SeuratObject::Idents(p3) <- "seurat_clusters"

  # Strict uniform standard applied to ALL clusters:
  # only.pos=TRUE, min.pct=.25, logfc.threshold=.25, Wilcoxon, return.thresh=.01.
  # Cluster 12 had no returned strict marker. This does not mean no expression.
  # Its relaxed comparison (only.pos=FALSE, min.pct=.10, logfc.threshold=0)
  # served only as a diagnostic and MUST NEVER enter cluster_markers_all.csv.
  if (recompute_evidence) {
    markers <- Seurat::FindAllMarkers(p3, assay = "RNA", only.pos = TRUE,
                                     min.pct = .25, logfc.threshold = .25,
                                     test.use = "wilcox", return.thresh = .01)
    relaxed12 <- Seurat::FindMarkers(p3, assay = "RNA", ident.1 = "12",
                                     only.pos = FALSE, min.pct = .10, logfc.threshold = 0,
                                     test.use = "wilcox")
    core1_write(relaxed12, p, "cluster12_relaxed_DIAGNOSTIC_ONLY.csv",
                sub = "diagnostics", row.names = TRUE)
  } else {
    markers <- utils::read.csv(file.path(p$root, "results/tables/cluster_markers_all.csv"))
    # Relaxed cluster 12 results were not exported in the available checkpoint.
    # Preserve its historical role without inventing a saved diagnostic table.
  }
  core1_write(markers, p, "cluster_markers_all.csv")
  top_markers <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) |>
    dplyr::ungroup()
  top_markers$SYMBOL <- unname(AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
    keys = top_markers$gene, column = "SYMBOL", keytype = "TAIR", multiVals = "first"))
  top_markers$GENENAME <- unname(AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,
    keys = top_markers$gene, column = "GENENAME", keytype = "TAIR", multiVals = "first"))
  core1_write(top_markers, p, "top20_strict_markers_annotated.csv")

  # First small canonical panel -> expanded v2 (8 signatures) -> v3 (12).
  # Compare mean scores and first/second-score margins; avoid assigning weak
  # or negative top scores mechanically. V3 helped inspect clusters 0/1/5/8/12/13.
  # IMPORTANT REVIEW ITEM: the original v3 comments called several added genes
  # candidate/placeholder markers. Their biological provenance is not recorded.
  # The lists are preserved verbatim, not promoted to a validated reference atlas.
  # AddModuleScore used its default seed=1; make that recorded default explicit.
  if (recompute_evidence) {
    p3 <- Seurat::AddModuleScore(p3, features = inputs$marker_panel_v2, name = "Sig", seed = 1)
    panel_v3 <- lapply(inputs$marker_panel_v3, function(x) x[x %in% rownames(p3)])
    stopifnot(all(lengths(panel_v3) > 0L))
    p3 <- Seurat::AddModuleScore(p3, features = panel_v3, name = "SigV3", seed = 1)
  } else {
    # Restore already computed scores by cell identifier; no scoring is rerun.
    score_cols <- base::grep("^Sig(V3)?[0-9]+$", names(historical_md), value = TRUE)
    score_md <- historical_md[, score_cols, drop = FALSE]
    rownames(score_md) <- historical_md$cell
    p3 <- SeuratObject::AddMetaData(p3, metadata = score_md)
  }
  score_map <- data.frame(score_column = paste0("Sig", seq_along(inputs$marker_panel_v2)),
                          cell_type = names(inputs$marker_panel_v2))
  scores <- p3[[]] |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(score_map$score_column), base::mean), .groups = "drop")
  score_long <- scores |>
    tidyr::pivot_longer(cols = dplyr::all_of(score_map$score_column),
                        names_to = "score_column", values_to = "mean_score") |>
    dplyr::left_join(score_map, by = "score_column")
  top2 <- score_long |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::arrange(dplyr::desc(mean_score), .by_group = TRUE) |>
    dplyr::slice_head(n = 2) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()
  margins <- top2 |>
    dplyr::select(seurat_clusters, rank, cell_type, mean_score) |>
    tidyr::pivot_wider(names_from = rank, values_from = c(cell_type, mean_score), names_prefix = "rank") |>
    dplyr::mutate(margin = mean_score_rank1 - mean_score_rank2)
  core1_write(scores, p, "cluster_signature_scores.csv")
  core1_write(margins, p, "cluster_signature_margin.csv")

  # Manual reasoning incorporated score strength/margins and strict markers.
  # Strong initial calls: 2/3/7/15/16/17. Remaining marker evidence distinguished
  # proliferating clusters 9/11 from mature cell types. Conservative labels
  # remain for mixed/weak evidence: *_like, Phloem_related, Unresolved_transition.
  # Cluster 12 retains low confidence; cluster 8 is Epidermal_like, low confidence.
  # Scores do not automatically overwrite these reviewed labels.
  stopifnot(nrow(annotation) == 18L, !anyDuplicated(annotation$cluster))
  label <- annotation$annotation[match(as.character(p3$seurat_clusters), annotation$cluster)]
  confidence <- annotation$confidence[match(as.character(p3$seurat_clusters), annotation$cluster)]
  stopifnot(!anyNA(label), !anyNA(confidence))
  cell_annotation <- data.frame(cell_population = label, annotation_confidence = confidence,
                                 row.names = colnames(p3))
  # Earlier named-vector assignment failed because names were cluster IDs,
  # not barcodes. Attach metadata with explicit cell names instead.
  p3 <- SeuratObject::AddMetaData(p3, metadata = cell_annotation)
  core1_write(annotation, p, "cluster_annotation.csv")
  present <- inputs$validation_markers[inputs$validation_markers %in% rownames(p3)]
  core1_plot(Seurat::DotPlot(p3, features = present, group.by = "cell_population") +
              Seurat::RotatedAxis(), p, "02_canonical_marker_validation.png", width = 12, height = 7)
  core1_plot(Seurat::DimPlot(p3, group.by = "cell_population", label = TRUE),
              p, "02_population_umap.png", width = 12, height = 8)
  saveRDS(p3, file.path(p$out, "checkpoints/02_annotated.rds"))
  core1_session(p, "02")
  invisible(p3)
}
