# Project 3 Core 1: RNA preprocessing and clustering.
# Historical decision trail, reconstructed from .Rhistory and saved outputs.
# Final: 24,277 nuclei; PC1:25; resolution 0.5; 18 clusters.
# Source 00_core1_common.R first. Calling this function explicitly reruns stage 1.
run_core1_01 <- function(root = getwd(), output = "results/core1_reproduction") {
  core1_require(c("Seurat", "SeuratObject", "Matrix", "dplyr", "mclust", "ggplot2"))
  p <- core1_paths(root, output)
  samples <- core1_samples(p$root)

  # 1. Audit every matrix before interpreting QC. Feature IDs come from column 1,
  # not symbols. Barcode prefixes are added only when merging biological samples.
  counts <- lapply(seq_len(nrow(samples)), function(i) {
    stem <- file.path(p$root, "data/rna", samples$sample[i], samples$stem[i])
    x <- Seurat::ReadMtx(mtx = paste0(stem, "_matrix.mtx.gz"),
                        cells = paste0(stem, "_barcodes.tsv.gz"),
                        features = paste0(stem, "_features.tsv.gz"),
                        feature.column = 1, cell.column = 1)
    stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)),
              all(is.finite(x@x)), all(x@x >= 0), all(x@x == round(x@x)))
    x
  })
  names(counts) <- samples$sample
  audit <- data.frame(sample = samples$sample,
                      features = vapply(counts, nrow, numeric(1)),
                      nuclei = vapply(counts, ncol, numeric(1)))
  core1_write(audit, p, "input_matrix_audit.csv")
  # Merge retains a union for exploratory RNA analysis. Stage 3 MUST intersect
  # the original four feature lists before formal DE: absent features are not 0.
  original_features <- core1_features(p$root, samples)
  common_features <- Reduce(base::intersect, original_features)
  writeLines(common_features, file.path(p$out, "diagnostics/common_features.txt"))

  objects <- lapply(seq_along(counts), function(i) {
    obj <- SeuratObject::CreateSeuratObject(counts = counts[[i]], project = samples$sample[i])
    obj$sample <- samples$sample[i]
    obj$condition <- samples$condition[i]
    obj$replicate <- paste0("rep", samples$replicate[i])
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^ATMG")
    obj[["percent.cp"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^ATCG")
    obj
  })
  names(objects) <- samples$sample

  # 2. Inspect distributions and tails before choosing thresholds. The recorded
  # analysis inspected summaries plus 0.5/1/5% and 95/99/99.5% quantiles.
  # For these nuclei, mt/cp percentages were recorded as diagnostics; a generic
  # mitochondrial/chloroplast cutoff was not mechanically imposed.
  qc_distribution <- dplyr::bind_rows(lapply(seq_along(objects), function(i) {
    md <- objects[[i]][[]]
    dplyr::bind_rows(lapply(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.cp"), function(metric) {
      q <- stats::quantile(md[[metric]], probs = c(.005, .01, .05, .5, .95, .99, .995))
      data.frame(sample = samples$sample[i], metric = metric,
                 quantile = names(q), value = unname(q))
    }))
  }))
  core1_write(qc_distribution, p, "qc_distribution_quantiles.csv")
  qc_objects <- lapply(objects, function(obj) {
    md <- obj[[]]
    cells <- rownames(md)[md$nFeature_RNA >= 800 & md$nFeature_RNA <= 15000 & md$nCount_RNA <= 50000]
    base::subset(obj, cells = cells)
  })
  # Seurat dimensions may be returned as double; numeric(1) accepts both
  # integer and double dimensions without changing the counts or QC rules.
  qc_result <- data.frame(sample = samples$sample,
                          before = vapply(objects, ncol, numeric(1)),
                          after = vapply(qc_objects, ncol, numeric(1)))
  qc_result$removed <- qc_result$before - qc_result$after
  qc_result$percent_removed <- 100 * qc_result$removed / qc_result$before
  core1_write(qc_result, p, "qc_retention.csv")
  stopifnot(sum(qc_result$after) == 24277L)
  p3 <- base::merge(x = qc_objects[[1]], y = qc_objects[-1],
                    add.cell.ids = samples$sample, project = "Arabidopsis_root_multiome_RNA")

  # 3. Preserve the original LogNormalize/vst workflow; no regression covariates
  # were added. Stress is the scientific target, so Harmony/batch correction was
  # not applied automatically. Sample composition is checked instead.
  p3 <- Seurat::NormalizeData(p3, normalization.method = "LogNormalize", scale.factor = 10000)
  p3 <- Seurat::FindVariableFeatures(p3, selection.method = "vst", nfeatures = 2000)
  p3 <- Seurat::ScaleData(p3, features = SeuratObject::VariableFeatures(p3))
  p3 <- Seurat::RunPCA(p3, features = SeuratObject::VariableFeatures(p3), npcs = 50, seed.use = 42)
  core1_plot(Seurat::ElbowPlot(p3, ndims = 50), p, "01_pca_elbow.png")

  # 4. ElbowPlot suggested PC20 as a candidate; it did not settle the choice.
  # We tested 10/15/20/25/30 PCs at fixed resolution 0.5 and compared pairwise ARI.
  # Recorded ARI: PC20/25=0.754838; PC25/30=0.854386; PC10/25=0.461993.
  # PC25 was the manually selected balance after this robustness comparison,
  # not an algorithm that selects the maximum ARI or an elbow-only rule.
  pc_values <- c(10L, 15L, 20L, 25L, 30L)
  assignments <- data.frame(cell = colnames(p3))
  final <- NULL
  for (n_pc in pc_values) {
    candidate <- Seurat::FindNeighbors(p3, dims = seq_len(n_pc))
    candidate <- Seurat::FindClusters(candidate, resolution = .5, random.seed = 0)
    assignments[[paste0("PC", n_pc)]] <- as.character(SeuratObject::Idents(candidate))
    if (n_pc == 25L) final <- candidate
  }
  pc_names <- paste0("PC", pc_values)
  ari <- outer(pc_names, pc_names, Vectorize(function(a, b) {
    mclust::adjustedRandIndex(assignments[[a]], assignments[[b]])
  }))
  dimnames(ari) <- list(pc_names, pc_names)
  core1_write(assignments, p, "pc_cluster_assignments.csv")
  core1_write(ari, p, "pc_ari_matrix.csv", row.names = TRUE)
  final$cluster_pc25 <- as.character(SeuratObject::Idents(final))
  stopifnot(length(unique(final$cluster_pc25)) == 18L)
  final <- Seurat::RunUMAP(final, dims = 1:25, seed.use = 42)
  for (group in c("seurat_clusters", "sample", "condition")) {
    core1_plot(Seurat::DimPlot(final, reduction = "umap", group.by = group,
                              label = group == "seurat_clusters"), p, paste0("01_umap_", group, ".png"))
  }

  # 5. Clusters 6/14 were enriched in particular replicates. Inspect their
  # sample fractions AND QC, rather than deleting them or correcting all cells.
  # Current frozen counts: cluster 6 = 223/1096/8/22; 14 = 83/587/4/4.
  composition <- table(cluster = final$seurat_clusters, sample = final$sample)
  core1_write(as.data.frame.matrix(composition), p, "cluster_sample_counts.csv", row.names = TRUE)
  core1_write(as.data.frame.matrix(100 * prop.table(composition, margin = 1)),
              p, "cluster_sample_percent.csv", row.names = TRUE)
  cluster_qc <- final[[]] |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::summarise(n = dplyr::n(), median_nCount = stats::median(nCount_RNA),
                     median_nFeature = stats::median(nFeature_RNA),
                     median_percent_mt = stats::median(percent.mt),
                     median_percent_cp = stats::median(percent.cp), .groups = "drop")
  core1_write(cluster_qc, p, "cluster_qc_summary.csv")
  core1_write(dplyr::filter(cluster_qc, seurat_clusters %in% c("6", "14")),
              p, "cluster_6_14_qc.csv", sub = "diagnostics")
  saveRDS(final, file.path(p$out, "checkpoints/01_clustered.rds"))
  core1_session(p, "01")
  invisible(final)
}
