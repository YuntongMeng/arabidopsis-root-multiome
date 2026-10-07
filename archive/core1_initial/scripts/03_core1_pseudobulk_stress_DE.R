# Project 3 Core 1: biological-sample pseudobulk DE (stress vs control).
# Final: 13 eligible populations, four libraries each, 32,833 common features.
# Frozen per-population DEG counts (up/down):
# Atrichoblast 687 (488/199); Columella 151 (55/96); Cortex 259 (86/173);
# Dividing_meristematic 410 (323/87); Dividing_S_phase 214 (115/99);
# Endodermis 104 (16/88); Endodermis_like 504 (82/422);
# Epidermal_like 535 (275/260); Initial_like 443 (284/159);
# LRC_like 524 (155/369); Pericycle 663 (136/527);
# Trichoblast 739 (251/488); Unresolved_transition 96 (50/46).
# Source 00_core1_common.R first. Default restores frozen final DE results.
# recompute=TRUE explicitly repeats aggregation/DE; it was NOT used for cleanup.
run_core1_03 <- function(root = getwd(), output = "results/core1_reproduction",
                         checkpoint = "results/core1_reproduction/checkpoints/02_annotated.rds",
                         recompute = FALSE) {
  core1_require(c("Matrix", "dplyr"))
  p <- core1_paths(root, output)
  if (!recompute) {
    saved <- readRDS(file.path(p$root, "results/core1_frozen/core1_existing_results.rds"))
    de_all <- dplyr::bind_rows(saved$de_results)
    de_summary <- utils::read.csv(file.path(p$root, "results/core1_frozen/de_summary.csv"))
    core1_write(de_all, p, "pseudobulk_DE_all.csv")
    core1_write(de_summary, p, "pseudobulk_DE_summary.csv")
    saveRDS(saved[c("pb_counts_common", "pb_meta", "eligible_populations", "common_features", "de_results")],
            file.path(p$out, "checkpoints/03_pseudobulk_DE.rds"))
    core1_session(p, "03_restore")
    return(invisible(list(de_all = de_all, de_summary = de_summary)))
  }
  core1_require(c("Seurat", "SeuratObject", "DESeq2"))
  p3 <- readRDS(file.path(p$root, checkpoint))
  samples <- core1_samples(p$root)
  # Nuclei are not biological replicates. Sum raw RNA counts separately for
  # cell_population x sample, leaving two control and two stress libraries.
  sample_counts <- table(factor(p3$cell_population), factor(p3$sample, levels = samples$sample))
  min_counts <- apply(sample_counts, 1, min)
  eligible <- names(min_counts)[min_counts >= 25L]
  # Endodermis meets the boundary at exactly 25 nuclei in control_rep1.
  # Excluded: Atrichoblast_like (min 4), Phloem (8), Phloem_related (8), Xylem (8).
  # Keep their annotations/descriptive evidence, but exclude formal DE.
  expected <- readLines(file.path(p$root, "results/core1_frozen/eligible_populations.txt"))
  stopifnot(length(eligible) == 13L, setequal(eligible, expected))
  core1_write(as.data.frame.matrix(sample_counts), p, "population_sample_counts.csv", row.names = TRUE)
  core1_write(data.frame(cell_population = names(min_counts), min_nuclei = min_counts,
                         eligible = names(min_counts) %in% eligible), p, "pseudobulk_eligibility.csv")
  cells <- colnames(p3)[p3$cell_population %in% eligible]
  primary <- base::subset(p3, cells = cells)
  primary <- SeuratObject::JoinLayers(primary, assay = "RNA")
  pb_counts <- Seurat::AggregateExpression(primary, assays = "RNA",
    group.by = c("cell_population", "sample"), return.seurat = FALSE, slot = "counts")$RNA

  # AggregateExpression converts underscores INSIDE group values to hyphens.
  # Original matching of eligible underscore labels therefore selected no cols.
  # Robust repair: build a key from the original population/sample metadata,
  # applying the same transformation only to the expected aggregation key.
  # Do not parse population names by naively splitting on underscores/hyphens.
  lookup <- primary[[]] |>
    dplyr::select(cell_population, sample, condition, replicate) |>
    dplyr::distinct()
  lookup$sample_id <- paste(gsub("_", "-", lookup$cell_population, fixed = TRUE),
                            gsub("_", "-", lookup$sample, fixed = TRUE), sep = "_")
  stopifnot(!anyDuplicated(lookup$sample_id), setequal(colnames(pb_counts), lookup$sample_id))
  pb_meta <- as.data.frame(lookup[match(colnames(pb_counts), lookup$sample_id), , drop = FALSE])
  rownames(pb_meta) <- pb_meta$sample_id
  stopifnot(identical(colnames(pb_counts), rownames(pb_meta)), nrow(pb_meta) == 52L)

  # Initial Cortex trial used merge's UNION: 126,394 features, 106,442 after
  # row-total >=10 filtering. This was exploratory and is NOT formal DE.
  # Different input feature universes can create artificial structural zeroes.
  # Intersect the four ORIGINAL feature lists before any population's filter.
  common_features <- Reduce(base::intersect, core1_features(p$root, samples))
  stopifnot(length(common_features) == 32833L, all(common_features %in% rownames(pb_counts)))
  pb_counts_common <- pb_counts[common_features, , drop = FALSE]
  stopifnot(all(pb_counts_common@x == round(pb_counts_common@x)))
  core1_write(pb_meta, p, "pseudobulk_metadata.csv")
  writeLines(common_features, file.path(p$out, "diagnostics/common_features.txt"))

  de_results <- list()
  for (pop in eligible) {
    cols <- rownames(pb_meta)[pb_meta$cell_population == pop]
    meta <- pb_meta[cols, , drop = FALSE]
    stopifnot(nrow(meta) == 4L, all(table(meta$condition) == 2L))
    pop_counts <- pb_counts_common[, cols, drop = FALSE]
    # Historical low-count prefilter: sum raw counts over the four libraries >=10.
    pop_counts <- pop_counts[Matrix::rowSums(pop_counts) >= 10, , drop = FALSE]
    stopifnot(all(Matrix::colSums(pop_counts) > 0))
    meta$condition <- factor(meta$condition, levels = c("control", "stress"))
    # Original formal model was ~condition. rep1/rep2 are sample labels; a paired
    # design was not established. Do not silently add a replicate covariate.
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(as.matrix(pop_counts)),
                                        colData = meta, design = ~ condition)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("condition", "stress", "control"), alpha = .05)
    tab <- as.data.frame(res)
    tab$gene <- rownames(tab)
    tab$cell_population <- pop
    de_results[[pop]] <- tab
  }
  de_all <- dplyr::bind_rows(de_results)
  # Formal threshold: padj<.05 AND |log2FC|>=1; positive means stress-up.
  # Preserve rows with NA padj in all-results; they are not significant.
  de_summary <- de_all |>
    dplyr::group_by(cell_population) |>
    dplyr::summarise(genes_tested = dplyr::n(),
      DEGs = sum(!is.na(padj) & padj < .05 & abs(log2FoldChange) >= 1, na.rm = TRUE),
      up = sum(!is.na(padj) & padj < .05 & log2FoldChange >= 1, na.rm = TRUE),
      down = sum(!is.na(padj) & padj < .05 & log2FoldChange <= -1, na.rm = TRUE), .groups = "drop")
  core1_write(de_all, p, "pseudobulk_DE_all.csv")
  core1_write(de_summary, p, "pseudobulk_DE_summary.csv")
  saveRDS(list(pb_counts_common = pb_counts_common, pb_meta = pb_meta,
                common_features = common_features, eligible_populations = eligible,
                de_results = de_results), file.path(p$out, "checkpoints/03_pseudobulk_DE.rds"))
  core1_session(p, "03_recomputed")
  invisible(list(de_all = de_all, de_summary = de_summary))
}
