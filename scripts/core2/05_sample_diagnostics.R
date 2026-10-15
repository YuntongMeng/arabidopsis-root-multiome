# Core 2: sample diagnostics within five RNA-annotated populations.
# source() only defines this function. The function call runs the analysis.
# Uses raw peak pseudobulk counts; preserves the four biological samples.
# No cells are filtered, no samples are removed, and no differential peak
# results are selected or interpreted here. Models are saved for the next step.

core2_sample_diagnostics <- function(
    analysis_dir,
    populations = c(
      "Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage",
      "Pericycle", "Root_cap"
    ),
    min_nuclei = 50L,
    min_peak_count = 10L,
    min_peak_samples = 2L,
    show_plots = TRUE) {

  packages <- c("DESeq2", "SummarizedExperiment", "ggplot2", "patchwork", "ggrepel")
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing R package: ", package)
    }
  }

  plain_ids <- function(x) {
    x <- as.character(x)
    attributes(x) <- NULL
    x
  }

  input_path <- file.path(analysis_dir, "atac_pseudobulk_unfiltered.rds")
  stopifnot(file.exists(input_path))
  input <- readRDS(input_path)
  counts <- input$counts
  metadata <- input$metadata
  stopifnot(
    all(c("sample", "condition", "rna_population", "n_nuclei") %in%
          colnames(metadata)),
    !anyDuplicated(plain_ids(colnames(counts))),
    !anyDuplicated(plain_ids(rownames(metadata))),
    setequal(plain_ids(colnames(counts)), plain_ids(rownames(metadata)))
  )

  sample_order <- c(
    "control_rep1", "control_rep2", "stress_rep1", "stress_rep2"
  )
  expected_condition <- c("control", "control", "stress", "stress")
  output_dir <- file.path(analysis_dir, "sample_diagnostics")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_rows <- list()
  pca_plots <- list()
  correlation_plots <- list()

  for (population in populations) {
    message("Checking population: ", population)
    md <- metadata[
      which(metadata$rna_population == population), , drop = FALSE
    ]
    stopifnot(
      nrow(md) == 4L,
      !anyDuplicated(as.character(md$sample)),
      setequal(as.character(md$sample), sample_order),
      all(is.finite(md$n_nuclei)),
      all(md$n_nuclei >= min_nuclei)
    )
    md <- md[base::match(sample_order, md$sample), , drop = FALSE]
    stopifnot(identical(as.character(md$condition), expected_condition))
    md$pseudobulk_id <- plain_ids(rownames(md))
    m <- as.matrix(counts[, md$pseudobulk_id, drop = FALSE])
    colnames(m) <- sample_order
    rownames(md) <- sample_order
    md$condition <- factor(md$condition, levels = c("control", "stress"))
    stopifnot(
      all(is.finite(m)), all(m >= 0), all(m == round(m)),
      max(m) <= .Machine$integer.max
    )
    storage.mode(m) <- "integer"
    keep <- rowSums(m >= min_peak_count) >= min_peak_samples
    stopifnot(sum(keep) >= 3L)
    m <- m[keep, , drop = FALSE]

    population_dir <- file.path(output_dir, population)
    dir.create(population_dir, recursive = TRUE, showWarnings = FALSE)
    model_path <- file.path(population_dir, "deseq2_model.rds")
    reuse_path <- model_path
    if (!file.exists(reuse_path) && population == "Endodermis_lineage") {
      reuse_path <- file.path(analysis_dir, "Endodermis_deseq2_model.rds")
    }

    if (file.exists(reuse_path)) {
      dds <- readRDS(reuse_path)
      previous_counts <- DESeq2::counts(dds)
      if (!identical(dim(previous_counts), dim(m)) ||
          !identical(plain_ids(rownames(previous_counts)), plain_ids(rownames(m))) ||
          !identical(plain_ids(colnames(previous_counts)), sample_order) ||
          !all(previous_counts == m) ||
          !identical(as.character(dds$condition), expected_condition) ||
          !identical(all.vars(DESeq2::design(dds)), "condition") ||
          is.null(DESeq2::sizeFactors(dds))) {
        stop("Saved model does not match current inputs: ", reuse_path)
      }
      message("Reusing matching fitted model.")
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = m, colData = md, design = ~ condition
      )
      dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
    }
    saveRDS(dds, model_path)

    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    vst_matrix <- SummarizedExperiment::assay(vsd)
    stopifnot(all(is.finite(vst_matrix)))
    correlation <- stats::cor(vst_matrix, method = "pearson")
    stopifnot(all(is.finite(correlation)))
    pca <- DESeq2::plotPCA(
      vsd, intgroup = c("condition", "sample"),
      ntop = min(500L, nrow(vsd)), returnData = TRUE
    )
    percent_var <- 100 * attr(pca, "percentVar")
    diagnostics <- data.frame(
      sample = sample_order,
      n_nuclei = md$n_nuclei,
      counts_after_peak_filter = colSums(m),
      size_factor = as.numeric(DESeq2::sizeFactors(dds)),
      zero_peak_percent = 100 * colMeans(m == 0L),
      normalized_total = colSums(DESeq2::counts(dds, normalized = TRUE)),
      row.names = NULL
    )
    utils::write.csv(diagnostics, file.path(population_dir, "sample_metrics.csv"),
                     row.names = FALSE)
    utils::write.csv(correlation, file.path(population_dir, "sample_correlations.csv"))
    utils::write.csv(pca, file.path(population_dir, "sample_pca_coordinates.csv"),
                     row.names = FALSE)

    title <- gsub("_", " ", population)
    limit <- max(abs(c(pca$PC1, pca$PC2))) * 1.3
    if (limit == 0) limit <- 1
    p_pca <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color = condition)) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(ggplot2::aes(label = sample), seed = 42,
        size = 3, max.overlaps = Inf, show.legend = FALSE) +
      ggplot2::coord_fixed(ratio = 1, xlim = c(-limit, limit), ylim = c(-limit, limit), expand = FALSE) +
      ggplot2::labs(title = title,
        x = sprintf("PC1: %.1f%%", percent_var[1]), y = sprintf("PC2: %.1f%%", percent_var[2])) +
      ggplot2::theme_classic(base_size = 11)

    cor_df <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
    colnames(cor_df) <- c("sample_x", "sample_y", "correlation")
    cor_df$sample_x <- factor(cor_df$sample_x, levels = sample_order)
    cor_df$sample_y <- factor(cor_df$sample_y, levels = rev(sample_order))
    p_cor <- ggplot2::ggplot(
      cor_df, ggplot2::aes(x = sample_x, y = sample_y, fill = correlation)
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.3f", correlation),
                     color = correlation >= 0.75), size = 3.5
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "black", "TRUE" = "white"), guide = "none"
      ) +
      ggplot2::scale_fill_gradient(
        low = "#F7FBFF", high = "#08519C", limits = c(0, 1), name = "Pearson r"
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
      )

    ggplot2::ggsave(file.path(population_dir, "sample_PCA.png"),
                    p_pca, width = 9, height = 6, dpi = 300)
    ggplot2::ggsave(file.path(population_dir, "sample_correlations.png"),
                    p_cor, width = 8, height = 7, dpi = 300)
    pca_plots[[population]] <- p_pca
    correlation_plots[[population]] <- p_cor
    summary_rows[[population]] <- data.frame(
      population = population,
      n_peaks = nrow(m),
      control_r = correlation["control_rep1", "control_rep2"],
      stress_r = correlation["stress_rep1", "stress_rep2"],
      PC1_percent = percent_var[1],
      PC2_percent = percent_var[2],
      row.names = NULL
    )
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  utils::write.csv(summary, file.path(output_dir, "diagnostic_summary.csv"),
                   row.names = FALSE)
  combined_pca <- patchwork::wrap_plots(pca_plots, ncol = 2) +
    patchwork::plot_annotation(caption = "Independent PCA per population; equal x/y unit scale within each panel; panel ranges differ.")
  combined_cor <- patchwork::wrap_plots(correlation_plots, ncol = 2)
  ggplot2::ggsave(file.path(output_dir, "all_populations_PCA.png"),
                  combined_pca, width = 16, height = 14, dpi = 250,
                  limitsize = FALSE)
  ggplot2::ggsave(file.path(output_dir, "all_populations_correlations.png"),
                  combined_cor, width = 16, height = 16, dpi = 250,
                  limitsize = FALSE)
  saveRDS(list(
    populations = populations, min_nuclei = min_nuclei,
    min_peak_count = min_peak_count, min_peak_samples = min_peak_samples,
    design = "~ condition", size_factor_method = "poscounts",
    PCA_ntop = 500L, VST_blind = TRUE,
    role = "Sample diagnostics; no sample exclusion; no selected differential peaks"
  ), file.path(output_dir, "parameters.rds"))
  writeLines(capture.output(utils::sessionInfo()),
              file.path(output_dir, "sessionInfo.txt"))
  if (show_plots) print(combined_pca)
  display <- summary
  display$control_r <- round(display$control_r, 3)
  display$stress_r <- round(display$stress_r, 3)
  display$PC1_percent <- round(display$PC1_percent, 1)
  display$PC2_percent <- round(display$PC2_percent, 1)
  print(display, row.names = FALSE)
  invisible(summary)
}
