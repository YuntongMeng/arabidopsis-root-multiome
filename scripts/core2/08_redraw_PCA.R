# Correct PCA display using saved coordinates only: no count/model/PCA refit.
# source() defines the function; call core2_redraw_PCA(analysis_dir) to draw.

core2_redraw_PCA <- function(analysis_dir, show_plot = TRUE, output_dir = NULL) {
  for (package in c("ggplot2", "patchwork", "ggrepel")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  }
  diagnostics_dir <- file.path(analysis_dir, "sample_diagnostics")
  summary_path <- file.path(diagnostics_dir, "diagnostic_summary.csv")
  stopifnot(file.exists(summary_path))
  summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(summary$population))
  out_dir <- if (is.null(output_dir)) file.path(analysis_dir, "final_figures", "PCA_equal_units") else output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()
  samples <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
  for (i in seq_len(nrow(summary))) {
    population <- summary$population[i]
    pca_path <- file.path(diagnostics_dir, population, "sample_pca_coordinates.csv")
    stopifnot(file.exists(pca_path))
    d <- utils::read.csv(pca_path, stringsAsFactors = FALSE)
    stopifnot(
      all(c("PC1", "PC2", "condition", "sample") %in% names(d)),
      nrow(d) == 4L, setequal(d$sample, samples),
      all(is.finite(d$PC1)), all(is.finite(d$PC2))
    )
    # The same numerical range on x/y makes each panel square AND gives equal
    # physical length to one PC1 unit and one PC2 unit. Ranges differ by panel.
    limit <- max(abs(c(d$PC1, d$PC2))) * 1.3
    if (limit == 0) limit <- 1
    p <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = condition)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
      ggplot2::geom_point(size = 3.5) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample), seed = 42, size = 3,
        max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(control = "#F8766D", stress = "#00BFC4"),
        breaks = c("control", "stress"), drop = FALSE
      ) +
      ggplot2::coord_fixed(
        ratio = 1, xlim = c(-limit, limit), ylim = c(-limit, limit),
        expand = FALSE
      ) +
      ggplot2::labs(
        title = gsub("_", " ", population),
        x = sprintf("PC1: %.1f%%", summary$PC1_percent[i]),
        y = sprintf("PC2: %.1f%%", summary$PC2_percent[i])
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    plots[[population]] <- p
    ggplot2::ggsave(file.path(out_dir, paste0(population, ".png")),
                    p, width = 8, height = 7, dpi = 300)
    ggplot2::ggsave(file.path(out_dir, paste0(population, ".pdf")),
                    p, width = 8, height = 7)
  }
  combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "ATAC pseudobulk: sample relationships within RNA populations",
      caption = paste(
        "Each panel is an independent PCA; coordinates are not a shared space.\n",
        "Within each panel, x and y use equal unit scales. Axis ranges differ by panel.\n",
        "Points are biological samples. PCA uses the 500 most variable retained peaks."
      )
    )
  combined <- combined & ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(out_dir, "all_populations_PCA_equal_units.png"),
                  combined, width = 15, height = 17, dpi = 250, limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, "all_populations_PCA_equal_units.pdf"),
                  combined, width = 15, height = 17, limitsize = FALSE)
  if (show_plot) print(combined)
  invisible(list(output_dir = out_dir, combined = combined, panels = plots))
}
