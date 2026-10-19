# Redraw saved exploratory effect-size panels only; never fits a model.
core2_plot_response_context <- function(response_dir,
    populations = c("Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage", "Pericycle", "Root_cap")) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  gene_summary <- utils::read.csv(file.path(response_dir, "RNA_ATAC_gene_summary.csv"), stringsAsFactors = FALSE)
  long <- utils::read.csv(file.path(response_dir, "RNA_ATAC_peak_context.csv"), stringsAsFactors = FALSE)
  for (population in populations) {
    g <- gene_summary[gene_summary$population == population, , drop = FALSE]
    d <- long[long$population == population & long$ATAC_in_model &
                long$peak %in% g$representative_peak, , drop = FALSE]
    d <- d[paste(d$gene_id, d$peak) %in% paste(g$gene_id, g$representative_peak), , drop = FALSE]
    if (!nrow(d)) next
    d$label <- paste0(d$gene_name, "\n", d$gene_id)
    dd <- rbind(data.frame(label = d$label, modality = "RNA", log2FC = d$RNA_log2FC,
                           SE = d$RNA_lfcSE, padj = d$RNA_padj),
                data.frame(label = d$label, modality = "Nearby ATAC peak", log2FC = d$ATAC_log2FoldChange,
                           SE = d$ATAC_lfcSE, padj = d$ATAC_padj))
    dd$modality <- factor(dd$modality, c("RNA", "Nearby ATAC peak"))
    p <- ggplot2::ggplot(dd, ggplot2::aes(x = modality, y = log2FC, color = modality)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = log2FC - 1.96 * SE, ymax = log2FC + 1.96 * SE), width = .12) +
      ggplot2::geom_point(size = 3) + ggplot2::facet_wrap(~ label, nrow = 1) +
      ggplot2::scale_color_manual(values = c("RNA" = "#008B8B", "Nearby ATAC peak" = "#B97725")) +
      ggplot2::labs(title = paste(gsub("_", " ", population), "- exploratory RNA / ATAC context"),
        subtitle = "Top RNA responses; representative ATAC peak chosen by position, not significance",
        x = NULL, y = "Estimated stress / control log2 fold change",
        caption = paste(sprintf("Displayed %d of %d RNA-selected genes; omitted genes have no eligible primary ATAC peak (see table).", nrow(d), nrow(g)),
                        "Bars: model estimate +/- 1.96 SE, not simultaneous intervals. RNA and ATAC retained nuclei differ.",
                        "Nearby peak is not a confirmed target link; RNA significance does not imply ATAC significance.", sep = "\n")) +
      ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
        plot.caption = ggplot2::element_text(hjust = 0, size = 8))
    ggplot2::ggsave(file.path(response_dir, paste0(population, "_RNA_ATAC_context.png")),
                     p, width = 11, height = 6, dpi = 220)
  }
  invisible(response_dir)
}

# Supplemental Core 2 review. Sourcing defines a function; it runs no analysis.
# Original models/results and all Core 1 files are read-only inputs.
core2_supplemental_review <- function(
    task_dir, project_dir,
    analysis_dir = file.path(task_dir, "outputs", "core2_analysis"),
    output_dir = file.path(task_dir, "outputs", "core2_supplement_2026-10-15"),
    exploratory_count = 5L, top_rna_genes = 3L) {
  packages <- c("DESeq2", "SummarizedExperiment", "S4Vectors", "rtracklayer",
                "GenomicRanges", "IRanges", "ggplot2")
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  }
  stopifnot(exploratory_count > 0L, exploratory_count < 10L,
            exploratory_count == as.integer(exploratory_count),
            top_rna_genes >= 1L, top_rna_genes == as.integer(top_rna_genes))
  task_dir <- normalizePath(task_dir, mustWork = TRUE)
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  if (file.exists(output_dir)) stop("Refusing to overwrite existing supplement: ", output_dir)
  populations <- c("Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage",
                   "Pericycle", "Root_cap")
  samples <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
  conditions <- c("control", "control", "stress", "stress")
  ids <- function(x) unname(as.character(x))
  write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "NA")
  bind <- function(x) { y <- do.call(rbind, x); rownames(y) <- NULL; y }
  source_csv <- function(path) {
    if (!file.exists(path)) path <- paste0(path, ".gz")
    stopifnot(file.exists(path))
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  model_paths <- file.path(analysis_dir, "sample_diagnostics", populations, "deseq2_model.rds")
  result_paths <- file.path(analysis_dir, "differential_accessibility", populations,
                            "all_peak_results.csv")
  rna_paths <- file.path(project_dir, "results", "core1", "tables",
                        paste0("DE_pca__", populations, ".csv.gz"))
  pb_path <- file.path(analysis_dir, "atac_pseudobulk_unfiltered.rds")
  common_path <- file.path(task_dir, "outputs", "core2_common_peaks.rds")
  gtf_path <- file.path(task_dir, "work", "annotation", "Arabidopsis_thaliana.TAIR10.62.gtf.gz")
  inputs <- c(pb_path, model_paths, result_paths, rna_paths, common_path, gtf_path)
  stopifnot(all(file.exists(inputs)))
  hashes_before <- tools::md5sum(inputs)
  pb <- readRDS(pb_path)
  stopifnot(setequal(ids(colnames(pb$counts)), ids(rownames(pb$metadata))))
  dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir)
  diagnostic_dir <- file.path(output_dir, "model_diagnostics")
  sensitivity_dir <- file.path(output_dir, paste0("sensitivity_count", exploratory_count))
  response_dir <- file.path(output_dir, "RNA_ATAC_response")
  for (path in c(diagnostic_dir, sensitivity_dir, response_dir)) dir.create(path)

  extract_results <- function(dds, population) {
    r <- as.data.frame(DESeq2::results(dds, contrast = c("condition", "stress", "control"),
      alpha = .05, pAdjustMethod = "BH", independentFiltering = TRUE))
    norm <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
    raw <- as.data.frame(DESeq2::counts(dds))
    names(norm) <- paste0("normalized_", samples)
    names(raw) <- paste0("raw_", samples)
    x <- data.frame(peak = ids(rownames(r)), r, norm, raw, population = population,
                    row.names = NULL, check.names = FALSE)
    x$FDR_significant <- !is.na(x$padj) & x$padj < .05
    x$priority_change <- x$FDR_significant & is.finite(x$log2FoldChange) &
      abs(x$log2FoldChange) >= 1
    x$replicate_pattern <- "Groups_overlap"
    nn <- as.matrix(norm)
    x$replicate_pattern[pmin(nn[, 3], nn[, 4]) > pmax(nn[, 1], nn[, 2])] <-
      "Both_stress_above_both_control"
    x$replicate_pattern[pmax(nn[, 3], nn[, 4]) < pmin(nn[, 1], nn[, 2])] <-
      "Both_stress_below_both_control"
    x
  }
  models <- lapply(model_paths, readRDS); names(models) <- populations
  primary <- lapply(result_paths, source_csv); names(primary) <- populations
  raw_matrices <- list()
  validation <- list()
  for (population in populations) {
    dds <- models[[population]]
    md <- pb$metadata[which(pb$metadata$rna_population == population), , drop = FALSE]
    stopifnot(nrow(md) == 4L, !anyDuplicated(md$sample), setequal(md$sample, samples))
    md <- md[match(samples, md$sample), , drop = FALSE]
    stopifnot(all(md$n_nuclei >= 50L), identical(as.character(md$condition), conditions))
    m <- as.matrix(pb$counts[, rownames(md), drop = FALSE])
    colnames(m) <- samples
    stopifnot(all(is.finite(m)), all(m >= 0), all(m == round(m)))
    storage.mode(m) <- "integer"
    raw_matrices[[population]] <- m
    m10 <- m[rowSums(m >= 10L) >= 2L, , drop = FALSE]
    stopifnot(identical(ids(colnames(dds)), samples),
      identical(as.character(dds$condition), conditions),
      identical(all.vars(DESeq2::design(dds)), "condition"),
      identical(ids(rownames(dds)), ids(rownames(m10))),
      identical(unname(DESeq2::counts(dds)), unname(m10)))
    old <- primary[[population]]
    current <- extract_results(dds, population)
    current <- current[match(old$peak, current$peak), , drop = FALSE]
    stopifnot(!anyNA(current$peak), !anyDuplicated(old$peak))
    fields <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                paste0("normalized_", samples), paste0("raw_", samples))
    equal <- vapply(fields, function(field) {
      isTRUE(all.equal(old[[field]], current[[field]], tolerance = 1e-7,
                       check.attributes = FALSE))
    }, logical(1))
    stopifnot(all(equal))
    validation[[population]] <- data.frame(population = population,
      original_counts_match_pseudobulk = TRUE, original_statistics_reproduced = TRUE)
  }

  # Ten full-size figures: five p distributions and five dispersion figures.
  # Each dispersion figure has a separate count/proportion panel; overplotting
  # in the scatter is never used to estimate the numerical floor proportion.
  min_disp <- 1e-8 # default used by the recorded DESeq2 1.52.0 fits
  diag_rows <- list(); strata_rows <- list()
  all_means <- unlist(lapply(models, function(x) S4Vectors::mcols(x)$baseMean))
  all_disp <- unlist(lapply(models, function(x) {
    mm <- S4Vectors::mcols(x); c(mm$dispGeneEst, mm$dispFit, mm$dispersion)
  }))
  xlim <- range(all_means[is.finite(all_means) & all_means > 0])
  ylim <- range(c(min_disp, all_disp[is.finite(all_disp) & all_disp > 0]))
  plot_device <- function(stem, draw) {
    grDevices::pdf(paste0(stem, ".pdf"), width = 10.5, height = 8, useDingbats = FALSE)
    tryCatch(draw(), finally = grDevices::dev.off())
    grDevices::png(paste0(stem, ".png"), width = 2100, height = 1600, res = 200)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  for (population in populations) {
    message("Supplemental diagnostic figures: ", population)
    dds <- models[[population]]; r <- primary[[population]]
    p <- r$pvalue[is.finite(r$pvalue)]
    mc <- S4Vectors::mcols(dds)
    valid <- is.finite(mc$dispGeneEst)
    floor <- valid & mc$dispGeneEst <= min_disp * 1.01
    n_floor <- sum(floor); n_est <- sum(valid)
    fit <- attr(DESeq2::dispersionFunction(dds), "fitType")
    diag_rows[[population]] <- data.frame(population = population,
      model_peaks = nrow(dds), finite_pvalues = length(p), raw_p_below_001 = sum(p < .01),
      raw_p_below_005 = sum(p < .05), fraction_p_below_001 = mean(p < .01),
      FDR_below_005 = sum(r$FDR_significant), dispersion_fit_type = fit,
      min_dispersion_boundary = min_disp, finite_initial_estimates = n_est,
      initial_at_boundary = n_floor, percent_initial_at_boundary = 100 * n_floor / n_est,
      median_final_dispersion = stats::median(mc$dispersion, na.rm = TRUE))
    h <- graphics::hist(p, breaks = seq(0, 1, by = .05), plot = FALSE)
    plot_device(file.path(diagnostic_dir, paste0(population, "_pvalues")), function() {
      graphics::par(mar = c(5, 5, 5.5, 2), mfrow = c(1, 1))
      graphics::plot(h, freq = FALSE, col = "#91BEE2", border = "white",
        main = paste(gsub("_", " ", population), "- raw p values"),
        xlab = "Raw p value", ylab = "Density", xlim = c(0, 1),
        ylim = c(0, max(1.2, h$density * 1.12)))
      graphics::abline(h = 1, col = "#D55E00", lty = 2, lwd = 2)
      graphics::mtext(sprintf("%s finite tests; p < 0.01: %s (%.2f%%); FDR < 0.05: %s",
        format(length(p), big.mark = ","), sum(p < .01), 100 * mean(p < .01),
        sum(r$FDR_significant)), side = 3, line = .5, cex = .85)
      graphics::mtext("Dashed line: uniform reference, not a required shape or a calibration test.",
        side = 1, line = 3.8, cex = .75)
    })
    plot_device(file.path(diagnostic_dir, paste0(population, "_dispersion")), function() {
      graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(4.5, 1.65))
      graphics::par(mar = c(4, 5, 4.5, 2))
      DESeq2::plotDispEsts(dds, main = paste(gsub("_", " ", population), "- dispersion"),
        xlab = "Mean normalized peak counts (log scale)", ylab = "Dispersion (log scale)",
        genecol = "#77777740", fitcol = "#D55E00", finalcol = "#0072B270",
        cex = .28, legend = FALSE, xlim = xlim, ylim = ylim)
      graphics::legend("topright", c("Initial peak-wise estimate", "Fitted trend", "Final estimate"),
        col = c("#777777", "#D55E00", "#0072B2"), pch = 16, bty = "n", cex = .8)
      graphics::mtext(paste("Trend fit:", fit, "| Same axis limits across the five dispersion figures"),
        side = 3, line = .4, cex = .75)
      graphics::par(mar = c(3.6, 5, 2.8, 2))
      graphics::plot.new(); graphics::plot.window(xlim = c(0, 100), ylim = c(0, 1), xaxs = "i")
      pct <- 100 * n_floor / n_est
      graphics::rect(0, .25, pct, .65, col = "#777777", border = NA)
      graphics::rect(pct, .25, 100, .65, col = "#9DC9E2", border = NA)
      graphics::axis(1, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"))
      graphics::title(main = "Initial estimates at the numerical boundary (1e-8)", cex.main = 1)
      graphics::legend("top", legend = c(sprintf("At boundary: %s (%.1f%%)", n_floor, pct),
        sprintf("Above boundary: %s (%.1f%%)", n_est - n_floor, 100 - pct)),
        fill = c("#777777", "#9DC9E2"), bty = "n", horiz = TRUE, cex = .85)
      graphics::mtext("Initial estimates only. Boundary does not mean no signal or a failed model; final estimates differ.",
        side = 1, line = 2.1, cex = .73)
    })
    # Count stratification is numeric supplementary data, not another crowded plot.
    breaks <- c(-Inf, 10, 50, 100, Inf)
    bins <- cut(mc$baseMean, breaks = breaks, right = FALSE,
                labels = c("<10", "10_to_<50", "50_to_<100", ">=100"))
    for (bin in levels(bins)) {
      ix <- which(bins == bin)
      ixr <- match(rownames(dds)[ix], r$peak)
      strata_rows[[paste(population, bin)]] <- data.frame(population = population,
        mean_count_bin = bin, model_peaks = length(ix),
        finite_pvalues = sum(is.finite(r$pvalue[ixr])),
        raw_p_below_001 = sum(r$pvalue[ixr] < .01, na.rm = TRUE),
        initial_at_boundary = sum(floor[ix]),
        percent_initial_at_boundary = if (length(ix)) 100 * mean(floor[ix]) else NA_real_)
    }
  }
  diagnostics <- bind(diag_rows)
  write_csv(diagnostics, file.path(diagnostic_dir, "diagnostic_summary.csv"))
  write_csv(bind(strata_rows), file.path(diagnostic_dir, "count_stratified_diagnostics.csv"))

  # One predeclared broader filter, not a search for a more significant result.
  sensitivity <- list(); comparisons <- list(); shared_rows <- list(); sample_rows <- list()
  for (population in populations) {
    message("Broader count filter: ", population)
    m <- raw_matrices[[population]]
    m <- m[rowSums(m >= exploratory_count) >= 2L, , drop = FALSE]
    md <- data.frame(sample = samples, condition = factor(conditions, c("control", "stress")),
                      row.names = samples)
    dds <- DESeq2::DESeqDataSetFromMatrix(m, md, ~ condition)
    dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
    pop_dir <- file.path(sensitivity_dir, population); dir.create(pop_dir)
    saveRDS(dds, file.path(pop_dir, "deseq2_model.rds"))
    r <- extract_results(dds, population)
    write_csv(r, file.path(pop_dir, "all_peak_results.csv"))
    write_csv(r[r$FDR_significant, , drop = FALSE], file.path(pop_dir, "FDR_significant_peak_results.csv"))
    sensitivity[[population]] <- r
    old <- primary[[population]]
    ix <- match(old$peak, r$peak)
    stopifnot(!anyNA(ix), identical(
      unname(as.matrix(old[, paste0("raw_", samples)])),
      unname(as.matrix(r[ix, paste0("raw_", samples)]))))
    for (setting in c("primary_count10", paste0("exploratory_count", exploratory_count))) {
      z <- if (setting == "primary_count10") old else r
      comparisons[[paste(population, setting)]] <- data.frame(population = population,
        setting = setting, min_count = if (setting == "primary_count10") 10L else exploratory_count,
        min_samples = 2L, peaks_in_model = nrow(z), finite_pvalues = sum(is.finite(z$pvalue)),
        adjusted_pvalues = sum(is.finite(z$padj)), FDR_significant = sum(z$FDR_significant),
        significant_higher_2fold = sum(z$FDR_significant & z$log2FoldChange >= 1, na.rm = TRUE),
        significant_lower_2fold = sum(z$FDR_significant & z$log2FoldChange <= -1, na.rm = TRUE))
    }
    shared_rows[[population]] <- data.frame(population = population,
      primary_peaks = nrow(old), additional_peaks = nrow(r) - nrow(old),
      original_significant_retained = sum(old$FDR_significant & r$FDR_significant[ix]),
      primary_only_significant = sum(old$FDR_significant & !r$FDR_significant[ix]),
      newly_significant = sum(r$FDR_significant & !r$peak %in% old$peak[old$FDR_significant]),
      median_abs_log2FC_change_shared = stats::median(abs(r$log2FoldChange[ix] - old$log2FoldChange)),
      max_abs_log2FC_change_shared = max(abs(r$log2FoldChange[ix] - old$log2FoldChange)),
      shared_raw_counts_identical = TRUE)
    sample_rows[[population]] <- data.frame(population = population, sample = samples,
      primary_size_factor = DESeq2::sizeFactors(models[[population]]),
      exploratory_size_factor = DESeq2::sizeFactors(dds), row.names = NULL)
  }
  comparison <- bind(comparisons)
  write_csv(comparison, file.path(sensitivity_dir, "threshold_comparison.csv"))
  write_csv(bind(shared_rows), file.path(sensitivity_dir, "shared_peak_comparison.csv"))
  write_csv(bind(sample_rows), file.path(sensitivity_dir, "size_factor_comparison.csv"))

  # RNA-only selection, frozen before looking up ATAC: top three significant
  # genes per population, ordered by RNA padj, |log2FC|, then gene ID.
  selected_rows <- list()
  for (i in seq_along(populations)) {
    population <- populations[i]
    r <- source_csv(rna_paths[i])
    stopifnot(all(c("gene", "padj", "log2FoldChange", "cell_population", "branch") %in% names(r)),
      !anyDuplicated(r$gene), all(r$branch == "pca"), all(r$cell_population == population))
    r <- r[which(is.finite(r$padj) & r$padj < .05 & is.finite(r$log2FoldChange) &
                   abs(r$log2FoldChange) >= 1), , drop = FALSE]
    r <- r[order(r$padj, -abs(r$log2FoldChange), r$gene), , drop = FALSE]
    r <- head(r, top_rna_genes)
    stopifnot(nrow(r) > 0L)
    selected_rows[[population]] <- data.frame(population = population, gene_id = r$gene,
      RNA_rank = seq_len(nrow(r)), RNA_log2FC = r$log2FoldChange,
      RNA_lfcSE = r$lfcSE, RNA_pvalue = r$pvalue, RNA_padj = r$padj,
      RNA_branch = "pca", stringsAsFactors = FALSE)
  }
  selected <- bind(selected_rows)
  common <- readRDS(common_path)
  common_ids <- paste0(as.character(GenomicRanges::seqnames(common)), "-",
                       GenomicRanges::start(common), "-", GenomicRanges::end(common))
  stopifnot(!anyDuplicated(common_ids))
  ann <- rtracklayer::import(gtf_path, format = "gtf")
  genes <- ann[which(ann$type == "gene" & as.character(GenomicRanges::seqnames(ann)) %in% as.character(1:5))]
  stopifnot(!anyDuplicated(as.character(genes$gene_id)))
  k <- match(selected$gene_id, as.character(genes$gene_id))
  stopifnot(!anyNA(k))
  selected$gene_name <- as.character(genes$gene_name[k])
  missing_name <- is.na(selected$gene_name) | selected$gene_name == ""
  selected$gene_name[missing_name] <- selected$gene_id[missing_name]
  write_csv(selected, file.path(response_dir, "RNA_selected_genes.csv"))
  long_rows <- list(); gene_rows <- list()
  for (i in seq_len(nrow(selected))) {
    s <- selected[i, , drop = FALSE]; g <- genes[k[i]]
    window <- GenomicRanges::promoters(g, upstream = 2000L, downstream = 500L)
    GenomicRanges::start(window) <- pmax(1L, GenomicRanges::start(window))
    body <- IRanges::overlapsAny(common, g, ignore.strand = TRUE)
    promoter <- IRanges::overlapsAny(common, window, ignore.strand = TRUE)
    hits <- which(body | promoter)
    tss <- if (as.character(GenomicRanges::strand(g)) == "+")
      GenomicRanges::start(g) else GenomicRanges::end(g)
    map <- data.frame(peak = common_ids[hits], overlaps_gene_body = body[hits],
      overlaps_TSS_window = promoter[hits], TSS_distance_bp = pmax(
        GenomicRanges::start(common)[hits] - tss, tss - GenomicRanges::end(common)[hits], 0))
    r <- primary[[s$population]]
    ri <- match(map$peak, r$peak)
    if (nrow(map)) {
      at <- r[ri, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
                    "FDR_significant", "replicate_pattern", paste0("normalized_", samples),
                    paste0("raw_", samples)), drop = FALSE]
      names(at) <- paste0("ATAC_", names(at))
      d <- cbind(s[rep(1L, nrow(map)), , drop = FALSE], map, at)
      d$ATAC_in_model <- !is.na(ri)
      d$direction_relation <- "Not_modelled_or_missing_estimate"
      ok <- d$ATAC_in_model & is.finite(d$ATAC_log2FoldChange)
      d$direction_relation[ok] <- ifelse(sign(d$RNA_log2FC[ok]) == sign(d$ATAC_log2FoldChange[ok]),
                                        "Same_estimated_direction", "Different_estimated_direction")
      long_rows[[i]] <- d
    }
    present <- which(!is.na(ri))
    # Representative display peak is selected by promoter overlap, then distance,
    # then coordinate ID. No ATAC p value or effect size enters this choice.
    representative <- if (length(present)) {
      present[order(!map$overlaps_TSS_window[present], map$TSS_distance_bp[present], map$peak[present])][1]
    } else NA_integer_
    gene_rows[[i]] <- data.frame(s, common_peaks = nrow(map), modelled_peaks = length(present),
      nearby_FDR_peaks = sum(r$FDR_significant[ri], na.rm = TRUE),
      representative_peak = if (is.na(representative)) NA_character_ else map$peak[representative],
      status = if (!nrow(map)) "No_common_peak_match" else if (!length(present))
        "Matching_peaks_not_modelled" else "Existing_ATAC_results_available")
  }
  long <- bind(long_rows); gene_summary <- bind(gene_rows)
  write_csv(long, file.path(response_dir, "RNA_ATAC_peak_context.csv"))
  write_csv(gene_summary, file.path(response_dir, "RNA_ATAC_gene_summary.csv"))
  response_summary <- bind(lapply(populations, function(population) {
    g <- gene_summary[gene_summary$population == population, , drop = FALSE]
    d <- long[long$population == population & long$ATAC_in_model, , drop = FALSE]
    data.frame(population = population, selected_RNA_genes = nrow(g),
      genes_with_common_peak = sum(g$common_peaks > 0),
      genes_with_modelled_peak = sum(g$modelled_peaks > 0),
      genes_with_nearby_FDR_peak = sum(g$nearby_FDR_peaks > 0),
      modelled_gene_peak_pairs = nrow(d),
      same_estimated_direction_pairs = sum(d$direction_relation == "Same_estimated_direction"),
      different_estimated_direction_pairs = sum(d$direction_relation == "Different_estimated_direction"))
  }))
  write_csv(response_summary, file.path(response_dir, "response_summary.csv"))
  core2_plot_response_context(response_dir, populations)

  rna_diag <- source_csv(file.path(project_dir, "results", "core1", "tables", "pseudobulk_sample_diagnostics.csv"))
  rna_diag <- rna_diag[rna_diag$branch == "pca" & rna_diag$population %in% populations, , drop = FALSE]
  atac_md <- pb$metadata[pb$metadata$rna_population %in% populations, , drop = FALSE]
  count_scope <- merge(rna_diag[, c("population", "sample", "n_cells")],
    atac_md[, c("rna_population", "sample", "n_nuclei")],
    by.x = c("population", "sample"), by.y = c("rna_population", "sample"), all = TRUE)
  names(count_scope)[3:4] <- c("RNA_nuclei_existing_analysis", "ATAC_nuclei_paired_subset")
  stopifnot(nrow(count_scope) == 20L, !anyNA(count_scope))
  write_csv(count_scope, file.path(response_dir, "retained_nucleus_scope.csv"))

  hashes_after <- tools::md5sum(inputs)
  stopifnot(identical(hashes_before, hashes_after))
  write_csv(data.frame(path = inputs, md5_before = unname(hashes_before),
    md5_after = unname(hashes_after), unchanged = unname(hashes_before == hashes_after)),
    file.path(output_dir, "original_input_integrity.csv"))
  write_csv(bind(validation), file.path(output_dir, "primary_model_validation.csv"))
  saveRDS(list(populations = populations, samples = samples, design = "~ condition",
    primary_filter = "counts >= 10 in >= 2 samples", exploratory_count = exploratory_count,
    min_samples = 2L, alpha = .05, FDR_scope = "BH within each population and analysis",
    min_dispersion_boundary = min_disp, original_inputs_unchanged = TRUE,
    RNA_rule = paste("pca branch, padj < .05, |log2FC| >= 1; top", top_rna_genes,
      "per population ordered by RNA padj, descending |log2FC|, gene ID"),
    mapping_rule = "Gene body OR strand-aware TSS upstream 2000/downstream 500 bp",
    representative_rule = "Modelled peaks: promoter overlap first, then TSS distance, then peak ID",
    no_new_RNA_fit = TRUE, no_nucleus_filter_change = TRUE, no_samples_removed = TRUE,
    no_Control1_attribution_experiment = TRUE), file.path(output_dir, "parameters.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  writeLines(c("# Core 2 supplemental review", "",
    "Original counts>=10 results remain primary; counts>=5 is a declared broader sensitivity analysis.",
    "The same four biological samples remain the inference units; no nucleus is a biological replicate.",
    "model_diagnostics contains ten standalone figures (each PDF and PNG); dispersion figures include numerical boundary proportions.",
    "A uniform p-value reference is not a requirement; histograms alone do not establish miscalibration.",
    "Initial dispersion boundary estimates are not final estimates and do not prove model failure or Control 1 causality.",
    "RNA response genes were selected from existing primary-branch RNA statistics before ATAC lookup.",
    "RNA and ATAC use different QC-retained nucleus sets; no new RNA model was fitted.",
    "All mapped peaks are retained, including unmodelled ones. Proximity is not regulatory target validation.",
    "Root_cap remains exploratory. FDR is population-wise, not experiment-wide or across both filters.",
    "No Control 1 equal-nucleus sampling, sample deletion, TSS sensitivity or upstream reconstruction was performed.",
    "This local revision has not been committed or pushed."), file.path(output_dir, "README.md"))
  print(comparison, row.names = FALSE)
  print(response_summary, row.names = FALSE)
  invisible(list(output_dir = output_dir, diagnostics = diagnostics,
                 sensitivity = comparison, RNA_ATAC = response_summary))
}
