# Core 2: export differential accessibility from existing pseudobulk models.
# source() only defines the function. No models are fitted or samples removed.
# Significance tests H0: log2FC = 0. The twofold threshold below is a separate
# point-estimate prioritization criterion, not a test of |log2FC| > 1.
# BH correction is performed within each population, not across the experiment.

core2_export_differential_accessibility <- function(
    analysis_dir,
    populations = c(
      "Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage",
      "Pericycle", "Root_cap"
    ),
    alpha = 0.05,
    priority_abs_log2fc = 1) {

  stopifnot(
    requireNamespace("DESeq2", quietly = TRUE),
    is.finite(alpha), alpha > 0, alpha < 1,
    is.finite(priority_abs_log2fc), priority_abs_log2fc > 0,
    length(populations) > 0, !anyDuplicated(populations)
  )
  plain_ids <- function(x) {
    x <- as.character(x)
    attributes(x) <- NULL
    x
  }
  sample_order <- c(
    "control_rep1", "control_rep2", "stress_rep1", "stress_rep2"
  )
  expected_condition <- c("control", "control", "stress", "stress")
  out_dir <- file.path(analysis_dir, "differential_accessibility")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summaries <- list()

  for (population in populations) {
    message("Exporting population: ", population)
    model_path <- file.path(
      analysis_dir, "sample_diagnostics", population, "deseq2_model.rds"
    )
    stopifnot(file.exists(model_path))
    dds <- readRDS(model_path)
    stopifnot(
      identical(plain_ids(colnames(dds)), sample_order),
      identical(as.character(dds$condition), expected_condition),
      identical(all.vars(DESeq2::design(dds)), "condition"),
      !is.null(DESeq2::sizeFactors(dds))
    )
    results <- DESeq2::results(
      object = dds,
      contrast = c("condition", "stress", "control"),
      alpha = alpha,
      pAdjustMethod = "BH",
      independentFiltering = TRUE
    )
    r <- as.data.frame(results)
    normalized <- DESeq2::counts(dds, normalized = TRUE)[
      rownames(r), sample_order, drop = FALSE
    ]
    raw <- DESeq2::counts(dds, normalized = FALSE)[
      rownames(r), sample_order, drop = FALSE
    ]
    stopifnot(all(is.finite(normalized)), all(is.finite(raw)))
    norm_table <- as.data.frame(normalized)
    raw_table <- as.data.frame(raw)
    colnames(norm_table) <- paste0("normalized_", sample_order)
    colnames(raw_table) <- paste0("raw_", sample_order)
    table <- data.frame(
      peak = plain_ids(rownames(r)), r, norm_table, raw_table,
      row.names = NULL, check.names = FALSE
    )
    table$population <- population
    table$analysis_status <- if (population == "Root_cap") {
      "Exploratory_low_sample_reproducibility"
    } else {
      "Replicate_and_region_review_required"
    }
    table$FDR_significant <- !is.na(table$padj) & table$padj < alpha
    table$priority_change <- table$FDR_significant &
      is.finite(table$log2FoldChange) &
      abs(table$log2FoldChange) >= priority_abs_log2fc
    table$category <- "Not_selected"
    table$category[
      table$priority_change & table$log2FoldChange >= priority_abs_log2fc
    ] <- "Higher_in_stress"
    table$category[
      table$priority_change & table$log2FoldChange <= -priority_abs_log2fc
    ] <- "Lower_in_stress"
    table$replicate_pattern <- "Groups_overlap"
    control_min <- pmin(normalized[, 1], normalized[, 2])
    control_max <- pmax(normalized[, 1], normalized[, 2])
    stress_min <- pmin(normalized[, 3], normalized[, 4])
    stress_max <- pmax(normalized[, 3], normalized[, 4])
    table$replicate_pattern[stress_min > control_max] <-
      "Both_stress_above_both_control"
    table$replicate_pattern[stress_max < control_min] <-
      "Both_stress_below_both_control"
    table$direction_supported <-
      (table$category == "Higher_in_stress" &
         table$replicate_pattern == "Both_stress_above_both_control") |
      (table$category == "Lower_in_stress" &
         table$replicate_pattern == "Both_stress_below_both_control")
    table <- table[order(table$padj, na.last = TRUE), , drop = FALSE]
    rownames(table) <- NULL

    population_dir <- file.path(out_dir, population)
    dir.create(population_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(table, file.path(population_dir, "all_peak_results.csv"),
                     row.names = FALSE)
    utils::write.csv(
      table[table$FDR_significant, , drop = FALSE],
      file.path(population_dir, "FDR_significant_peak_results.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      table[table$priority_change, , drop = FALSE],
      file.path(population_dir, "priority_peak_results.csv"),
      row.names = FALSE
    )
    summaries[[population]] <- data.frame(
      population = population,
      peaks_in_model = nrow(table),
      peaks_with_adjusted_p = sum(!is.na(table$padj)),
      FDR_significant = sum(table$FDR_significant),
      higher_in_stress_2fold = sum(table$category == "Higher_in_stress"),
      lower_in_stress_2fold = sum(table$category == "Lower_in_stress"),
      priority_with_direction_support = sum(table$direction_supported),
      analysis_status = unique(table$analysis_status),
      row.names = NULL
    )
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  utils::write.csv(summary, file.path(out_dir, "differential_accessibility_summary.csv"),
                   row.names = FALSE)
  saveRDS(list(
    populations = populations, contrast = c("condition", "stress", "control"),
    null_hypothesis = "log2FC = 0", alpha = alpha,
    priority_abs_log2fc = priority_abs_log2fc,
    fold_change_type = "Unshrunken model estimate",
    FDR_scope = "BH within each population, not across all populations",
    replicate_pattern = "Descriptive comparison of all four normalized sample counts",
    samples_removed = character(0),
    Root_cap_status = "Exploratory_low_sample_reproducibility"
  ), file.path(out_dir, "export_parameters.rds"))
  writeLines(c(
    "# Core 2 differential accessibility exports", "",
    "All results use raw peak pseudobulk models with design ~ condition.",
    "Contrast: stress / control. Positive log2FC means relatively higher accessibility in stress.",
    "Each column/sample is a biological replicate, not an individual nucleus.",
    paste0("BH-adjusted p < ", alpha, " is evaluated separately within each population."),
    "This is not a single experiment-wide FDR guarantee over all populations.",
    paste0("Priority peaks additionally have estimated |log2FC| >= ", priority_abs_log2fc, "."),
    "That effect-size filter is not a hypothesis test that the true change exceeds the threshold.",
    "Not_selected and missing adjusted p values do not establish absence of biological change.",
    "All modeled peaks and the four raw and normalized counts are preserved.",
    "Direction support is descriptive; it neither replaces the model nor filters the saved priority list.",
    "Root_cap has weak within-condition sample reproducibility and is labeled exploratory.",
    "The other populations also require replicate/region review before biological claims.",
    "These are genomic regions, not measured gene expression or confirmed regulatory targets."
  ), file.path(out_dir, "README.md"))
  display <- summary[, c(
    "population", "peaks_in_model", "FDR_significant",
    "higher_in_stress_2fold", "lower_in_stress_2fold"
  ), drop = FALSE]
  print(display, row.names = FALSE)
  invisible(summary)
}
