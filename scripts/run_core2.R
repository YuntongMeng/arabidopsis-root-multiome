#!/usr/bin/env Rscript
# Explicit runner. "help" needs no packages, checkpoints or workspace settings.
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "help"
allowed <- c("help", "select", "embeddings", "annotation", "markers", "pseudobulk", "diagnostics", "da", "p1", "pca")
if (!stage %in% allowed) stop("Unknown stage: ", stage)
if (stage == "help") {
  cat("Core 2 stages (run explicitly; no implicit all):\n",
      "select      ATAC thresholds + exact RNA QC/singlet pairing; merge\n",
      "embeddings  TF-IDF, LSI, UMAP and clustering\n",
      "annotation  Exact RNA identity matching + gene activity\n",
      "markers     Optional marker dot plot and SCR coverage\n",
      "pseudobulk  Sum raw peak counts by sample x RNA population\n",
      "diagnostics Fit/reuse four-sample DESeq2 models; sample diagnostics\n",
      "da          Export differential accessibility from saved models\n",
      "p1          Match all P1 candidates to nearby common peaks\n",
      "pca         Redraw saved PCA coordinates with equal x/y unit scale\n",
      "Set P3_ROOT, P3_CORE2_WORK and optionally P3_R_LIB/P3_CORE2_ANALYSIS.\n",
      "Default destination: workspace/outputs/core2_rebuild (separate from reviewed run).\n",
      "See docs/core2/running.md. No analysis is run by help.\n", sep = "")
} else {
  root <- normalizePath(Sys.getenv("P3_ROOT", getwd()), mustWork = TRUE)
  code <- file.path(root, "scripts", "core2")
  for (file in c("00_common.R", "01_pair_and_merge.R", "02_embeddings.R",
                 "03_annotation_and_activity.R", "04_pseudobulk.R",
                 "05_sample_diagnostics.R", "06_export_differential_accessibility.R",
                 "07_P1_contextualization.R", "08_redraw_PCA.R")) {
    sys.source(file.path(code, file), envir = .GlobalEnv)
  }
  config <- core2_config()
  switch(stage,
    select = core2_pair_and_merge(config),
    embeddings = core2_embeddings(config),
    annotation = core2_annotation_and_activity(config),
    markers = core2_marker_and_coverage_plots(config),
    pseudobulk = core2_pseudobulk(config),
    diagnostics = core2_sample_diagnostics(config$analysis),
    da = core2_export_differential_accessibility(config$analysis),
    p1 = core2_P1_contextualization(config$work, config$root, config$analysis),
    pca = core2_redraw_PCA(config$analysis, show_plot = FALSE))
}
