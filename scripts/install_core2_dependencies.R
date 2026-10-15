#!/usr/bin/env Rscript
# Explicit installation only; default/help does not install or run analysis.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] == "help") {
  cat("Usage: Rscript scripts/install_core2_dependencies.R install\n",
      "Optional P3_R_LIB selects the installation library.\n",
      "Installs missing packages; does not lock recorded versions or install upstream tools.\n", sep = "")
} else {
  stopifnot(identical(args[1], "install"))
  library <- Sys.getenv("P3_R_LIB", "")
  if (nzchar(library)) {
    dir.create(library, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(library, .libPaths()))
  }
  cran <- c("Seurat", "SeuratObject", "Signac", "Matrix", "ggplot2", "patchwork", "ggrepel", "BiocManager")
  bioc <- c("GenomicRanges", "IRanges", "S4Vectors", "Rsamtools", "rtracklayer", "DESeq2", "GenomeInfoDb", "SummarizedExperiment")
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) utils::install.packages(missing, repos = "https://cloud.r-project.org")
  missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)
  stopifnot(all(vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)))
  cat("Core 2 R dependencies available. Review session versions before reproducing results.\n")
}
