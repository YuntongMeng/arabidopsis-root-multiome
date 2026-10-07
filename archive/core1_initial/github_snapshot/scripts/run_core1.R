# Full Core 1 reproduction from original matrices.
# Rscript scripts/run_core1.R (run from the project root)
# source() defines the wrapper only; call run_core1_full() explicitly in RStudio.
run_core1_full <- function(root = getwd()) {
  root <- normalizePath(root, mustWork = TRUE)
  output <- file.path("results", paste0("core1_rerun_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  pipeline <- new.env(parent = globalenv())
  for (script in c("00_core1_common.R", "01_core1_rna_preprocessing_clustering.R",
                   "02_core1_cell_population_annotation.R", "03_core1_pseudobulk_stress_DE.R",
                   "04_core1_P1_contextualization.R")) {
    sys.source(file.path(root, "scripts", script), envir = pipeline)
  }
  message("Stage 1/4: RNA preprocessing and clustering")
  invisible(pipeline$run_core1_01(root = root, output = output))
  invisible(gc())
  message("Stage 2/4: markers, signatures and reviewed annotation")
  invisible(pipeline$run_core1_02(root = root, output = output,
    checkpoint = file.path(output, "checkpoints/01_clustered.rds"), recompute_evidence = TRUE))
  invisible(gc())
  message("Stage 3/4: biological-sample pseudobulk DE")
  invisible(pipeline$run_core1_03(root = root, output = output,
    checkpoint = file.path(output, "checkpoints/02_annotated.rds"), recompute = TRUE))
  invisible(gc())
  message("Stage 4/4: P1 candidate and GO contextualization")
  invisible(pipeline$run_core1_04(root = root, output = output,
    de_file = file.path(output, "tables/pseudobulk_DE_all.csv"), use_frozen_mapping = TRUE))
  message("Core 1 completed. Outputs: ", file.path(root, output))
  invisible(output)
}
if (sys.nframe() == 0L) run_core1_full()
