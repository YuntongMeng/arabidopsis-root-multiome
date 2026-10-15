# Synthetic checks only. Run from the project root; no project data are loaded.
root <- Sys.getenv("P3_ROOT", getwd())
source(file.path(root, "scripts", "core2", "00_common.R"))
library <- Sys.getenv("P3_R_LIB", "")
if (nzchar(library)) .libPaths(c(library, .libPaths()))

metadata <- data.frame(nuclear_fragments = c(1000, 999, 1500, 1500, 1500, 1500),
  TSS.enrichment = c(1.1, 2, 1, 2, 2, 2), nucleosome_signal = c(1, 1, 1, 2, 1, 1),
  row.names = paste0("b", 1:6))
calls <- data.frame(cell = paste0("control_rep1_b", 1:6), sample = "control_rep1",
  qc_pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  scDblFinder_class = c("singlet", "singlet", "singlet", "singlet", "doublet", NA), scDblFinder_score = 0)
selection <- core2_selection(metadata, "control_rep1", calls)
stopifnot(identical(which(selection$pass_paired_qc), 1L),
          identical(which(selection$pass_basic_atac), c(1L, 5L, 6L)))

core2_require("Matrix")
counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3, 0, 4, 5, 6), nrow = 2), sparse = TRUE)
rownames(counts) <- c("p1", "p2")
cells <- c("c1", "c2", "c3", "c4")
attr(cells, ".match.hash") <- "fixture"
colnames(counts) <- cells
metadata <- data.frame(sample = c("s1", "s1", "s2", "s2"),
  condition = c("control", "control", "stress", "stress"),
  rna_population = c("A", "B", "A", "A"), nuclear_fragments = c(10, 20, 30, 40),
  percent_organelle = 0, row.names = c("c1", "c2", "c3", "c4"))
pb <- core2_pseudobulk_counts(counts, metadata[c(4, 1, 3, 2), ])
stopifnot(ncol(pb$counts) == 3L, all(as.numeric(pb$counts[, "s2__A"]) == c(5, 10)),
          sum(pb$counts) == sum(counts), !("s2__B" %in% colnames(pb$counts)))
bad <- try(core2_align_metadata(counts, metadata[-1, ]), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
cat("Synthetic selection/alignment/pseudobulk checks: PASS\n")
