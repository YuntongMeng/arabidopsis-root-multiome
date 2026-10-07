# Core 1 full rerun verification

The user ran all four stages from the original RNA matrices on 2026-10-05. The completed run directory was `results/core1_rerun_20261005_230515/`; its tables, figures and diagnostics are published under `results/core1/`.

| Check | Result |
| --- | --- |
| Retained nuclei | 24,277 |
| PC pairwise ARI matrix | Matches original matrix |
| Reviewed cluster annotation table | Matches original; manual labels reused |
| 13-population DE summary: genes_tested, DEGs, up, down | All match |
| P1 candidate summary: denominators, significant counts, directions and percentages | All match |
| Complete GO × population matrix | All match |
| Stage checkpoints | 01_clustered.rds, 02_annotated.rds, 03_pseudobulk_DE.rds created locally |

Table comparisons matched rows by their identifiers and used a numerical tolerance of 1e-9 relative / 1e-10 absolute. Checks above establish reproducibility of the specified key outputs; every gene-level statistic was not separately compared.

The initial run encountered a strict `vapply(..., integer(1))` type error because Seurat dimensions returned double. The four dimension calls were changed to `numeric(1)`; the subsequent complete run finished successfully.

Clustering, strict markers, module scores and pseudobulk DE were recomputed. Reviewed manual annotation labels and the resolved GO mapping were reused. Annotation provenance limitations remain documented in the analysis summary. Large full Seurat checkpoints remain local; the lightweight archived result bundle is included.
