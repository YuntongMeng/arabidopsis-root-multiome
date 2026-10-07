# Project 3 — Core 1: corrected RNA analysis under osmotic stress

[Project 3 overview](../../README.md) · [Core 2 — planned](../core2/README.md) · [Core 3 — planned](../core3/README.md)

Core 1 studies wild-type root-tip transcription under sorbitol osmotic stress and contextualizes Project 1 BRL3 × drought candidates. The first version mistakenly included ATAC peaks in its RNA assay. Its scripts and results are preserved as historical evidence; the corrected analysis starts with gene-expression rows only.

## Key findings

- Each sample supplies 32,833 gene features; zero coordinate-shaped peak IDs remain after input filtering.
- Gene-only QC and per-capture scDblFinder retain 21,550 predicted singlets.
- Unintegrated PCA is primary; sample-based Harmony is a sensitivity branch. At 25 dimensions there are 18 and 17 clusters respectively, yielding nine and ten DE-eligible broad populations under the same >=50 nuclei per biological sample rule. These counts are observations, not required answers.
- Pericycle has 19 significant P1 candidates in both branches (120/117 fitted candidates respectively). P1 candidates are enriched relative to ordinary testable genes, but no population shows BH-significant extra enrichment relative to the P1 drought-DEG control.
- Seven GO views summarize 13 unique overlapping candidate genes, not new P3 GO enrichment tests. Only three terms have nonzero P3 DEG overlap.

## Analysis overview

| Script | Purpose |
| --- | --- |
| [00](../../scripts/core1/00_common.R) | Paths, shared helpers and sample metadata |
| [01](../../scripts/core1/01_input_audit.R) | Feature-type audit and gene-only matrices |
| [02](../../scripts/core1/02_qc_doublets.R) | Distribution-driven per-sample QC and predicted doublets |
| [03](../../scripts/core1/03_embeddings.R) | PCA/Harmony and dimension robustness grid |
| [04](../../scripts/core1/04_branch_evidence.R) | Replicate mixing, condition geometry and sourced marker evidence |
| [05](../../scripts/core1/05_annotation_pseudobulk.R) | Reviewed labels, eligibility and raw-count DESeq2 pseudobulk |
| [06](../../scripts/core1/06_P1_context.R) | Candidate Fisher tests, drought controls and GO overlap |
| [07](../../scripts/core1/07_compare_validate.R) | Historical comparison and structural checks |

The former four large scripts are split for teaching and inspection. Sourcing defines functions; invoking each stage recomputes its work. The [runner](../../scripts/run_core1_corrected.R) makes this explicit. Visualization code is [09_render_figures.R](../../exploration/core1/09_render_figures.R); optional historical diagnostics are separate.

## Running the workflow

Restore the four GEO inputs using [data sources](../../docs/data_sources.md), install dependencies, and start in the project root:

```sh
Rscript scripts/install_core1_dependencies.R
Rscript scripts/run_core1_corrected.R all
```

`P3_ROOT` and `P3_R_LIB` are optional overrides. Reviewed annotation is bound to the partition MD5; a changed partition stops DE for renewed marker review. Current package versions are recorded, but cross-version bitwise reproduction is not guaranteed. No old cluster/population/DE count is hardcoded as a success criterion.

## Outputs and decision trail

[Tables](../../results/core1/tables/) · [Main figures](../../results/core1/figures/) · [Supplementary diagnostics](../../results/core1/diagnostics/) · [Manual rerun validation](../../docs/core1/manual_rerun_validation.md)

[Analysis story](../../docs/core1/analysis_story.md), [decision log](../../docs/core1/decision_log.md), and [version comparison](../../docs/core1/version_comparison.md) record the error, alternatives and updated conclusions. [The initial archive](../../archive/core1_initial/README.md) and frozen outputs are historical, not recommended run paths.

## Interpretation boundary

Two biological replicates per condition limit inference. Broad lineage, low-confidence and cell-state labels remain visible. Harmony improves replicate mixing while attenuating condition-associated geometry; this is not proof of biological overcorrection. Raw counts remain unchanged for DE. No verified empty-droplet input is available, so ambient correction is not forced. P1/P3 overlap gives stress context, not cell-specific BRL3 interaction validation; zero GO overlap is not absence of a biological program.
