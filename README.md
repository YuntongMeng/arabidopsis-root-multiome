# Arabidopsis root-tip multiome under osmotic stress

Project 3 investigates cell-population-specific osmotic-stress responses in Arabidopsis root tips and places Project 1 BRL3 × drought interaction candidates into an independent WT cellular context.

**Current release: Core 1 (RNA) only.** The complete RNA workflow has been rerun from the original input matrices, and its key outputs reproduce the archived analysis.

## Three-core project structure

| Core | Scope in this release | Status |
| --- | --- | --- |
| Core 1 | RNA QC, robust clustering, conservative population annotation, biological-sample pseudobulk DE, and P1 candidate/GO contextualization | Complete; full rerun passed |
| Core 2 | Reserved for the next core of the original project plan; detailed scope to be restored from the agreed plan before implementation | Planned; no analysis released |
| Core 3 | Reserved for the final core of the original project plan; detailed scope to be restored from the agreed plan before implementation | Planned; no analysis released |

The two future-core descriptions above deliberately do not infer missing planning details. This snapshot contains only Core 1 analysis code and results.

## Core 1 workflow and decisions

Raw RNA matrices → input audit → distribution-based QC → LogNormalize/vst/PCA → compare PC10/15/20/25/30 by pairwise ARI → PC1:25 clustering → evidence-guided population annotation → population × biological-sample pseudobulk → stress/control DESeq2 → P1 contextualization.

The scripts preserve the decision trail, including the PC20-to-PC25 adjustment, cluster 6/14 sample enrichment checks, Seurat v5 JoinLayers repair, cluster 12 diagnostic-only relaxed markers, conservative annotation labels, pseudobulk name-matching repair, and the four-sample shared feature universe. Ambiguous R functions use explicit package namespaces.

- **24,277 nuclei**, **18 transcriptional clusters**, **17 annotated populations**.
- **13 eligible populations** with ≥25 nuclei in each of the four biological samples.
- **32,833 common features** before per-population count prefiltering.
- Formal DE: `padj < 0.05` and `abs(log2FoldChange) >= 1`; positive means stress-up.
- P1 input: **137 BRL3 × drought interaction candidates**; 129 occur in at least one P3 population result table.
- Pericycle has the largest P1 candidate response fraction: **23/118 = 19.49%**, with 6 up and 17 down.
- Seven overlapping P1 GO annotations are contextualized; Pericycle appears in all seven, with three cold-response genes downregulated.

P1/P3 differ in genotype and stress conditions. These results provide **contextualization**, not validation of BRL3 × drought interaction effects. Shared GO genes mean the seven terms are not seven independent programs. Some expanded annotation-signature markers still require provenance review; this limitation is documented and is not removed by a successful computational rerun.

## Results and analysis notes

- [Core 1 workflow README](README_core1.md)
- [Decision trail and complete numerical summary](docs/core1_analysis_summary.md)
- [Full rerun verification](docs/core1_rerun_verification.md)
- [Annotation review checklist](docs/core1_manual_review.md)
- [Formal DE summary](results/core1/tables/pseudobulk_DE_summary.csv)
- [P1 candidate contextualization](results/core1/tables/P1_P3_contextualization_summary.csv)
- [Complete GO × population matrix](results/core1/tables/P1_GO_population_matrix_complete.csv)

## Reproduction

Generated figures are distributed with the complete release artifact package.

The scripts are under `scripts/`. Four main stage files and a shared helper define functions without starting an analysis when sourced. `scripts/run_core1.R` provides the full workflow entry point.

1. Clone this repository and install the packages listed in the archived session information. Package versions from the successful rerun are preserved in `results/core1/diagnostics/`.
2. Download `core1-complete-artifacts.zip` from the [Core 1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1) and extract its contents into the repository root. This supplies the archived cell metadata, all-results tables, lightweight checkpoint and generated figures. Large artifacts are distributed with the release rather than repeated in the main source tree.
3. Obtain the original RNA files from [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495), using the exact filenames and paths in [data sources](docs/data_sources.md). The original matrices are not uploaded here.
4. From the project root, run:

```r
source("scripts/run_core1.R")
run_core1_full()
```

Or from a terminal with R installed:

```sh
Rscript scripts/run_core1.R
```

Each full run creates a new timestamped results directory. Clustering, strict markers, signature scores and DE are recomputed; reviewed manual labels and the resolved GO mapping are reused. A historical-cluster correspondence check stops annotation if the new partition or cluster IDs differ. Full Seurat checkpoints are generated locally and excluded from Git.

## Cross-project context

- [Project 1: Arabidopsis drought RNA-seq](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq)
- [Project 2: Arabidopsis root single-cell context](https://github.com/YuntongMeng/arabidopsis-root-singlecell)

Project 1 candidate and enriched GO inputs are fixed in `data/project1/`; Core 1 annotation and symbol-to-TAIR inputs are fixed in `data/core1_reference/`.
