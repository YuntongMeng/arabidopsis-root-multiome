# Project 3 — Core 1: RNA analysis under osmotic stress

[Project 3 overview](../../README.md) · [Core 2 — planned](../core2/README.md) · [Core 3 — planned](../core3/README.md)

Core 1 investigates cell-population-specific transcriptional responses in
wild-type Arabidopsis root tips under sorbitol osmotic stress. It also asks where
Project 1 BRL3 genotype × drought interaction candidates occur among those
responses. The full RNA workflow has been rerun from the original matrices;
its checked key outputs agree with the archived analysis.

## Key biological findings

- **24,277 nuclei**, **18 clusters**, and **17 annotated populations** were retained.
- **13 populations** had at least 25 nuclei in every biological sample and were
  eligible for formal pseudobulk differential-expression analysis.
- **129 of 137 Project 1 candidates** occurred in at least one population result
  table. Pericycle had the largest significant-response fraction,
  **23/118 (19.49%)**, comprising 6 upregulated and 17 downregulated candidates.
- Pericycle appeared in all seven selected Project 1 GO annotations. Its three
  cold-response overlap genes were all downregulated.

These are osmotic-stress cellular contexts, not validation of the BRL3 × drought
interaction. The seven GO annotations overlap in gene membership.

## Analysis overview

| Script | Purpose | Main checked result |
| --- | --- | --- |
| [01 — RNA preprocessing and clustering](../../scripts/01_core1_rna_preprocessing_clustering.R) | Audit matrices, inspect QC distributions, compare PC choices, and inspect sample composition | 24,277 nuclei; PC1:25; 18 clusters |
| [02 — Population annotation](../../scripts/02_core1_cell_population_annotation.R) | Join RNA layers, calculate strict markers and signature scores, and apply conservative reviewed labels | 17 populations; cluster 12 relaxed markers retained as diagnostics only |
| [03 — Pseudobulk stress DE](../../scripts/03_core1_pseudobulk_stress_DE.R) | Aggregate by population × biological sample, enforce eligibility, and fit DESeq2 models | 13 populations; 32,833 shared input features |
| [04 — Project 1 contextualization](../../scripts/04_core1_P1_contextualization.R) | Summarize the 137 candidates and seven existing GO annotations in P3 | Pericycle 23/118 candidates; complete 7 × 13 GO matrix |

The [shared helper](../../scripts/00_core1_common.R) supports file handling;
[run_core1.R](../../scripts/run_core1.R) runs the complete workflow. Sourcing
these files defines functions. Calling `run_core1_full()` starts the analysis.
Ambiguous non-base functions use explicit package namespaces.

## Inputs and provenance

The four RNA matrices come from GSE235495: control_rep1, control_rep2,
stress_rep1, and stress_rep2. Each condition has two biological samples.
Exact GEO sample IDs, matrix filenames, and required directories are listed in
[data sources](../../docs/data_sources.md).

Fixed inputs include the Project 1 interaction and GO tables,
reviewed population labels, annotation signature versions, and resolved GO
symbol-to-TAIR mappings. Expanded signature-marker provenance remains partly
unresolved; see the [review checklist](../../docs/core1_manual_review.md).

## Main analysis

### QC, clustering, and annotation

Matrix auditing precedes filtering. QC uses the inspected count and feature
distributions: `800 <= nFeature_RNA <= 15000` and `nCount_RNA <= 50000`.
Mitochondrial and chloroplast percentages remain diagnostic metrics without a
mechanical organelle cutoff. LogNormalize (scale factor 10,000), 2,000 vst
variable features, scaling, and PCA precede clustering at resolution 0.5.

PC20 was an initial candidate from the elbow plot. PC10/15/20/25/30 partitions
were compared using pairwise ARI before selecting PC1:25. Harmony was not
applied automatically because condition-associated variation is part of the
research question. Clusters 6 and 14 were checked for replicate enrichment.

Annotation uses reference/canonical signature scores with manual interpretation.
Dataset-specific strict markers and a canonical DotPlot provide supporting
checks. Low-confidence groups retain conservative labels. Seurat v5 RNA layers
are joined before formal marker testing; cluster 12 has no markers under the
same strict criteria used for other clusters. Its relaxed marker test is a
separate diagnostic.

### Biological-sample pseudobulk differential expression

Raw counts are summed by population × biological sample. Formal analysis
requires at least 25 nuclei in each of the four samples. The feature universe
is restricted to the **32,833 features shared across all four original samples**,
avoiding structural zeros introduced by merged union features. Genes with
summed counts below 10 in a population are removed before model fitting.

DESeq2 uses `~ condition` with the stress/control contrast. The two replicate
numbers are not treated as paired blocks. Significant genes meet
`padj < 0.05` and `abs(log2FoldChange) >= 1`.

| Population | DEGs | Up | Down |
| --- | ---: | ---: | ---: |
| Atrichoblast | 687 | 488 | 199 |
| Columella | 151 | 55 | 96 |
| Cortex | 259 | 86 | 173 |
| Dividing_S_phase | 214 | 115 | 99 |
| Dividing_meristematic | 410 | 323 | 87 |
| Endodermis | 104 | 16 | 88 |
| Endodermis_like | 504 | 82 | 422 |
| Epidermal_like | 535 | 275 | 260 |
| Initial_like | 443 | 284 | 159 |
| LRC_like | 524 | 155 | 369 |
| Pericycle | 663 | 136 | 527 |
| Trichoblast | 739 | 251 | 488 |
| Unresolved_transition | 96 | 50 | 46 |

The [formal summary table](../../results/core1/tables/pseudobulk_DE_summary.csv)
records the number of genes in each result table as well as DEG counts.

### Project 1 candidate and GO context

The formal candidate input is the existing 137-gene BRL3 × drought interaction
set. Population-specific denominators count candidates present in that
population's filtered DE result table, including rows with NA adjusted p-values.
The 129-gene total is a distinct count across populations.

Seven existing Project 1 enriched GO terms are mapped to TAIR genes and joined
to P1 candidates that are significant P3 stress DEGs. The complete matrix
includes all 13 eligible populations, including zero-overlap combinations.
This is descriptive contextualization; no new GO enrichment test is performed.

## Biological interpretation

Pericycle provides a prominent context for the Project 1 candidate set under
WT osmotic stress. Its candidate response fraction is the largest among the
eligible populations, and all seven selected GO annotations have downregulated
overlap genes there. This prioritizes a cellular setting for follow-up; it does
not demonstrate that a BRL3-dependent drought mechanism occurs in Pericycle.

Responses also occur in other root-tip populations. DEG counts describe genes
passing the selected thresholds and depend on expression, effect size,
dispersion, and sample representation. They alone do not rank the biological
strength of stress responses.

## Robustness checks and decision trail

The scripts preserve the PC comparison, cluster 6/14 replicate checks, RNA-layer
repair, strict-versus-relaxed marker distinction, pseudobulk name-matching
repair, and shared-feature correction. Full details and all population-level
P1 results are in the [analysis decision trail](../../docs/core1_analysis_summary.md).

The [full rerun verification](../../docs/core1_rerun_verification.md) checks retained
nuclei, the PC ARI matrix, population DE summaries, P1 summaries, GO matrix,
and final labels against the archive. It does not claim a row-by-row comparison
of every gene's fitted statistics. Computational agreement does not resolve
annotation-marker provenance.

## Running the complete workflow

1. Restore the original RNA matrices using the exact paths in
   [data sources](../../docs/data_sources.md).
2. Download `core1-complete-artifacts.zip` from the
   [Core 1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1).
   Restore its `data/` and `results/` contents into the repository root. The
   release is the original verified snapshot; keep the current repository
   scripts and documentation when extracting it.
3. Install the R packages recorded in the successful run's
   [diagnostics](../../results/core1/diagnostics).
4. Start R in the project root and run:

```r
source("scripts/run_core1.R")
run_core1_full()
```

Alternatively:

```sh
Rscript scripts/run_core1.R
```

Each full run creates a timestamped results directory. It recomputes clustering,
strict markers, signature scores, and DE while reusing reviewed labels and the
resolved GO mapping. A historical-partition correspondence check stops
annotation if clusters or their numbering change. Full Seurat checkpoints are
created locally and excluded from Git.

## Outputs

- [results/core1/tables](../../results/core1/tables): formal DE, population and
  Project 1 contextualization summaries.
- [results/core1/diagnostics](../../results/core1/diagnostics): sample-composition
  checks, eligibility, shared features, package and session records.
- [results/core1_frozen](../../results/core1_frozen): archived original results.
- [Core 1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1):
  complete gene-level tables, all generated figures, archived cell metadata,
  and the lightweight checkpoint.

## Interpretation boundary

Project 3 Core 1 uses WT root tips under sorbitol osmotic stress; Project 1
uses a BRL3 × drought design in bulk roots. Candidate and GO overlap therefore
adds context rather than validation. Zero overlap does not establish absence
of a biological program. Uncertain signature provenance and low-confidence
labels remain documented.
