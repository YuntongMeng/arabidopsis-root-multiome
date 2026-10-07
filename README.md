# Arabidopsis root-tip multiome stress contextualization

Project 3 reanalyzes wild-type Arabidopsis root-tip multiome RNA data to ask how
Project 1 BRL3 genotype × drought interaction candidates relate to population-level
stress responses. Project 2 supplies normal-root cellular context; Project 3 adds
condition-resolved context. These datasets do not establish a cell-specific BRL3 interaction.



Related repositories: [Project 1](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq) · [Project 2](https://github.com/YuntongMeng/arabidopsis-root-singlecell).

## Analysis status

Core 1 was rerun on 2026-10-06 after discovering that the original RNA input
included ATAC peaks. The current analysis uses gene-expression rows only and
retains 21,550 predicted singlets. Unintegrated clustering is primary; sample-based
Harmony is retained as a sensitivity comparison. Annotation remains provisional.
The user manually reran stages 01–07; all eight structural checks passed.
Core 2 and Core 3 are planned and have not started. Biological labels remain provisional.

The current version is a documented working analysis, not a claim of a uniquely
correct solution. Earlier assumptions, errors, alternatives and remaining uncertainty
are recorded separately from the executable current workflow.

## Project cores

| Core | Status | Overview |
| --- | --- | --- |
| Core 1 | Corrected RNA workflow; manual rerun checked | [RNA analysis](cores/core1/README.md) |
| Core 2 | Planned; not started | [Core 2](cores/core2/README.md) |
| Core 3 | Planned; not started | [Core 3](cores/core3/README.md) |

## Inputs and provenance

- GEO GSE235495: four processed multiome feature matrices, two samples per condition.
- Only feature rows whose third column is exactly `Gene Expression` enter RNA analysis.
- Raw matrices and sample metadata: `data/rna/`.
- Project 1 contextualization inputs: `data/project1/`.
- Reviewed annotation, Shahan 2022 Data S2 and extracted panels: `data/core1_corrected_reference/`.
- Other existing references remain under `data/reference/` and `data/core1_reference/`.

No verified unfiltered empty-droplet matrix is available locally. Ambient RNA
was assessed qualitatively; correction was not forced.

## Running the complete workflow

From the project root:

```bash
Rscript scripts/run_core1_corrected.R help
```

This displays the available stages without rerunning biological analysis.
Run a stage explicitly: `audit`, `qc`, `embeddings`, `evidence`, `de`, `context`,
`compare`, `diagnostics`, or `figures`. `all` explicitly recomputes stages 01–07
and renders figures. Historical `diagnostics` are optional and excluded from `all`.
The entry filename is retained for continuity; its implementation now lives in `scripts/core1/`.

Set `P3_ROOT` when running outside the project root. `P3_R_LIB` can specify a local
R package library. The default root is the current working directory and standard R package libraries
are used. Optional machine settings belong in an untracked `.Renviron`; personal
absolute paths are not part of the executable defaults. Install missing dependencies
with `Rscript scripts/install_core1_dependencies.R`. This is not a version lock;
reviewed session versions are recorded alongside the results.
Core packages include Seurat, Matrix, scDblFinder, Harmony and DESeq2. Stage-specific
versions are recorded in `results/core1/diagnostics/`.

Sourcing the stage files defines functions. The reviewed annotation is bound to
the saved partition MD5: changed clustering requires another evidence review before DE.
Do not rerun archived scripts to generate current results.

## Main analysis

1. Audit feature types and retain gene-only RNA counts.
2. Inspect gene-only QC distributions, select sample-specific thresholds and call doublets per capture.
3. Compare unintegrated and Harmony embeddings, PC stability and annotation evidence.
4. Aggregate raw counts by population × biological sample and run DESeq2 for eligible populations.
5. Contextualize the 137 Project 1 candidates and seven existing GO terms.

The main analysis uses the unintegrated branch. Harmony results and all diagnostic
comparisons remain available; their existence does not imply that integration is mandatory.

## Robustness checks

QC sensitivity, PC/ARI comparisons, replicate mixing, canonical markers, condition
preservation, population eligibility and branch-level DE concordance are retained.
Some sensitivity checks change only retention or eligibility, not the whole downstream workflow.
No cluster count, population count or DEG total is asserted as a required answer.

## Outputs

- [Current result tables](results/core1/tables/) (large tables published as `.csv.gz`).
- [Current figures](results/core1/figures/).
- Saved checkpoints are generated locally and excluded from Git.
- [Environment and diagnostic records](results/core1/diagnostics/).
- [Analysis story](docs/core1/analysis_story.md): start here to follow our reasoning.
- [Decision log](docs/core1/decision_log.md): parameters, evidence and uncertainty.
- [Version comparison](docs/core1/version_comparison.md): what changed and what remains supported.
- [Detailed comparison report](docs/core1/comparison_report.md): full numerical review.
- [Manual review checklist](docs/core1/manual_review.md).
- [Manual rerun validation](docs/core1/manual_rerun_validation.md).
- [Figure review and rendering guide](docs/core1/figure_review.md): 22 main figures and two supplementary ARI diagnostics.

## Repository structure

```text
scripts/core1/          Current Core 1 workflow functions
scripts/run_core1_corrected.R  Explicit stage runner
results/core1/          Current tables, figures, checkpoints and diagnostics
docs/core1/             Reasoning, decisions, comparison and review guide
exploration/core1/      Additional diagnostics and reference-extraction tools
archive/core1_initial/  Earlier scripts, documents and intermediate result versions
results/core1_frozen/   Unchanged original frozen analysis, for comparison only
data/                  Original inputs and reference provenance
```

The historical archive README explains validity and location changes. Historical
text retains its original assumptions and paths; those paths are not current run instructions.
The frozen analysis stays at its original location for traceability.

See [data sources and file layout](docs/data_sources.md) for restoring the GEO matrices.

## Interpretation boundary

Two biological replicates per condition limit inference. Broad lineage labels
and unresolved states must be distinguished. Predicted doublets are not a ground-truth
classification; ambient contamination is not ruled out. Sample-based Harmony may
attenuate condition geometry, but these comparisons cannot disentangle every
technical and biological contribution.

P1 candidates show stress-DEG enrichment relative to the ordinary gene background,
but no population passes BH correction for additional enrichment relative to the
P1 noninteraction drought-DEG control. The seven GO views summarize overlapping
candidate gene sets; they are not new P3 GO enrichment tests.

## Data and attribution statement

This is an independent computational reanalysis for bioinformatics training and
reproducibility. Experimental design, plant material and primary data generation
belong to the original studies. Cross-project contextualization is not an independent
experimental replication or mechanistic validation.

## References

- [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495).
- [Shahan et al. (2022), Arabidopsis root atlas](https://doi.org/10.1016/j.devcel.2022.01.008).
- [scDblFinder method and assumptions](https://plger.github.io/scDblFinder/articles/scDblFinder.html).
