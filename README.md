# Arabidopsis root-tip multiome stress contextualization

Project 3 reanalyzes wild-type Arabidopsis root-tip multiome RNA and ATAC data to ask how
Project 1 BRL3 genotype × drought interaction candidates relate to population-level
stress responses. Project 2 supplies normal-root cellular context; Project 3 adds
condition-resolved context. These datasets do not establish a cell-specific BRL3 interaction.



Related repositories: [Project 1](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq) · [Project 2](https://github.com/YuntongMeng/arabidopsis-root-singlecell).

## Analysis status

**Project complete within the agreed scope as of 2026-10-15: Core 1 + Core 2.**
The planned Core 3 has been cancelled, not left pending. The same-nucleus
RNA–ATAC comparison is already included in Core 2. See the
[final scope decision](docs/project_scope.md).

Core 1 was rerun on 2026-10-06 after discovering that the original RNA input
included ATAC peaks. The current analysis uses gene-expression rows only and
retains 21,550 predicted singlets. Unintegrated clustering is primary; sample-based
Harmony is retained as a sensitivity comparison. Annotation remains provisional.
The user manually reran stages 01–07; all eight structural checks passed.
Core 2 now includes common-region ATAC recounting, QC and exact-nucleus RNA context
for 14,146 paired nuclei. Five population-level pseudobulk comparisons identify one
candidate peak at within-population FDR <0.05; none of the 137 P1 candidates has
significant matched nearby ATAC support under the fixed mapping rule. Results are
limited and labels remain provisional. The project concludes at Core 2;
there is no remaining Core 3 deliverable.

The current version is a documented working analysis, not a claim of a uniquely
correct solution. Earlier assumptions, errors, alternatives and remaining uncertainty
are recorded separately from the executable current workflow.

## Project cores

| Core | Status | Overview |
| --- | --- | --- |
| Core 1 | Corrected RNA workflow; manual rerun checked | [RNA analysis](cores/core1/README.md) |
| Core 2 | ATAC / paired RNA analysis complete; assembled downstream workflow manually rerun and matched | [ATAC and paired context](cores/core2/README.md) |

## Inputs and provenance

- GEO GSE235495: four processed multiome feature matrices, two samples per condition.
- Only feature rows whose third column is exactly `Gene Expression` enter RNA analysis.
- Raw matrices and sample metadata: `data/rna/`.
- Project 1 contextualization inputs: `data/project1/`.
- Reviewed annotation, Shahan 2022 Data S2 and extracted panels: `data/core1_corrected_reference/`.
- Other existing references remain under `data/reference/` and `data/core1_reference/`.

No verified unfiltered empty-droplet matrix is available locally. Ambient RNA
was assessed qualitatively; correction was not forced.

## Running Core 1

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

## Core 2 workflow and current findings

[Start with Core 2](cores/core2/README.md) for the conclusions, figures and interpretation boundaries.
Four peak sets were union/reduced and recounted from reconstructed fragments into
27,677 common nuclear intervals. ATAC QC and corrected RNA singlet matching retain
1,239 / 3,817 / 4,841 / 4,249 nuclei across C1/C2/S1/S2. Shared capture + barcode
identifiers preserve exact same-nucleus pairing; Gene Expression and peak features
remain separate assays.

The single differential candidate is chr4:5,457,000–5,457,860 in Atrichoblast_lineage
(estimated log2FC 1.860; population-wise adjusted p 0.004526). It is near AT4G08570,
not a validated regulatory link. All other tested populations have no significant
peaks at the current criterion, with Root_cap explicitly exploratory. Of 137 P1
candidates, 120 match common regions; none has a matched FDR-significant region.
Lack of significance does not establish absence of stress effects.

Control 1 contains a pre-existing low-nuclear-signal, organelle-rich subgroup.
Independent full-fragment recounting and matrix comparison found no error in the
checked quantities; the specific experimental cause remains unresolved.
[Evidence and audit limits](docs/core2/control1_review.md).

```sh
Rscript scripts/run_core2.R help
python3 scripts/verify_core2.py
```

These inspect stages and validate the published record without rerunning biological
analysis. [Explicit execution instructions](docs/core2/running.md) require a local
upstream workspace. The previously interactive stages have now been assembled into
scripts and manually rerun by the user from the four existing unfiltered QC objects.
All ten recorded comparisons matched, including raw pseudobulk counts, reviewed
ATAC assignments and five population differential statistics (tolerance 1e-7).
FASTQ reconstruction and QC metric calculation were not repeated. Saved-coordinate
PCA figures were redrawn with equal x/y unit scale; their independent panels are
not directly comparable coordinate spaces.

No new ATAC GO analysis is reported. The project ends at Core 2; Core 3 is cancelled. [Submission status](docs/core2/review_before_push.md). See the [directory guide](docs/core2/directory_map.md) for current scripts, reviewed results and local-only data.

## Core 1 main analysis

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

## Core 2 outputs

- [Result tables](results/core2/tables/), [main figures](results/core2/figures/) and [diagnostics/provenance](results/core2/diagnostics/).
- [Conclusions](docs/core2/analysis_story.md), [methods](docs/core2/methods.md), [decision trail](docs/core2/decision_log.md) and [validation limits](docs/core2/validation.md).
- Local fragments, reads, runtimes and RDS checkpoints are excluded from Git. Large CSV results are losslessly compressed.

## Repository structure

```text
scripts/core2/          Core 2 explicit stages and recorded upstream reconstruction
scripts/run_core2.R    Explicit Core 2 runner; help performs no analysis
results/core2/         ATAC tables, figures, diagnostic/provenance records
docs/core2/            ATAC conclusions, methods, audit, running and review
data/core2_reference/  Common BED, reviewed cluster binding and provenance
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
- [Liu et al. (2024), original root-tip multiome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11199978/).
- [scDblFinder method and assumptions](https://plger.github.io/scDblFinder/articles/scDblFinder.html).
