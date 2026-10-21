# Arabidopsis root-tip multiome: osmotic-stress context for Project 1 genes

Project 1 found 137 genes whose drought response differs between
BRL3-overexpressing and wild-type roots. Here I asked whether those genes also
respond to osmotic stress in specific root-tip cell populations, and whether
nearby chromatin accessibility changes with them. I reanalyzed paired
single-nucleus RNA and ATAC data from control and sorbitol-treated wild-type
root tips (GEO GSE235495; Liu et al., 2024), with two biological replicates per
condition.

Related repositories: [Project 1](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq) · [Project 2](https://github.com/YuntongMeng/arabidopsis-root-singlecell). Part of [Arabidopsis root omics](https://github.com/YuntongMeng/arabidopsis-root-omics).

## Key findings

- **RNA (Core 1):** after removing ATAC peak rows that were mixed into the RNA
  input, 21,550 nuclei remained, and nine broad populations had enough nuclei
  for replicate-based pseudobulk DE. The Project 1 genes are enriched among
  osmotic-stress DEGs relative to the ordinary gene background in every
  eligible population, but not beyond Project 1's other drought-responsive
  genes.
- **ATAC (Core 2):** across 14,146 paired nuclei and five populations, one
  region passed within-population FDR < 0.05, in Atrichoblast_lineage near
  AT4G08570.
- **Linking the two:** 120 of the 137 Project 1 genes have a nearby common ATAC
  region, but none of those regions is significant. With two replicates per
  condition, this is weak evidence either way, not evidence that chromatin is
  unchanged.

## Analysis overview

The analysis has two parts:

| Part | What it does | Details |
| --- | --- | --- |
| Core 1 (snRNA-seq analysis) | Gene-only QC, predicted doublets, clustering, provisional annotation, pseudobulk DE and Project 1 context | [Core 1 page](cores/core1/README.md) |
| Core 2 (snATAC-seq and paired RNA–ATAC analysis) | Common-region ATAC recounting, QC, same-nucleus pairing with Core 1, pseudobulk differential accessibility and Project 1 neighbourhood context | [Core 2 page](cores/core2/README.md) |

My first version of Core 1 read the multiome matrices without separating gene
rows from peak rows, so ATAC peaks entered the RNA analysis. When I found this
I rebuilt Core 1 from gene-expression rows only; the first version is kept in
`archive/` and `results/core1_frozen/` for comparison. A planned Core 3 was
dropped because the paired RNA–ATAC comparison already belongs to Core 2. The
[project scope and history](docs/project_scope.md) has the dates and details.

## Inputs and provenance

- GEO GSE235495: four processed multiome matrices (control_rep1/2,
  stress_rep1/2). Only feature rows whose third column is exactly
  `Gene Expression` enter the RNA analysis; the `Peaks` rows are used for ATAC.
- ATAC fragments were reconstructed from the four SRA runs with Chromap and
  recounted on 27,677 common nuclear intervals.
- Project 1 inputs (137 interaction genes, seven GO terms, drought DEGs):
  `data/project1/` and `data/core1_corrected_reference/`.
- Annotation references (Shahan et al. 2022 Data S2 and extracted marker
  panels): `data/core1_corrected_reference/`. Core 2 reference intervals and
  reviewed cluster assignments: `data/core2_reference/`.

The GEO matrices, SRA reads, fragments and RDS checkpoints are not in Git.
[Data sources](docs/data_sources.md) explains how to restore them.

## Running the analysis

Both runners print their stages with `help` and run nothing until a stage is
named. From the project root:

```sh
Rscript scripts/install_core1_dependencies.R
Rscript scripts/run_core1_corrected.R help
Rscript scripts/run_core1_corrected.R all      # Core 1 stages 01–07 plus figures

Rscript scripts/run_core2.R help
python3 scripts/verify_core2.py                # check the published Core 2 record
python3 scripts/verify_core2_supplement.py     # check the supplementary analysis
```

Core 1 stages can also be run one at a time (`audit`, `qc`, `embeddings`,
`evidence`, `de`, `context`, `compare`, `figures`). Core 2 needs a local
upstream workspace; its stages and requirements are in
[Core 2 running instructions](docs/core2/running.md). `P3_ROOT`, `P3_R_LIB`
and `P3_CORE2_WORK` set the project root, an R library and the Core 2
workspace. Core packages include Seurat, Signac, scDblFinder, Harmony and
DESeq2; session versions are recorded under `results/core1/diagnostics/` and
`results/core2/diagnostics/`.

The reviewed annotation is tied to the MD5 of the saved partition, so a
changed clustering stops before DE until the annotation evidence is reviewed
again. Archived scripts should not be used to generate current results.

## Main analysis

### Core 1 (snRNA-seq analysis)

1. Audit feature types and keep gene-only RNA counts (32,833 genes per sample).
2. Set per-sample QC thresholds from the gene-only distributions and call
   doublets per capture with scDblFinder: 24,474 input nuclei, 430 removed by
   QC, 2,494 predicted doublets removed, 21,550 retained.
3. Compare unintegrated and Harmony embeddings. The unintegrated branch is
   primary and Harmony is a sensitivity branch.
4. Annotate broad populations from atlas signatures and canonical markers, sum
   raw counts by population × biological sample, and run DESeq2 for the nine
   eligible populations.
5. Test the 137 Project 1 genes for enrichment among stress DEGs, against both
   the ordinary gene background and Project 1's other drought DEGs, and
   summarize overlap with the seven Project 1 GO terms.

### Core 2 (snATAC-seq and paired RNA–ATAC analysis)

1. Merge the four published peak lists into 27,677 common intervals and
   recount reconstructed fragments.
2. Apply ATAC QC and keep nuclei that are also Core 1 singlets, matched by
   capture and barcode: 1,239 / 3,817 / 4,841 / 4,249 nuclei, 14,146 in total.
3. Cluster the ATAC data and describe each cluster by its RNA population
   composition.
4. Sum counts by RNA population × biological sample and test differential
   accessibility in five populations with DESeq2.
5. Map the 137 Project 1 genes to nearby common intervals with a fixed
   gene-body / TSS rule.

## Biological interpretation

The Project 1 genes do respond to osmotic stress in the root tip. In every
eligible population they are over-represented among stress DEGs compared with
ordinary genes; Pericycle, for example, has 19 of them among its DEGs in both
branches. But they are no more enriched than Project 1's other
drought-responsive genes. So this dataset links them to a general
water-stress response, not to the BRL3-dependent part of their behaviour.

Of the seven Project 1 GO terms, only jasmonic-acid metabolism, long-chain
fatty-acid metabolism and response to cold have any overlap with stress DEGs
here, and these terms share genes.

Chromatin accessibility changed little at the thresholds used. The one
significant region, in Atrichoblast_lineage, is higher in both stress
replicates than in both controls. It lies near AT4G08570, but proximity alone
does not make it a regulatory region for that gene.

## Robustness checks

- **Core 1:** QC sensitivity, PC/ARI stability, replicate mixing, canonical
  markers, condition preservation, population eligibility at 25/50/100 nuclei
  and DE concordance between the two branches. No cluster count, population
  count or DEG total is treated as a required answer.
- **Core 2:** a broader count filter (counts ≥ 5) beside the primary counts ≥ 10
  run, model diagnostics, an independent recount of Control 1's
  low-signal subgroup, and a read-only audit of 111 count, model and mapping
  checks. Control 1's subgroup is already present in the published matrix; its
  experimental cause is unknown.
- **Reruns:** I reran Core 1 stages 01–07 from the original inputs (all eight
  structural checks passed) and the assembled Core 2 downstream workflow from
  the saved QC objects (all ten comparisons matched). FASTQ processing was not
  repeated.

Details: [Core 1 decision log](docs/core1/decision_log.md),
[Core 2 validation](docs/core2/validation.md) and
[Core 2 supplementary review](docs/core2/supplemental_review.md).

## Interpretation boundary

There are two biological replicates per condition, and nuclei are not
replicates. Population labels are provisional; some are broad lineages or
mixed groups, and Root_cap is exploratory in Core 2. Predicted doublets are not
ground truth, and ambient RNA was not corrected because no empty-droplet matrix
is available. Harmony may weaken condition structure, and these comparisons
cannot separate every technical and biological contribution. FDR in Core 2 is
corrected within each population, not across all five.

This is wild-type data under sorbitol, not BRL3 plants under drought. It can
show whether the Project 1 genes respond to a related stress in particular
root-tip populations. It cannot test the BRL3 × drought interaction itself.
Project 2 adds where these genes are normally expressed; this project adds how
they respond to stress, but neither establishes a cell-specific BRL3 effect.

## Outputs and repository structure

- Core 1: [tables](results/core1/tables/) (large tables as `.csv.gz`),
  [figures](results/core1/figures/), [diagnostics](results/core1/diagnostics/).
- Core 2: [tables](results/core2/tables/), [figures](results/core2/figures/),
  [diagnostics and provenance](results/core2/diagnostics/),
  [supplement](results/core2/supplement/).
- Core 1 reasoning: [analysis story](docs/core1/analysis_story.md) (start
  here), [decision log](docs/core1/decision_log.md),
  [version comparison](docs/core1/version_comparison.md),
  [detailed comparison](docs/core1/comparison_report.md),
  [manual review order](docs/core1/manual_review.md),
  [manual rerun](docs/core1/manual_rerun_validation.md),
  [figure guide](docs/core1/figure_review.md).
- Core 2 reasoning: [conclusions](docs/core2/analysis_story.md),
  [methods](docs/core2/methods.md), [decision log](docs/core2/decision_log.md),
  [Control 1 review](docs/core2/control1_review.md),
  [directory guide](docs/core2/directory_map.md).

```text
scripts/core1/                 Core 1 stage functions
scripts/run_core1_corrected.R  Core 1 stage runner
scripts/core2/                 Core 2 stages and recorded upstream reconstruction
scripts/run_core2.R            Core 2 stage runner
results/core1/                 Core 1 tables, figures and diagnostics (checkpoints stay local)
results/core2/                 Core 2 tables, figures, diagnostics and supplement
docs/core1/, docs/core2/       Reasoning, methods, decisions and validation
data/                          Inputs, references and provenance
exploration/core1/             Extra diagnostics and reference-extraction tools
archive/core1_initial/         First Core 1 version (mixed input), for the record
results/core1_frozen/          Frozen first-version results, for comparison only
```

## Data and attribution

This is an independent computational reanalysis of public data. The
experimental design, plant material, sequencing and the conclusions of the
original publication belong to its authors (Liu et al., 2024). The scripts, statistics, figures
and interpretations here come from this reanalysis.

I did this as a bioinformatics training and reproducibility exercise, not as an
experimental replication of the original study.

## References

- [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495).
- [Liu et al. (2024), original root-tip multiome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11199978/).
- [Shahan et al. (2022), Arabidopsis root atlas](https://doi.org/10.1016/j.devcel.2022.01.008).
- [scDblFinder method and assumptions](https://plger.github.io/scDblFinder/articles/scDblFinder.html).
