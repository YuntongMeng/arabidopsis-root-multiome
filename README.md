# Arabidopsis root-tip multiome under osmotic stress

An independent reanalysis of publicly available *Arabidopsis thaliana* root-tip
multiome data to investigate cell-population-specific responses to osmotic stress.

Project 3 extends the biological questions developed in
[Project 1](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq) and
[Project 2](https://github.com/YuntongMeng/arabidopsis-root-singlecell).
Project 1 identified BRL3 genotype × drought interaction genes in bulk roots;
Project 2 placed those candidates into a normal wild-type root atlas. Project 3
adds a wild-type root-tip osmotic-stress context using GSE235495.

## Key biological findings

The findings currently available come from **Core 1 — RNA analysis**:

- QC retained **24,277 nuclei**, yielding **18 transcriptional clusters** and
  **17 conservatively annotated populations**.
- **13 populations** met the requirement of at least 25 nuclei in each biological
  sample and entered formal stress-versus-control pseudobulk differential analysis.
- Of the 137 Project 1 interaction candidates, **129 occurred in at least one
  population's differential-expression result table**. Pericycle had the largest
  significant-response fraction: **23/118 candidates (19.49%)**, with 6 up and
  17 down under osmotic stress.
- Pericycle contributed downregulated candidates to all seven selected Project 1
  GO annotations; the cold-response annotation contained three such genes.

These results provide cellular context for the Project 1 candidates. They do
not validate BRL3-dependent drought effects: Project 3 uses wild-type plants,
root tips, and sorbitol osmotic stress. The GO terms share genes and should not
be treated as seven independent biological programs.

## Project 3 analysis overview

Project 3 is organized into three cores. **Only Core 1 has been conducted.**

| Core | Analysis | Status |
| --- | --- | --- |
| [Core 1](cores/core1/README.md) | RNA preprocessing, clustering, population annotation, pseudobulk stress DE, and Project 1 contextualization | Complete; full rerun verified against key archived outputs |
| [Core 2](cores/core2/README.md) | Reserved for the next planned core | Planned; not started |
| [Core 3](cores/core3/README.md) | Reserved for the final planned core | Planned; not started |

### Core 1 — RNA analysis

Core 1 asks which root-tip populations respond transcriptionally to osmotic
stress and where Project 1 interaction candidates occur within those responses.
Its four stages preserve the reasoning behind QC, PC selection, conservative
annotation, replicate eligibility, and differential-expression testing.

The [Core 1 analysis page](cores/core1/README.md) contains the script overview,
complete population DE summary, biological interpretation, robustness checks,
and reproduction instructions. The
[decision trail](docs/core1_analysis_summary.md) records intermediate method
adjustments and distinguishes formal results from diagnostics.

### Core 2 — Planned

This core has not been started. Its [reserved page](cores/core2/README.md) will
hold the research question, workflow, and results when the analysis is conducted.

### Core 3 — Planned

This core has not been started. Its [reserved page](cores/core3/README.md) will
hold the research question, workflow, and results when the analysis is conducted.

## Inputs and provenance

Core 1 uses four RNA count matrices from
[GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495):
two control and two sorbitol-treated wild-type root-tip biological samples.
The sample mapping and input filenames are documented in
[data sources](docs/data_sources.md) and
[sample metadata](data/rna/sample_metadata.tsv).

Project 1 supplies the fixed 137-gene interaction table and seven
positive-interaction GO enrichment terms in [data/project1](data/project1).
Annotation signatures, reviewed labels, and resolved symbol-to-TAIR mappings
are retained in [data/core1_reference](data/core1_reference).
Some expanded signature genes still require provenance review; the uncertainty
remains visible in the [annotation review notes](docs/core1_manual_review.md).

## Reproducing the analysis

The complete Core 1 workflow can be run from the repository root after restoring
the original matrices and the required archived inputs:

```r
source("scripts/run_core1.R")
run_core1_full()
```

Detailed input preparation, package records, and output instructions are on the
[Core 1 page](cores/core1/README.md#running-the-complete-workflow).
The [Core 1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1)
contains complete gene-level results, generated figures, diagnostics, and a
lightweight archived checkpoint. Raw matrices and full Seurat objects are
excluded. Core 2 and Core 3 have no executable workflows or results yet.

## Outputs and repository structure

```text
arabidopsis-root-multiome/
├── README.md                 # Project 3 overview
├── cores/
│   ├── core1/README.md       # Completed RNA analysis
│   ├── core2/README.md       # Planned; not started
│   └── core3/README.md       # Planned; not started
├── scripts/                 # Core 1 executable stages and full runner
├── data/                    # Sample metadata and fixed reference inputs
├── docs/                    # Decision trail, provenance and verification
└── results/
    ├── core1/               # Successful rerun tables and diagnostics
    └── core1_frozen/        # Archived original analysis
```

Scripts and inputs retain their project-root paths so the verified Core 1
workflow can be run with the same file layout. Figures and large result files
are distributed in the release package.

## Interpretation boundary

Core 1 measures WT osmotic-stress responses at biological-sample resolution.
Population labels are conservative summaries of signature and marker evidence;
low-confidence populations retain descriptive labels. A successful rerun
supports computational reproducibility of the checked outputs, while signature
provenance and biological interpretation still require their own evidence.

## Data and attribution statement

This repository contains an independent computational reanalysis of public
data. Experimental design, plant material, sample preparation, sequencing, and
primary data generation belong to the original study's authors. The scripts,
processing decisions, statistical analyses, and interpretations here describe
this reanalysis.

## References and related projects

- [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495): public input dataset for Project 3.
- [Project 1 — Arabidopsis drought RNA-seq reanalysis](https://github.com/YuntongMeng/arabidopsis-drought-rnaseq): BRL3 genotype × drought interaction candidates and GO inputs.
- [Project 2 — Arabidopsis root single-cell contextualization](https://github.com/YuntongMeng/arabidopsis-root-singlecell): normal WT root cellular context.
