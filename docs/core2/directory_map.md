# Core 2 directory guide

## Repository: current code and reviewed outputs

| Location | Purpose |
| --- | --- |
| [Core 2 overview](../../cores/core2/README.md) | Start here for conclusions and selected figures |
| [scripts/core2](../../scripts/core2/) | Numbered downstream functions; Source defines them only |
| [RStudio validation entry](../../scripts/core2/validate_in_rstudio.R) | Explicit fresh downstream run into a separate directory |
| [scripts/core2/upstream](../../scripts/core2/upstream/) | Recorded upstream recovery/QC code; inspect before explicit execution |
| [docs/core2](.) | Conclusions, methods, decisions, Control 1 review, validation and running guide |
| [results/core2/figures](../../results/core2/figures/) | Selected figures; corrected PCA uses equal x/y units |
| [results/core2/tables](../../results/core2/tables/) | Reviewed numerical results, including full compressed DA tables |
| [results/core2/diagnostics](../../results/core2/diagnostics/) | QC plots, audits, session information and user rerun evidence |
| [data/core2_reference](../../data/core2_reference/) | Common intervals, reviewed cluster binding and source metadata |

Core 1 analysis/results and `results/core1_frozen` are preserved. Core 3 was
cancelled on 2026-10-15; the project is complete at Core 2. Same-nucleus RNA
identity mapping is already included in Core 2. See the [scope decision](../project_scope.md).

## Local reconstruction workspace: preserve file locations

`P3_CORE2_WORK` points to a separate workspace. Its `outputs/` contains four
`*_upstream/` directories, `core2_qc/`, the reviewed `core2_analysis/`, and the
user's fresh `core2_validation_rstudio/`. Its `work/` contains FASTQ/SRA staging,
tools, reference annotation and the R library. These large local files are not
published in Git. RDS fragment references and upstream scripts depend on their
locations, so directory organization preserves them in place.

Root-level teaching scripts are historical working copies. Use the repository
scripts for future work; do not assume a teaching copy has the same Source
behavior. Old preparation notes/packages are grouped under local
`outputs/archive/2026-10-13_prior_preparation/` and marked superseded.

The repository is the current code/results record. The local review mirror and
ZIP are snapshots of the Core 2 additions/updates, not full repository backups.
The local sequencing workspace remains necessary for reconstruction.

## Local supplementary revision

`results/core2/supplement/` contains the ten full-size diagnostic figures and numeric summaries, separate count-filter sensitivity results, exploratory RNA/ATAC context, input-integrity checks and execution record. It is outside the main figure gallery. Large fitted RDS models remain in the external workspace `outputs/core2_supplement_2026-10-15/`. The private discussion PDF is temporary and is not included in the project, result snapshot, or Git. The user authorized publication of these revisions on 2026-10-18; the private review PDF remains excluded.
