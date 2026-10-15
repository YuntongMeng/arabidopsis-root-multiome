# Core 2: inspection and explicit execution

## Read the results without rerunning analysis

Start at [Core 2](../../cores/core2/README.md), [conclusions](analysis_story.md), [tables](../../results/core2/tables/) and [figures](../../results/core2/figures/). Large published CSVs are losslessly gzip-compressed. R reads them with `utils::read.csv("path/to/table.csv.gz")`; Python can use `gzip.open(path, "rt")`.

```sh
Rscript scripts/run_core2.R help
python3 scripts/verify_core2.py
```

Help lists stages. Verification checks the published recorded results and package, without fitting models or scanning fragments. Neither command downloads or reruns biological analysis.

## Configuration

Start in the repository root, or set `P3_ROOT` to it. Set `P3_CORE2_WORK` to the local reconstruction workspace. That external directory has the following layout:

```text
work/annotation/Arabidopsis_thaliana.TAIR10.62.gtf.gz
work/barcode_reference/atac_737K-arc-v1.txt.gz
work/barcode_reference/gex_737K-arc-v1.txt.gz
work/barcode_reference/sources.json
work/atac_env/bin/                         chromap, samtools, bgzip, tabix
work/tools/sratoolkit.3.4.1-mac-arm64/bin/   original SRA Toolkit location
work/raw_sra/<SRR>/<SRR>.sra               remaining samples
work/remaining_run_metadata.tsv           public run metadata
outputs/core2_common_peaks.rds
outputs/core2_common_peaks.bed
outputs/<control1|control2|stress1|stress2>_upstream/
outputs/core2_qc/<sample>/atac_qc_unfiltered.rds
outputs/core2_analysis/                   original reviewed checkpoints/models
outputs/core2_rebuild/                    default new downstream-run destination
```

Control 1's downloaded SRA was originally placed in `P3_ROOT/data/atac/raw/SRR24921616/`; its preview directory is beside the SRA. These local raw data are excluded from Git. Restore the original processed GEO inputs and corrected Core 1 output tables before reconstruction. The candidate list and annotation provenance are already in the project.

Optional `P3_R_LIB` selects a local R library; otherwise the existing workspace `work/R_library` is used when present. `P3_CORE2_ANALYSIS` selects an analysis directory. By default it is **workspace/outputs/core2_rebuild**, keeping the original reviewed analysis separate. The scripts do not install packages automatically. `Rscript scripts/install_core2_dependencies.R help` describes an explicit missing-package installer; use its `install` command only when needed. Recorded R sessions are in results/core2/diagnostics; this is not an environment lock.

## Reconstructed downstream steps

The user ran the original analysis in RStudio. The previously interactive steps are now organized into callable files. Source defines their functions; the runner invokes only the selected stage. Run one stage at a time and inspect it:

```sh
Rscript scripts/run_core2.R select
Rscript scripts/run_core2.R embeddings
Rscript scripts/run_core2.R annotation
Rscript scripts/run_core2.R markers       # optional explanatory plots
Rscript scripts/run_core2.R pseudobulk
Rscript scripts/run_core2.R diagnostics
Rscript scripts/run_core2.R da
Rscript scripts/run_core2.R p1
Rscript scripts/run_core2.R pca
```

There is intentionally no implicit `all` command. These stages compute or overwrite outputs in the explicitly selected analysis directory. Do not point a fresh reconstruction at the reviewed directory unless intentionally replacing it. `diagnostics` can reuse a fitted model only after counts/design/sample validation; otherwise it fits a model. `da` extracts results from fitted models; it does not refit. `p1` compares fixed gene neighborhoods against existing results. `pca` redraws saved coordinates only.

The reviewed cell-level ATAC partition is stored as `data/core2_reference/reviewed_cell_clusters.csv.gz`. On rebuilding, old descriptions are attached only if the same nuclei have exactly the reviewed cluster assignments. Otherwise descriptions remain unreviewed; cluster ID numbers alone never justify reusing annotations. RNA populations still require exact Core 1 cell matching. A changed partition needs renewed ATAC description review.

The user has now completed a **fresh real-data downstream run in RStudio**, starting from the four existing unfiltered QC objects. All ten comparisons matched the recorded run, including exact pseudobulk counts and reviewed cell-level ATAC assignments; differential-statistic tolerance was 1e-7. This does not include remapping FASTQ or recalculating QC metrics. Original published results were not overwritten. See results/core2/diagnostics/manual_rerun for evidence.

## RStudio entry for a separate validation run

Open the project in RStudio. Source the entry once, then invoke it explicitly:

```r
source("scripts/core2/validate_in_rstudio.R")
# Supply the existing local reconstruction workspace path as task_dir.
validation <- core2_validate_in_rstudio(
  project_dir = getwd(),
  task_dir = Sys.getenv("P3_CORE2_WORK"),
  validation_name = "core2_validation_rstudio_02",
  dry_run = TRUE
)
```

Set `P3_CORE2_WORK` to the actual workspace first, or supply its absolute path
directly as `task_dir`. `dry_run = TRUE` lists stages and destination only;
change it to `FALSE` to run. A new destination name avoids replacing the
completed validation run. Existing destinations are refused by default.
Sourcing alone starts no analysis. The repository entry only adapts path defaults
from the executed local entry, whose exact source is archived in
`results/core2/diagnostics/manual_rerun/executed_entry_point.R.txt`.

## Redraw only the reviewed PCA

In RStudio, configure `task_dir` as your existing local reconstruction workspace:

```r
source("scripts/core2/08_redraw_PCA.R")
core2_redraw_PCA(
  analysis_dir = file.path(task_dir, "outputs", "core2_analysis"),
  output_dir = file.path(getwd(), "results", "core2", "figures", "PCA_equal_units")
)
```

This uses saved coordinates and explained variances. It does not refit PCA, QC or DESeq2. This saved-coordinate redraw was performed during repository preparation, and the corrected images are included. The original stretched plots are not presented as the main final figures.

## QC and upstream provenance

The repository QC edition now also follows source-to-define behavior:

```r
source("scripts/core2/upstream/core2_batch_qc.R")
run_core2_batch_qc(task_dir = Sys.getenv("P3_CORE2_WORK"))
```

This explicit call scans fragments and calculates QC; it is not required to read completed results. Unlike the earlier teaching copy, simply sourcing the repository edition does not start a batch. NucleosomeSignal is retained to preserve the recorded metric implementation despite its deprecation warning.

`scripts/core2/upstream/` contains the recorded reconstruction/counter/audit code with path/tool configuration adapted for repository use. It is an advanced local recovery route, not a one-command clean installation. Restore SRA inputs, the four processed matrices, TAIR10 references, paired ARC barcode lists and their source records first. Public run metadata and whitelist provenance are included under data/core2_reference. The GTF is from Ensembl Plants release 62; record and verify its local checksum. Reference and artifact hashes are stored alongside the results.

Original tool versions were Chromap 0.3.2-r518, SRA Toolkit 3.4.1 and samtools/htslib 1.24 on Apple arm64. `P3_CORE2_TOOLS`, `P3_SRA_TOOLS`, `P3_RSCRIPT` and `P3_CXX` override tool paths where supported. Platform portability has not been established by a new upstream run. The remaining-sample pipeline reuses the Control 1 reference, so Control 1 comes first.

Explicit examples after restoring those inputs:

```sh
Rscript scripts/core2/upstream/core2_01_common_peaks.R
python3 scripts/core2/upstream/core2_control1_upstream.py all
python3 scripts/core2/upstream/core2_download_remaining.py
python3 scripts/core2/upstream/core2_run_remaining.py
python3 scripts/core2/upstream/core2_finalize_four_samples.py
```

Some helpers contain top-level execution and are command-line programs, not files to source for inspection. Download/export/mapping can take substantial disk space and time. The independent Control 1 audit refuses to overwrite its existing audit directory. Read the source before rerunning these recovery commands; no upstream rerun was performed during this preparation.

## Restore checkpoints

The local reviewed checkpoint is `workspace/outputs/core2_analysis/atac_annotated_checkpoint.rds`. It preserves sample, condition, barcode, RNA population/confidence, ATAC cluster and reviewed descriptions. Pseudobulk is `atac_pseudobulk_unfiltered.rds`; fitted models are in sample_diagnostics/<population>/deseq2_model.rds. RDS and fragments remain local, so the Git repository is a code/results record, not a complete sequencing-data backup.
