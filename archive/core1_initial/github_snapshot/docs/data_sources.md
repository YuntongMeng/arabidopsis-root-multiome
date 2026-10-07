# Core 1 inputs and artifact layout

RNA input accession: [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495). The sample/accession mapping is archived in `data/rna/sample_metadata.tsv`. This file list follows the four matrices used in the completed analysis; no ATAC inputs are required by this Core 1 workflow.

| Sample | GEO sample | Required file stem | Directory |
| --- | --- | --- | --- |
| control_rep1 | GSM7504014 | GSM7504014_Control1 | data/rna/control_rep1/ |
| control_rep2 | GSM7504015 | GSM7504015_Control2 | data/rna/control_rep2/ |
| stress_rep1 | GSM7504016 | GSM7504016_Case1 | data/rna/stress_rep1/ |
| stress_rep2 | GSM7504017 | GSM7504017_Case2 | data/rna/stress_rep2/ |

For each stem, place `*_matrix.mtx.gz`, `*_barcodes.tsv.gz` and `*_features.tsv.gz` in the listed directory. Feature column 1 and barcode column 1 are used. The original local download archive was named `GSE235495_RAW.tar`.

P1 candidate input: `data/project1/FULL_interaction_DEGs_padj0.05_log2FC1.tsv` (137 TAIR IDs).
P1 GO input: `data/project1/FULL_interaction_GO_enrichment_positive.tsv` (7 terms). These are existing Project 1 outputs; Core 1 does not rerun enrichment.

The main repository includes scripts, notes, fixed reference inputs and compact result tables. The release asset `core1-complete-artifacts.zip` contains the full curated publication snapshot, including all archived Core 1 results, successful-rerun tables/figures/diagnostics and the lightweight `core1_existing_results.rds`. Restore its data/ and results/ contents into the clone to supply the local artifacts needed by the cached-result paths and the historical-partition correspondence check. Keep the current repository scripts and documentation; the release archive preserves the original verified publication snapshot. The file checksums in docs/core1_file_checksums.sha256 describe that snapshot, not later documentation revisions. It does not contain original matrices or full Seurat objects.

The original 2026-09-30 full Seurat checkpoint is not distributed. It is unnecessary for the full `run_core1_full()` path, which produces stage 1 directly from the raw RNA matrices. The optional cached stage-2 entry point requires that original local checkpoint or an explicitly supplied new stage-1 checkpoint.
