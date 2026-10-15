# Core 2 reference records

- `common_peaks.bed`: 27,677 nuclear intervals, zero-based half-open. No peak calling is performed by reading this file.
- `reviewed_cell_clusters.csv.gz`: capture-prefixed cell IDs and reviewed ATAC cluster assignments. Binds descriptive labels to the recorded partition; changing cluster numbers/assignments prevents automatic label reuse.
- `atac_cluster_annotation_review.csv`: RNA compositions, partial ATAC marker support and limits for five descriptive clusters.
- `remaining_run_metadata.tsv`: public run metadata used to verify ATAC accession/capture identity. Copy to workspace/work/remaining_run_metadata.tsv when restoring the original upstream layout.
- `barcode_sources.json`: URLs, row counts and SHA-256 for the paired ARC whitelist mirror. Download both lists to workspace/work/barcode_reference with their recorded basenames and preserve row order. Save a copy of this JSON there as sources.json. The lists are a Teichmann-lab mirror, not a claimed direct 10x download.

The gene marker references remain in `data/core1_corrected_reference`; no new marker truth set was invented for ATAC. TAIR10 FASTA/GTF, SRA, whitelist payloads and runtime binaries remain external. [Restore/run instructions](../../docs/core2/running.md).
