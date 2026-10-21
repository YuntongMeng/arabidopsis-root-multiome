# Core 2 validation scope

Repository preparation date: 2026-10-13. Saved biological results were not refitted or reselected.

## Recorded original checks

- Four-stream FASTQ synchronization, lengths, barcode extraction/translation, alignment fragment ordering, coordinates, deduplication and tabix access were checked during upstream reconstruction.
- All four common-peak matrices have identical row coordinates/order and original per-capture barcode columns. Counter fixtures and sampled real-region recount checks were recorded.
- Control 1 received a separate full-fragment audit and full-matrix comparison. Scope and exclusions are in [Control 1 review](control1_review.md).
- The original QC, pairing, embeddings, annotation, pseudobulk, DESeq2 and P1 outputs were saved in local checkpoints and tables.

## Checks performed during repository preparation

- R scripts parse, Python scripts compile, and runner help executes without fitting or downloading anything.
- Synthetic cases check ATAC threshold boundaries, exact RNA call matching, RNA doublet/not-tested exclusion, metadata reordering with name-cache attributes, missing-cell rejection, pseudobulk count conservation and omission of unobserved combinations.
- Existing pseudobulk sample/population totals and fitted-model inputs are checked against saved counts; exported summaries are checked against full result tables.
- The published 14,146-nucleus metadata, common BED, retention totals, P1 mapping summaries and result flags are validated by `scripts/verify_core2.py`.
- Corrected PCA figures were rendered from **existing coordinates**, with equal x/y unit scale and independent-panel interpretation. No PCA/model recomputation was involved.
- Documentation links, publication file sizes, ignored raw data and Git whitespace/status are checked before publication.

Detailed checks are recorded in results/core2/diagnostics/repository_validation.json and results/core2/diagnostics/saved_model_input_validation.csv. Publication file hashes are in results/core2/diagnostics/publication_manifest.csv.

## Manual downstream rerun

After repository preparation, I ran the supplied RStudio entry against the real data. It rebuilt selection/merge, embeddings, RNA annotation + GeneActivity, pseudobulk, DESeq2 diagnostics/models, DA export, P1 contextualization and PCA. Outputs were written to a separate validation directory.

All ten comparisons matched: common peak rows, observed sample/population groups, every raw pseudobulk entry, reviewed cell-level ATAC assignments, differential statistics for five populations (numeric tolerance 1e-7), and P1 summary. The comparison does not assert bitwise equality of every embedding, GeneActivity entry or figure pixel. The optional marker/coverage plot stage was not invoked. Evidence: results/core2/diagnostics/manual_rerun.

## Independent scientific audit

The [read-only audit](../../results/core2/diagnostics/scientific_audit/README.md)
checked source feature types, all final raw count submatrices, exact paired
selection and labels, all 46 pseudobulk sums, five saved model inputs and
exported statistics, P1 interval arithmetic, fragment-to-nucleus mapping and
shared TSS coordinates. The count/model/mapping checklist has 111 passing
items. It did not refit models, remap FASTQ or independently recount TSS scores.
No principal scientific error was found within this scope; labels, causal
targets, experimental causes and biological truth remain limited as documented.

## What has not been established

- The assembled downstream workflow was rerun from existing QC objects, and all ten comparisons matched. This is not a fresh FASTQ-to-final-results run, and optional marker/coverage plots were not rerun by that entry point.
- Repository path adaptations have not been validated by fresh whole-upstream reconstruction on another platform.
- No new TSS audit or FASTQ remapping was performed during the independent Control 1 audit or this preparation.
- No ground-truth cell labels, experimental cause of Control 1's subgroup, causal peak target, independent experiment-wide discovery validation, new GO result or regulatory network is claimed.

Original vs repository source fingerprints document code provenance. Recorded outputs remain the authoritative evidence for the reported numbers; fresh results, changed annotations or different environments need renewed review.
