# Core 1 initial analysis — historical archive

These files preserve the earlier reasoning and implementation. They are not the
current workflow. Mixed RNA/ATAC input invalidates the original grouping and downstream
interpretation as formal results; see ../../docs/core1/analysis_story.md.

- `github_snapshot/`: the complete lightweight previously published layout, including its old runner and interpretation; historical only.
- `scripts/`: original 00–04 scripts; the old preprocessing entry is disabled.
- `docs/`: original summary, checksums, integrity notes and manual-review checklist.
- Earlier personal working notes remain local. The previously public README and four-script workflow are retained under `github_snapshot/`.
- `results/reproducibility/`: original full Seurat checkpoint and PC assignments (large RDS stays local and is not distributed).
- `results/tables/`, `results/figures/`: legacy output locations.
- `results/core1_rerun_20261005_*`: earlier rerun/intermediate output folders; not promoted as current results.
- `../../results/core1_frozen/`: untouched frozen original results, kept in their original location.

Historical documents may use the original project-root paths. Those strings describe
the earlier layout; they are not current runnable instructions. Scripts and documents
were moved without rewriting their analytical history. Original inputs were not moved.
The current comparison scripts explicitly use this archive for the old full checkpoint.
The actual current entry is ../../scripts/run_core1_corrected.R.

Public historical notes omit machine-specific private paths and local history metadata. Local originals remain unchanged; analytical parameters and results are preserved.
