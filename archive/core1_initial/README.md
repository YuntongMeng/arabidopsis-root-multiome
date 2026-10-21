# Core 1 initial analysis — historical archive

These files preserve the earlier reasoning and implementation. They are not the
current workflow. Mixed RNA/ATAC input invalidates the original grouping and downstream
interpretation as formal results; see ../../docs/core1/analysis_story.md.

- The complete first published version, including its README, runner and interpretation, is the [core1-v1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1). It is historical only.
- `scripts/`: original 00–04 scripts; the old preprocessing entry is disabled.
- `docs/`: original file checksums. The first version's summary, integrity notes and review checklist were working notes and are kept locally.
- Earlier personal working notes remain local.
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
