# Project 3: scope, status and history

Decision date: 2026-10-15. I decided to conclude Project 3 at Core 2
and publish the reviewed work to GitHub. Core 3 is cancelled, not deferred.

## Timeline

- **2026-10-06, Core 1 rebuilt.** I found that the original RNA input included
  ATAC peaks and rebuilt Core 1 from gene-expression rows only. It retains
  21,550 predicted singlets; unintegrated clustering is primary and Harmony is a
  sensitivity comparison. I reran stages 01–07 from the original inputs, and all
  eight structural checks passed.
- **2026-10-13, Core 2 checked.** The interactive Core 2 downstream steps were
  assembled into scripts and rerun from the four saved QC objects; all ten
  recorded comparisons matched. A read-only scientific audit passed 111 checks.
- **2026-10-15, project concluded.** Project 3 was published with Core 1 and
  Core 2, and Core 3 was cancelled (see below).
- **2026-10-18, supplement published.** Model diagnostics, a counts>=5
  sensitivity run and RNA-response examples were added in a
  [supplementary review](core2/supplemental_review.md). The primary results did
  not change.

The current version is a documented working analysis, not a claim of a single
correct solution. Earlier assumptions, errors, alternatives and remaining
uncertainty are recorded separately from the current workflow.

## Completed scope

- [Core 1](../cores/core1/README.md): corrected gene-only RNA QC, predicted-doublet handling, clustering, provisional annotation, population-level inference and Project 1 context.
- [Core 2](../cores/core2/README.md): common-region ATAC recounting, QC, exact-nucleus RNA pairing, ATAC clustering and descriptive annotation, replicate-aware pseudobulk differential accessibility and Project 1 neighborhood context.
- The assembled Core 2 downstream workflow was rerun manually from existing QC objects; all ten recorded comparisons matched. Full FASTQ-to-final reconstruction is not claimed.

## Why the project stops here

Same-nucleus RNA–ATAC identity matching, composition comparisons and marker
region checks already belong to Core 2. Relabeling those steps as Core 3 would
duplicate work. The limited differential evidence does not justify adding a
mechanistic story solely to extend the project. It also does not establish
that further questions are impossible or that stress effects are absent.

The final record keeps the one within-population FDR-significant ATAC candidate,
the absence of significant matched nearby ATAC support for the Project 1 genes under the fixed
rule, uncertain mixed-population labels, replicate limitations and the unresolved
experimental cause of Control 1's low-signal subgroup.

No new Harmony/LSI integration, equal-depth resampling, contamination correction,
GO analysis, motif analysis, peak-to-gene linkage, co-accessibility, trajectory
or regulatory network was added during closure. Methods I considered but did not
run are not presented as completed analyses.

## Record preservation

The old [Core 3 placeholder](../cores/core3/README.md) is retained only to explain
the cancelled plan and preserve older links. It is not an active project core.
Historical plans remain historical; no Core 3 workflow or result is implied.
Raw sequencing inputs, fragments, checkpoints and runtimes stay local; Core 1
analysis scripts, result tables and the frozen archive are preserved. This
closure updates scope and documentation, not biological results or thresholds.
