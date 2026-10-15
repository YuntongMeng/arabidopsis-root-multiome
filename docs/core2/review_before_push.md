# Core 2: submission record and inspection

The reviewed Core 2 analysis was prepared locally and manually validated in
RStudio. On 2026-10-15, the user explicitly authorized publication to GitHub and
cancelled Core 3. Project 3 concludes at Core 2; see the
[final scope](../project_scope.md). The local task workspace records the actual
push outcome and remote commit after upload. Git history records the commits.

## Scientific review remains visible

- [Conclusions](analysis_story.md): one candidate peak, limited power, no added significant nearby ATAC support for P1.
- [Control 1 review](control1_review.md): a pre-existing low-signal subgroup; specific experimental cause unresolved.
- [Cluster composition](../../results/core2/figures/02_RNA_composition_of_ATAC_clusters.png), [corrected PCA](../../results/core2/figures/PCA_equal_units/all_populations_PCA_equal_units.png), and [candidate counts](../../results/core2/figures/08_candidate_peak_sample_counts.png).
- Mixed cluster descriptions remain provisional, Root_cap exploratory, and nearby genes unproven targets. Core 3 is cancelled; the project concludes at Core 2.

## Reproducibility and file boundaries

Read [running.md](running.md), [validation.md](validation.md) and the
[directory guide](directory_map.md). All ten user-run downstream comparisons
matched. Upstream and QC reconstruction were not repeated. Raw sequencing
data, fragments, RDS checkpoints and runtimes remain local and excluded from
the commit. Core 1's scripts, results and frozen archive remain unchanged.

These commands inspect the local submission without biological analysis:

```sh
git status --short
git log -1 --oneline
git diff --check
Rscript scripts/run_core2.R help
python3 scripts/verify_core2.py
```

## Authorized publication procedure

Recheck the branch, working tree and remote state. If `main` still contains the
prepared commit and the remote has not diverged, upload the existing commit:

```sh
git push origin main
```

Do not force-push, silently discard later edits, or rerun the analysis as part
of uploading. A future remote or local change may require handling that new
state before push. Publication was explicitly authorized by the user on 2026-10-15; actual success must be checked against the remote commit.
