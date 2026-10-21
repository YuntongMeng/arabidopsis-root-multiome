# Core 1 corrected version: rerun and comparison with the initial version

2026-10-06 · Project 3 · Run completed. The summaries and eight structural
checks from my manual 01–07 rerun were checked. The final publication changes
are in `decision_log.md` and `manual_rerun_validation.md`.

## Summary

**These results are not a reproduction of the initial version, and the initial
biological conclusions have to be revised to match them.** This round started
again from the raw multiome matrices: Gene Expression filtering, QC,
scDblFinder, PCA/Harmony, markers/annotation, DESeq2, Project 1 Fisher tests
and GO. Frozen results were not reused as new results.

- **Still supported:** the population × biological sample pseudobulk design;
  the Project 3 osmotic-stress context for the Project 1 candidates;
  significant enrichment of the candidates in Pericycle relative to the
  testable genes of the same population; DE direction and a sizeable share of
  genes in several populations agree with the initial version.
- **Withdrawn or rewritten:** old cluster numbers and labels; "Pericycle has
  the highest overlap proportion"; "Pericycle is involved in all 7 GO terms";
  "GO-linked genes are all downregulated".
- **New interpretation limit:** the Project 1 interaction candidates are
  significantly enriched relative to the ordinary gene background, but neither
  branch shows extra enrichment relative to **Project 1 non-interaction drought
  DEGs**. This does not let me claim that Project 3 validates BRL3
  interaction specificity.
- **A finding that corrects an earlier review:** none of the 2,000 HVGs in the
  old checkpoint were coordinate peaks (**0**). The ATAC rows did contaminate
  the input, QC, normalization denominator and marker tests, but this
  checkpoint does not support "peaks entered the HVGs in large numbers and
  drove the PCA".

This round changed input handling, QC, doublets and annotation at the same
time, and added a Harmony branch. Differences cannot all be attributed to the
ATAC filtering alone. No one-factor-at-a-time ablation was done.

## 1. Input correction and QC

Each original sample contained 19,858 / 24,317 / 24,490 / 25,001 Peaks rows.
Rows were kept strictly when **column 3 of the features file equals Gene
Expression**, and matrix row order, IDs, dimensions, non-negative integer
counts and unique barcodes were checked one by one. Every sample has 32,833
genes, and the shared gene universe is also 32,833; after filtering there are
0 coordinate-style peak IDs. The rows were not accepted as genes just because
the count is close to the total gene number.

I first looked at gene-only quantiles, densities and scatter plots, then set
within-sample log1p-MAD thresholds: 3 MAD for low features/counts with a floor
of 500; 4 MAD for extremely high counts; 3 MAD for high mt. The resulting
minimum features are 500/500/583/602, and the maximum mt is
10.41/4.76/4.53/4.57%. cp is diagnostic only. The old mixed-feature limits of
800–15,000 features / 50,000 counts were not reused.

| Sample | Input nuclei | Removed by QC | Predicted doublets removed | Singlets retained |
| --- | --- | --- | --- | --- |
| control_rep1 | 3431 | 85 | 236 | 3110 |
| control_rep2 | 6631 | 180 | 666 | 5785 |
| stress_rep1 | 7827 | 94 | 917 | 6816 |
| stress_rep2 | 6585 | 71 | 675 | 5839 |

Of 24,474 input nuclei, QC removed 430 and scDblFinder removed another 2,494,
leaving **21,550**. Each capture was processed separately; the rate prior was
0.008 per 1,000 captured cells, and versions and random seeds are saved. These
are predicted doublets, not confirmed ones; the actual loading rate is
unknown, and homotypic doublets may remain.

The initial version retained 24,277 nuclei, of which 21,488 are shared with
the new version; 2,789 are unique to the old set and 62 to the new set.
Rebuilding the old barcode set with the original mixed-count rules matches the
frozen metadata exactly. The 800-gene floor and the mt 4-MAD alternative were
only checked for how many nuclei they retain; a full QC-parameter sensitivity
rerun is not claimed.

**Ambient RNA:** the local GEO `RAW.tar` contains only the same 12 processed
matrix/features/barcodes files, with no verified raw or empty droplets. So I
did not force SoupX, and cross-lineage expression was not treated as a
contamination fraction. The canonical DotPlot can be inspected, but ambient
RNA cannot be ruled out.
[GEO data source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495)

## 2. PC stability and the Harmony decision

Both branches use the gene-only singlets and the RNA assay, LogNormalize with
10,000, and an initial 2,000 HVGs; after removing 1 organelle HVG, the
embedding uses 1,999. 1,524 of these genes are shared with the old HVGs. Raw
counts keep all gene features, and per-sample count totals were confirmed to
be identical before and after building the final object.

Both branches were computed at 10/15/20/25/30/40/50 dimensions with
resolution = 0.5. I chose 25 PCs for local stability in the 20–30 range and
for parsimony: unintegrated ARI(20,25) = 0.809 and ARI(25,30) = 0.880; Harmony
0.886 and 0.903. 40–50 PCs form another stable region whose partition still
differs clearly from 25 PCs, so the whole PC range cannot be called stable.

At 25 PCs the unintegrated branch happens to give **18 clusters** again, and
Harmony gives **17**; this does not mean the old partition was reproduced. The
old/new unintegrated ARI on shared nuclei is **0.639**; new unintegrated vs.
Harmony is 0.622.

| Diagnostic | Unintegrated | Harmony |
| --- | --- | --- |
| Replicate mixing ratio within the same condition and the same higher-confidence atlas identity, group mean | 0.876 | 0.984 |
| Neighbourhood agreement of broad atlas identity | 0.485 | 0.485 |
| Fraction of same-condition neighbours | 0.859 | 0.600 |
| Spearman between raw and neighbour-smoothed expression of sparse stress markers | 0.229 | 0.184 |

I chose **unintegrated as the primary analysis and Harmony as a sensitivity
analysis**. Harmony improves mixing of the control replicates, but it does not
improve this identity-agreement metric to match, and it weakens
condition-associated structure. Here sample is nested in condition, and no
independent ground truth can show whether the removed variation is all batch
or all biology. The four stress diagnostic genes are sparse, so that
correlation is only supporting evidence. Low-dimensional separation by
condition does not by itself show that biological signal was preserved
correctly. Raw-count DE was run in both branches, and both are reported below.

## 3. Annotation: sourced, with uncertainty kept

The main panel comes from
[Shahan et al., 2022, Developmental Cell](https://doi.org/10.1016/j.devcel.2022.01.008):
the original **Data S2K** 14 classes with the top 50 signature genes each;
**S2F/G** canonical markers with their original sources and spatial
expression ranges; and **S2J** cell-cycle markers. The original xlsx, the
extracted panels, the citations and the reasoning for each new cluster are in
`data/core1_corrected_reference/`. This is not the old placeholder list under
a new name.

AddModuleScore is supporting evidence only and is judged together with strict
markers from this dataset and the canonical DotPlot. Strict markers use one
standard: Wilcoxon, min.pct = .25, logFC = .25, returning p < .01; the marker
tables also include adjusted p values. Cell-level marker p values do not
replace replicate-based DE. The canonical and atlas panels come from related
literature, so they are not fully independent validation.

- Broad lineages merge subgroups that share identity support but differ in
  developmental or condition state. For example, the control-enriched and
  stress-enriched unintegrated Pericycle subgroups are merged before DE.
- The Atrichoblast/Cortex/Endodermis/Trichoblast `_lineage` groups contain
  stronger and weaker subgroups; weaker subgroups are not upgraded to a
  definite mature identity.
- **Root_cap** has evidence such as SMB/BRN1/BRN2, but shared markers cannot
  reliably separate LRC from columella.
- **Meristematic_unresolved** and **Epidermal_initial_like** stay low
  confidence; **Dividing_meristematic** is a state with cell-cycle evidence,
  not a single lineage.
- Phloem, small unresolved groups and sample-specific states keep their
  descriptions; groups that do not meet the per-replicate threshold do not
  enter formal DE.

Cells from the old LRC_like group are largely spread across the new
Atrichoblast_lineage and Cortex_lineage. The old LRC-specific interpretation
therefore cannot be carried over. 88.6% of cells in the new Root_cap come from
the old Columella, which supports continuity at the root-cap level, but the
old fine labels are not forced back.

## 4. Pseudobulk and DE

Each population × biological sample sums **raw integer gene counts**. Within
the shared gene universe, genes are pre-filtered at row total ≥ 10 across the
four samples. DESeq2 uses `~condition`, stress vs. control, with padj < .05
and |log2FC| ≥ 1. rep1/rep2 are not treated as a paired design, and cells are
not treated as replicates.

The primary analysis has **9 eligible populations** and Harmony has **10**.
The smallest per-sample count in every selected population is ≥ 100, so the
≥ 25 / ≥ 50 / ≥ 100 eligibility rules give the same set in this run. With n = 2
per condition, conclusions remain exploratory. DE in low-confidence
populations can only be read under those conservative labels.

| Population | Primary DEGs | Harmony DEGs |
| --- | --- | --- |
| Atrichoblast_lineage | 1.56e+03 | 1.55e+03 |
| Cortex_lineage | 1.46e+03 | 846 |
| Dividing_meristematic | 395 | 388 |
| Endodermis_lineage | 524 | 439 |
| Epidermal_initial_like | 403 | 440 |
| Meristematic_unresolved | 504 | 316 |
| Pericycle | 637 | 630 |
| Procambium_like | NA | 239 |
| Root_cap | 148 | 155 |
| Trichoblast_lineage | 1.75e+03 | 1.53e+03 |

NA for Procambium_like in the primary analysis means the label exists only in
the Harmony branch, not that it has zero DEGs.

These populations are redefined, so a row cannot be treated as a strict
one-to-one match with an old group of the same or a similar name. All
fitted / valid-pvalue / finite-padj denominators, raw pseudobulk, sample
correlations, cell counts, size factors and full DE results are saved.

**Stability between the two new branches:**

| Population | log2FC Spearman | DEG Jaccard | Cell membership Jaccard |
| --- | --- | --- | --- |
| Atrichoblast_lineage | 0.946 | 0.717 | 0.827 |
| Cortex_lineage | 0.851 | 0.339 | 0.713 |
| Dividing_meristematic | 0.99 | 0.882 | 0.949 |
| Endodermis_lineage | 0.947 | 0.787 | 0.908 |
| Epidermal_initial_like | 0.968 | 0.809 | 0.819 |
| Meristematic_unresolved | 0.896 | 0.48 | 0.528 |
| Pericycle | 0.901 | 0.705 | 0.725 |
| Root_cap | 0.963 | 0.836 | 0.904 |
| Trichoblast_lineage | 0.988 | 0.825 | 0.858 |

Every gene significant in both branches has the same sign in both. But the DEG
Jaccard is only 0.339 for Cortex and 0.480 for Meristematic_unresolved, so
membership and the number of detected genes clearly depend on the branch.
Pericycle has 524 shared DEGs between branches, with log2FC Spearman = 0.901;
the dividing state and Root_cap are more stable. "Same direction" does not
mean every population is fully stable.

Among populations whose membership is close to the old version,
Dividing_meristematic has log2FC Spearman = 0.945 and DEG Jaccard = 0.709;
Pericycle 0.873 / 0.631; and the new Root_cap vs. the old Columella
0.900 / 0.738. The detailed comparison includes membership intersections, so
it does not compare names alone.

## 5. The 137 Project 1 candidates and Fisher enrichment

The main Fisher background is, for each population, the genes with a
**non-missing padj that actually enter significance testing**. The table is
Project 1 / non-Project 1 × stress DEG / non-DEG, with a one-sided enrichment
test and BH correction across populations within each branch. A sensitivity
version using all fitted genes, including independently filtered ones, is
saved separately; the denominators are not mixed.

| Population | Project 1 stress DEGs | Testable Project 1 genes | OR | Fisher BH |
| --- | --- | --- | --- | --- |
| Atrichoblast_lineage | 26 | 108 | 3.28 | 3.83e-06 |
| Cortex_lineage | 19 | 87 | 2.38 | 0.00177 |
| Dividing_meristematic | 14 | 73 | 9.13 | 2.25e-08 |
| Endodermis_lineage | 11 | 69 | 3.7 | 0.000667 |
| Epidermal_initial_like | 11 | 70 | 6.59 | 5.58e-06 |
| Meristematic_unresolved | 17 | 87 | 7.74 | 1.2e-08 |
| Pericycle | 19 | 97 | 5.51 | 8.48e-08 |
| Root_cap | 5 | 56 | 6.34 | 0.00177 |
| Trichoblast_lineage | 30 | 102 | 3.47 | 4.8e-07 |

All 9 primary populations, and all 10 Harmony populations, are enriched
relative to the ordinary background. This does not show precise cell-type
specificity, especially for low-confidence labels. Each Fisher family is
corrected separately; the alternative tests are not pooled into one global BH.

**Pericycle:** the primary analysis has 19 Project 1 stress DEGs (6 up, 13
down), 19/120 fitted = 15.83%; in the actual Fisher background it is 19/97,
OR = 5.51, BH = 8.48×10⁻⁸. Harmony also has 19 (6 up, 13 down), 19/117 fitted,
19/78 finite-padj, OR = 5.93, BH = 1.45×10⁻⁷. The old figure of 23/118 should
no longer be cited. The highest fitted overlap proportion in the primary
analysis is now Trichoblast_lineage, 30/115 = 26.09%, not Pericycle; ranking
by proportion is still not a test of enrichment specificity.

**A more suitable control:** Project 1 has 1,480 WT drought DEGs. Removing
those that overlap the 137 interaction candidates leaves 1,442 control genes;
also requiring a non-missing interaction padj ≥ .05 leaves 1,389 strict
control genes. Against both controls, **no population in either branch reaches
BH < .05 for extra enrichment of the interaction candidates**. For Pericycle,
the primary analysis against the first control gives OR = 1.39, BH = 0.437;
Harmony gives OR = 1.36, BH = 0.688.

So I keep the description "the candidates converge with an independent
osmotic-stress response" and drop "interaction-specific validation" and
"Pericycle is unique or strongest". A non-significant result does not show
that the two sets are equivalent; no expression-matched control was used, so
confounding by expression level and detection power remains.

## 6. The seven GO annotations

After remapping there are still 7 terms, 36 term–gene memberships and 13
distinct genes. These are the **candidate members listed in the original
Project 1 enrichment table**, not the full set of genes annotated to each GO
term. No new Project 3 GO enrichment was run.

In both branches only **3 of 7 terms** have at least one overlap meeting the
current Project 1 × Project 3 DEG criteria: JA metabolism, long-chain fatty
acid metabolism and response to cold. Cold acclimation, response to water,
water deprivation and acid response are all 0. This means no candidate overlap
at these thresholds, not that the program is absent.

In Pericycle, cold response drops from 3 genes in the old version to **1** in
the primary analysis and **2** in Harmony. The JA and long-chain fatty acid
terms each have 1 downregulated gene, but they share members and are not two
independent findings. In Atrichoblast_lineage the overlap for these two terms
is now upregulated, so the old statement "GO-linked genes are all down" no
longer holds.

The full 7 × population matrix includes zero rows and columns, and pairwise
gene-set overlap and Jaccard values are exported. The seven terms overlap
strongly to begin with, so they cannot be described as seven independent
programs that were each confirmed.

## 7. Which old conclusions still hold

| Old conclusion or method | Status now |
| --- | --- |
| Framework: Project 1 → population → stress pseudobulk → contextualization | Kept |
| Biological sample as the statistical replicate; gene-only shared universe | Kept and checked directly |
| PC25, 18 clusters | PC25 has new stability support; 18 is a coincidence in the unintegrated branch, and the partition differs |
| Old 13 populations and labels | Redone; formal DE uses 9/10 populations, low-confidence labels read conservatively |
| Pericycle has strong stress/Project 1 context | Still supported in both branches, but not unique, highest or interaction-specific |
| Pericycle 23/118 | Replaced by the new branch results, with clearly separated denominators |
| Pericycle in all 7 GO terms; all GO overlaps downregulated | Withdrawn; replaced by the 3-term result and directions above |
| Harmony does not need to be considered at all | Replaced by a measured sensitivity comparison; unintegrated is currently primary |
| Peaks entered the HVGs in large numbers and drove the PCA | Not supported by the old checkpoint; should not be stated as fact |

## 8. Files, reproduction and review

The full local checkpoints in `results/core1/checkpoints/` (gene-only counts,
singlets, both branch objects and pseudobulk; about 5.6 GB) are not in Git.
The repository has the scripts, all CSV tables, figures, panel sources, method
records, run logs and sessionInfo. The original matrix files stay where they
were.

**Entry point:** from the project root, `Rscript scripts/run_core1_corrected.R help`.
Stages can be run one at a time (audit / qc / embeddings / evidence / de /
context / compare / diagnostics) or explicitly with `all`. Ambiguous R
functions are always written as `package::function()`. Standard R libraries
are used by default; an untracked local `P3_R_LIB` setting can override this.
The common helper no longer contains personal absolute paths.

The manual annotation is the product of this evidence review and is tied to
the partition MD5. A future change in the partition requires a new review;
labels should not be applied without that check. The old mixed-input stage-1
entry point now warns against accidental use, and the old code and the frozen
archive are kept.

Checked: syntax of the 10 new R files; 8 structural invariants; 152 Fisher p
values recomputed by independent hypergeometric sums, and BH recomputed with an
independent sorting implementation; 27 frozen files match their existing
SHA-256 records. Figures were checked for readability.

Main methods and data sources:
[Shahan atlas and supplementary tables](https://pmc.ncbi.nlm.nih.gov/articles/PMC9014886/),
[scDblFinder documentation](https://plger.github.io/scDblFinder/articles/scDblFinder.html),
[Project 3 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495),
[root osmotic-stress marker study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10573044/).
