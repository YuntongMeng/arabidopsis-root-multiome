# Core 2: supplementary diagnostics and exploratory checks

Supplementary revision completed on 2026-10-15. The user reviewed the separate private PDF and authorized GitHub publication on 2026-10-18.
The original counts>=10 analysis remains primary. This supplement does not change
QC, cell identities, sample membership, FDR, or the two biological replicates per condition.
Project 3 still concludes at Core 2; further Control 1 attribution is paused.

## 1. Ten standalone diagnostic figures

The [model-diagnostic folder](../../results/core2/supplement/model_diagnostics/) has
one p-value histogram and one dispersion figure for each of five RNA populations.
Each is a full-size figure, available as PNG and PDF. They are supplementary
diagnostics, not ten additional main-result panels. Each dispersion figure includes
a separate proportion bar, calculated from the model rather than inferred from
overlapping points. Axis limits are shared across the five dispersion figures.

| RNA population | Raw p<.01 | Initial dispersion at boundary | Fit |
| --- | --- | --- | --- |
| Atrichoblast_lineage | 203 | 81.3% | local |
| Cortex_lineage | 10 | 80.7% | parametric |
| Endodermis_lineage | 8 | 81.4% | local |
| Pericycle | 5 | 81.9% | local |
| Root_cap | 1 | 81.3% | local |

The numerical boundary is 1e-8, the lower bound of the initial peak-wise estimate.
It is not a zero-count category, a final-estimate category or proof of model failure.
Boundary proportions are similar across these populations, even though their
replicate correlations and p-value distributions differ. Thus this proportion
alone does not explain their different results or identify Control 1 as the cause.

The p histograms show fewer small p values, especially outside Atrichoblast.
This is compatible with conservative behavior or weak information, but does not
on its own prove inflated dispersion, invalid calibration, or a particular sample's
causal role. Discrete counts, limited replication and differing signal strengths
also matter. Local dispersion fitting is a supported fallback, not an error.
[Count-stratified numeric diagnostics](../../results/core2/supplement/model_diagnostics/count_stratified_diagnostics.csv)
are retained without adding another crowded multi-panel plot.

## 2. A declared broader peak-count filter

The primary prefilter is counts>=10 in at least two of the four sample pseudobulks.
The additional analysis uses **counts>=5 in at least two samples**, chosen as a
simple broader setting before inspecting its results. It refits DESeq2 with the
same design ~condition, poscounts size factors, stress/control contrast, alpha=.05,
population-wise BH and independent filtering. No requirement was added that a
peak have counts in all four samples; condition-specific regions remain eligible.

| RNA population | Primary >=10: peaks | Sensitivity >=5: peaks | Primary FDR<.05 | Sensitivity FDR<.05 |
| --- | --- | --- | --- | --- |
| Atrichoblast_lineage | 27634 | 27674 | 1 | 1 |
| Cortex_lineage | 27384 | 27631 | 0 | 0 |
| Endodermis_lineage | 25508 | 27396 | 0 | 0 |
| Pericycle | 23757 | 27297 | 0 | 1 |
| Root_cap | 2487 | 14664 | 0 | 0 |

The original Atrichoblast region 4-5457000-5457860 remains significant
(primary log2FC 1.86037, padj .004526; sensitivity log2FC 1.86044, padj .004335).
Pericycle gains one FDR-significant region in the sensitivity run:
2-3313938-3315145, log2FC .74298 (about 1.67-fold), padj .005350.
It does **not** meet the unchanged estimated twofold priority criterion.
It was already present in the primary model (log2FC .73160, padj .999926), so this
is a filter-sensitive inference after refitting, not simply a newly eligible peak.
Its standard error changes from .22399 to .14280; the change is not explained by
BH adjustment alone. The four raw counts stay 71 / 197 / 147 / 146.
Both stress normalized counts exceed both control values in both fits, a descriptive
pattern that does not remove the sensitivity limitation.

All retained peaks' four raw counts are identical between the two settings.
FDR is controlled separately within each population/run, not jointly across both
settings. The supplement does not replace the primary result or create a new
twofold candidate. Root_cap remains exploratory despite its larger eligible peak set.
[Threshold table](../../results/core2/supplement/sensitivity_count5/threshold_comparison.csv)
and [shared-peak changes](../../results/core2/supplement/sensitivity_count5/shared_peak_comparison.csv)
preserve both outcomes.

## 3. Existing Core 1 RNA responses versus nearby ATAC results

This adds response context, not another round of identity annotation.
For each population, select the first three primary/unintegrated (pca branch)
RNA genes with padj<.05 and |estimated log2FC|>=1, ordered by RNA padj,
descending |log2FC|, then gene ID. No ATAC statistic enters selection.
This gives 15 population-gene examples (14 distinct gene IDs), including both
up- and down-regulated RNA responses. This small set is illustrative, not a
genome-wide estimate of RNA/ATAC concordance.

Mapping reuses the declared positional rule: common peak overlapping gene body
OR strand-aware TSS window, 2000 bp upstream / 500 bp downstream. All matched
peaks, including those excluded from the ATAC model, are retained in the long table.
For illustration only, select one modelled peak by promoter overlap first, then
TSS distance, then coordinate ID; no ATAC p value or effect size enters the choice.
Effect-size panels use +/-1.96 model SE as descriptive intervals, not simultaneous
intervals or a gene-level ATAC test.

All 15 examples have nearby common peaks; 13 have at least one peak in the
primary ATAC model. None has a nearby primary FDR-significant peak. Two Root_cap
examples have only peaks excluded by the primary count filter; this is not a zero
ATAC response. RNA and ATAC estimated directions can agree or differ, and a
non-significant ATAC result is not evidence of no accessibility change.
RNA response is not a guaranteed positive control for ATAC change. Proximity
does not establish a regulatory target link.

Existing RNA statistics use Core 1's RNA-QC singlets; ATAC uses the smaller paired
ATAC/RNA-QC subset. The captures can be matched at nucleus level, but the two
existing differential analyses do not use exactly the same retained nucleus sets.
No new RNA fitting, correlation test or causal regulatory analysis was performed.
[Gene selection](../../results/core2/supplement/RNA_ATAC_response/RNA_selected_genes.csv),
[all nearby peak results](../../results/core2/supplement/RNA_ATAC_response/RNA_ATAC_peak_context.csv),
and [retained-nucleus scopes](../../results/core2/supplement/RNA_ATAC_response/retained_nucleus_scope.csv)
make these boundaries explicit.

## 4. What remains unchanged and what is paused

- TSS enrichment >1, nuclear fragments >=1000 and nucleosome signal <2 remain unchanged.
- No nucleus is treated as an independent biological replicate; no sample is removed.
- The original one-candidate result and P1 comparison are preserved.
- No additional Control 1 equal-nucleus sampling, RNA quality rerun, FASTQ remapping,
  GO analysis or Core 3 work was started.
- High prefilter organelle content is compatible with organelle carryover or
  deficient nuclear signal, but does not establish isolated organelle-only droplets
  or a specific laboratory mistake. Final C1 organelle median is 4.67%, not 73.27%.

## Reproducibility and local status

The new script is [09_supplemental_review.R](../../scripts/core2/09_supplemental_review.R).
Sourcing defines a function only; explicit invocation runs the supplement into a
new output directory and refuses overwrite. Five original models' raw inputs and
saved statistics were checked against the pseudobulk and tables (tolerance 1e-7).
The 18 consumed input files, including original models, primary results, Core 1
tables and references, have identical before/after hashes.
Large supplemental model RDS files remain local. A separate private Chinese PDF
review is held only in temporary storage; it is not part of the repository.

Run `python3 scripts/verify_core2.py` for the original results and
`python3 scripts/verify_core2_supplement.py` for the supplement. Neither refits models.
The 2026-10-15 analysis execution did not commit or push. On 2026-10-18 the user authorized publication of the reviewed supplement; Git history records the publication commit.
