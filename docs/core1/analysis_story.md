# Core 1: from the first reading of the data to the current analysis

This is the main thread for following the analysis. The current workflow has
been run in full; reorganizing the directories did not recompute any biological
result. "Current" means the choice the evidence supports now, not that every
question has a single answer.

## 1. How I first read the data

I first used the multiome feature matrix as RNA input, reading it by feature ID
without separating Gene Expression rows from Peaks rows in the third column.
QC, clustering, annotation and pseudobulk were built on that input. I also
looked at coordinate-style markers but did not trace where they came from soon
enough. That code and reasoning stay in `archive/core1_initial`, and the
original results stay in `results/core1_frozen`. Old conclusions should not be
cited as current results.

## 2. What changed my mind

Coordinate-style markers and unusual feature counts led me to re-audit the raw
features files. Each sample has 32,833 gene features plus 19,858–25,001 ATAC
peaks. The new reader first aligns matrix rows with feature rows, then keeps
only rows whose third column is exactly `Gene Expression`. No coordinate peak
IDs remain in the RNA data.

This also corrected my view of what the error did: none of the old 2,000 HVGs
were coordinate peaks. The mixed input certainly affected QC, the normalization
denominator and marker tests, but I cannot claim that peaks directly drove the
old PCA. The new workflow changed several steps at once, without a
one-factor-at-a-time ablation, so not every difference can be attributed to
the ATAC rows alone.

## 3. Why I redid QC and added doublet detection

Once the peaks were removed, the old nCount/nFeature thresholds no longer
measured the same thing. I re-examined each sample's distributions, quantiles
and scatter plots. The tails are continuous, with no single obvious cut-off.
Control samples have lower RNA complexity, and control_rep1 has a higher
mitochondrial baseline. I therefore used within-sample MAD thresholds with a
floor of 500, instead of reusing the old upper limits. The rules and
sensitivity checks are in `decision_log.md`.

QC removed 430 of 24,474 input nuclei. scDblFinder, run per sample/capture,
then flagged 2,494 predicted doublets, leaving 21,550 nuclei. These are
algorithmic predictions; the loading-rate prior and missed homotypic doublets
remain uncertain.

The name GEO `RAW.tar` does not mean it contains raw droplets. After checking
the local files I found no empty-droplet matrix, so I did not force ambient
RNA correction.

## 4. Why I kept both unintegrated and Harmony embeddings

My earlier reasoning, "this is a stress study, so it must not be integrated",
was too absolute. But "Harmony is required" is not established either. The
same gene-only singlets go into two embedding branches, and downstream DE uses
the raw RNA counts. Harmony improves replicate mixing, but canonical identity
consistency does not improve with it, and condition-associated structure is
weakened. So the unintegrated branch is primary and Harmony is a sensitivity
comparison. This does not prove that every Harmony change is over-correction.

The PC grid is locally stable around 20–30 PCs, and I used 25 PCs for the
comparisons that follow. A second stable region at 40–50 PCs gives a partition
clearly different from the 25-PC one, so I do not claim stability across the whole range. The
new 18 and 17 clusters are run outputs, not targets set to match the old
numbers.

## 5. Why I re-annotated instead of reusing old cluster names

The new partition was checked again against Shahan 2022 Data S2 signatures,
canonical markers and dataset-specific markers. The original literature
citations and spatial domains are kept in the reference tables. Evidence from
the same source is not fully independent validation. I kept broad lineages
where possible: Root_cap is not forced into LRC versus columella, and the
meristematic and some initial-like groups stay low confidence.
Dividing_meristematic is a cell state, not a single lineage.

Many cells from the old LRC_like group now fall into Atrichoblast and Cortex,
so the old LRC-specific interpretation is withdrawn. Root_cap has strong cell
continuity with the old Columella group, but that is not enough to restore
every fine label. The manual mapping is tied to the MD5 of the reviewed
partition, so any future clustering change requires a new review of the
annotation evidence.

## 6. How DE matches biological replicates

Raw gene counts are summed by population × biological sample, with at least 50
cells required per sample. The formal gene background is the 32,833 genes
shared by all four samples. DESeq2 uses `~condition` and does not invent a
pairing. DE is run separately in both branches, and differences caused by
changed membership are reported as sensitivity.

The eligible populations are the same at 25, 50 and 100 cell thresholds. This
only shows that eligibility is stable in this run. Two biological replicates
per condition still limit inference, and low-confidence labels cannot support
fine cell-type conclusions.

## 7. Why the Project 1 contextualization adds background tests

Overlap between the 137 interaction candidates and stress DEGs does not by
itself show a cell-specific interaction. Fisher tests with BH correction show
that the candidates are enriched in every eligible population relative to the
ordinary stress-DEG background. When the comparison set is the Project 1
non-interaction drought DEGs instead, no population in either branch reaches
BH-significant extra enrichment. I therefore keep the link between the Project
1 and Project 3 stress contexts, but this does not show Pericycle specificity
or validate the BRL3 interaction.

Quote characters in a Project 1 TSV caused a parsing warning. After fixing the
quote/comment options I reran the contextualization, and the earlier run with
the warning is superseded. This fix is recorded in the decision log rather than
removed from the history.

## 8. What the seven GO terms can show

The seven GO rows contain only 13 distinct candidate genes, and the gene sets
overlap. In both branches only three terms have a non-zero DEG overlap. This is
context for an existing candidate set, not a genome-wide Project 3 GO
enrichment. Zero overlap does not show that a biological program is absent,
and shared genes are not independent lines of evidence.

## 9. Where to read next

Start here, then read `decision_log.md` for the reasoning behind parameters,
and use `comparison_report.md` with the current figures and tables to check
the numbers. `version_comparison.md` shows how far the old conclusions still
hold, and `manual_review.md` gives a review order. New evidence can change the
current choices, but the reasons and version differences should be recorded
rather than erasing earlier judgments.
