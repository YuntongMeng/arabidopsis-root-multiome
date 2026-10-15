# Core 2: what the analysis supports

## Scientific conclusion

The ATAC reanalysis provides limited condition-associated evidence. One Atrichoblast-lineage peak passes the current within-population FDR criterion, while the other four evaluated populations yield no significant regions. None of the 137 Project 1 interaction candidates has a nearby matched peak that passes this criterion. These observations describe the available evidence; they do not establish that stress has no chromatin effects or that P1 candidates lack stress functions.

The purpose was also to learn and document the complete ATAC analysis chain. We established common measurement intervals, reconstructed fragments, calculated fragment-level QC, matched nuclei to corrected RNA calls, inspected embeddings and population identity, and compared biological-sample pseudobulks. Completing that chain does not guarantee strong biological discoveries.

## Why the original matrices could not simply be merged

The four published peak lists contained 19,858 / 24,317 / 24,490 / 25,001 regions. Exact coordinate-string intersection was zero, although many intervals overlapped. We used union/reduce of those existing intervals and recounted fragments in 27,677 common nuclear regions. Expanding or renaming matrix rows alone would not recover fragments outside the original regions. No new peak calling was performed.

The RNA assay uses only Gene Expression rows. The ATAC assay uses only peak-region counts. This preserves the correction established in Core 1 and prevents mixed gene/peak features from being interpreted as RNA expression.

## QC, pairing and population interpretation

| Sample | Initial | Basic ATAC QC | Paired RNA singlets |
| --- | ---: | ---: | ---: |
| control_rep1 | 3,431 | 1,417 | 1,239 |
| control_rep2 | 6,631 | 4,508 | 3,817 |
| stress_rep1 | 7,827 | 5,772 | 4,841 |
| stress_rep2 | 6,585 | 4,925 | 4,249 |

The 14,146 retained nuclei carry capture identifiers and raw barcodes. Matching uses sample + barcode, because a raw barcode can recur in different captures. RNA-derived predicted doublet calls were reused from the corrected gene-only Core 1. The 196 ATAC-basic-pass nuclei without a tested RNA doublet call were excluded from paired analysis, not relabeled as doublets.

Five ATAC clusters were observed at the chosen working resolution; this number is not a target or the number of RNA cell populations. The descriptions below summarize partial evidence:

| ATAC cluster | RNA composition and provisional description |
| --- | --- |
| 0 | Atrichoblast 54.2%, Cortex 19.5%; Atrichoblast-dominant mixed |
| 1 | Endodermis 69.5%; Endodermis-dominant |
| 2 | Pericycle 66.8%; Pericycle-dominant |
| 3 | Meristematic unresolved 35.8%, dividing meristematic 26.3%, epidermal initial-like 20.3%; mixed state-related group |
| 4 | Root cap 27.3%, Cortex 20.7%, Atrichoblast 19.0%; mixed |

Marker-region ATAC summaries and SCR coverage provide partial supporting checks, not independent ground truth. An ATAC cluster's description does not replace each nucleus's RNA population. Differential testing uses the RNA-defined populations rather than pretending the five ATAC clusters are pure cell types. The same UMAP positions can be colored by sample, condition, cluster or RNA identity without recomputing the coordinates.

## Reproducibility depends on population

Within-population pseudobulk VST correlations were:

| Population | Control replicate Pearson r | Stress replicate Pearson r |
| --- | ---: | ---: |
| Atrichoblast_lineage | 0.948 | 0.972 |
| Cortex_lineage | 0.732 | 0.930 |
| Endodermis_lineage | 0.735 | 0.847 |
| Pericycle | 0.762 | 0.703 |
| Root_cap | 0.370 | 0.191 |

Control 1 is not uniformly discordant across every population, and weak reproducibility is not restricted to Control 1. Root_cap is explicitly exploratory. Its small testable peak set and weak repeatability make the zero-significant result particularly uninformative about absence of biology. A visually separated four-point PCA is not by itself proof of a treatment effect.

## The one differential-accessibility candidate

In Atrichoblast_lineage, `4-5457000-5457860` has estimated log2FC 1.860 and adjusted p 0.004526. Its DESeq2-normalized counts are 25.5 / 21.2 / 83.5 / 83.5 in the four samples. Both stress replicates are above both controls; this region is not driven solely by an unusually low Control 1 value.

AT4G08570 is the nearest annotated gene at a minimum TSS distance of 900 bp. AT4G08580 and PUM15 also appear in the nearby-region table. The latter's TSS is 6,566 bp away even though part of its gene body enters the displayed ±5 kb neighborhood. Distance, gene-body overlap and TSS windows are distinct definitions. No regulatory target or causal effect on expression has been established.

This remains a candidate from five population-wise analyses, not an experiment-wide validated discovery. No significance claim is made from a point-estimate fold-change alone.

## P1 comparison

All 137 P1 candidates occur in the nuclear reference annotation; 120 have one or more common peaks overlapping the gene body or the strand-aware TSS window (2 kb upstream / 500 bp downstream). The other 17 have no matching region under this rule, not missing annotation.

| Population | P1 genes with common peak | With ≥1 peak in model | With ≥1 FDR-significant peak |
| --- | ---: | ---: | ---: |
| Atrichoblast_lineage | 120 | 120 | 0 |
| Cortex_lineage | 120 | 120 | 0 |
| Endodermis_lineage | 120 | 118 | 0 |
| Pericycle | 120 | 119 | 0 |
| Root_cap | 120 | 19 | 0 |

The unit in this summary is a candidate gene, not a peak. A nearby peak p value is not a gene-level p value. P1 tests BRL3 genotype × drought interaction; P3 compares wild-type sorbitol stress with control. Matching effect signs would not validate the same contrast. Core 1's RNA context remains valid within its own limitations; this ATAC comparison does not overturn it.

## Where this iteration stops

We retain the single candidate, the sparse/null comparison results and the uncertainty. No threshold relaxation was used to produce a longer candidate list. No new ATAC GO analysis was performed. Same-nucleus RNA–ATAC comparison is already part of Core 2; The user cancelled Core 3 on 2026-10-15 and concluded the project at Core 2; see the [scope decision](../project_scope.md). Further mechanistic analysis would require its own question and validation, not merely another stage name.

[Methods](methods.md) · [Control 1](control1_review.md) · [Decisions](decision_log.md) · [Validation](validation.md)
