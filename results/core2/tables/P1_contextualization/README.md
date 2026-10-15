# Project 1 candidate contextualization in Core 2

This describes genomic proximity to the original Project 1 interaction candidates.
Mapping rule: overlap the gene body OR the strand-aware TSS window (upstream 2000 bp; downstream 500 bp).
The two overlap flags are retained. A peak may match multiple candidate genes.
This proximity window is separate from the earlier GeneActivity 500-bp upstream window.
Peaks, counts, models, p values, and primary selection thresholds are unchanged.
Counts in P1_summary refer to candidate genes with at least one matching peak.
with_FDR_peak does not mean the gene itself has a gene-level adjusted p value.
Multiple nearby peaks are not treated as independent gene-level tests.
P1 interaction log2FC is not the same contrast as Core 2 stress/control log2FC.
Matching signs cannot by themselves establish cross-project agreement.
No_common_peak_match does not prove a gene is inactive or chromatin is closed.
Modelled_peaks_without_FDR_support does not prove absence of a stress response.
Reference_annotation_missing and Outside_nuclear_analysis are retained explicitly.
Root_cap remains exploratory because of weak within-condition reproducibility.
No GO analysis, target-gene validation, or new enrichment/significance test is performed.
