# Control 1: evidence, explanation and remaining uncertainty

## What is explained

The unusually low unfiltered Control 1 median and retention rate are associated with a large subgroup of low nuclear ATAC barcodes. The same subgroup is already low in the original published peak matrix. This finding predates our common-peak reconstruction and argues against the checked recounting steps creating the subgroup.

| Control 1 descriptive group | Barcodes | Published ATAC median | Nuclear fragment median | Organelle median % |
| --- | ---: | ---: | ---: | ---: |
| Nuclear fragments <500 | 1,822 | 233 | 70 | 90.69 |
| Nuclear fragments ≥500 | 1,609 | 5,722 | 3,570 | 6.98 |

The 500 cutoff here is a diagnostic grouping, not our final 1,000-fragment filter. The low group exceeds half the original 3,431 nuclei and pulls the unfiltered median down. Its organelle fraction is computed before generating nuclear-only fragment copies.

Original ATAC count medians were 797 / 3,437 / 3,820 / 3,979 across C1/C2/S1/S2. Rank correlations between original ATAC counts and common-peak recounts were .9920 / .9984 / .9987 / .9989. Different peak sets/count units need not have identical numerical counts; the diagnostic asks whether the same barcodes stay relatively low or high.

## Independent audit

A separate Python implementation scanned all **49,108,250 original fragment rows**, without using the original Signac counting helper. It recomputed fragment totals, organelle fraction, nuclear FRiP, nucleosome signal and peak-region counts. All eight checked metrics agreed (only floating-point differences up to approximately 5.7e-14), and the full 27,677 ×3,431 matrix agreed, including 3,884,627 nonzero entries and the remaining zeros. The nuclear fragment copy and recorded file hashes also agreed.

See the [original Chinese audit record](control1_audit_original_zh.md) and [audit evidence](../../results/core2/diagnostics/control1_audit/). That audit record describes the stage when filtering/doublet selection was still pending; the later paired selection is documented in Core 2 results.

## What is not explained

The audit did not rerun TSS enrichment or remap FASTQ. It does not rule out every upstream assumption or identify which experimental preparation, capture, library or sequencing factor produced the subgroup. Organelle-rich low nuclear signal is an observed association, not a demonstrated causal laboratory diagnosis.

After filtering, retained Control 1 nuclei have higher per-nucleus fragment medians in many populations, but fewer nuclei contribute to pseudobulk. This does not prove higher absolute accessibility or best overall quality. Endodermis has C1/C2/S1/S2 nuclear counts 109/473/575/528 and median nuclear fragments 3715/2295/2045/2057. DESeq2 scaling compensates gross library scale, not sampling noise, subpopulation composition or every technical difference.

Atrichoblast control-replicate correlation is .948, while Endodermis is .735; Root_cap and Pericycle also show weak stress-replicate agreement. It is inaccurate to attribute all instability to Control 1. Control 1 was retained, not silently dropped to leave a single control replicate. The sole significant Atrichoblast region has both stress values above both controls.

The defensible explanation is therefore: **a pre-existing low-nuclear-signal subgroup explains much of the original Control 1 QC difference; checked recounting errors were not found; specific experimental causes and residual population-dependent differences remain unresolved.**
