# Core 2 methods and fixed working parameters

## Data and reconstruction

Processed mixed-feature matrices were restored from GEO GSE235495. Capture/run mapping is control_rep1/SRR24921616, control_rep2/SRR24921615, stress_rep1/SRR24921614 and stress_rep2/SRR24921613. Gene Expression rows supply Core 1; Peaks rows supply the original interval definitions. Fragments were unavailable among the locally restored processed inputs, so public SRA technical and biological read streams were recovered.

The four FASTQ streams have lengths 8 / 24 / 50 / 50 bp. Barcode extraction uses the reverse complement of I2 bases 9–24 (one-based), with paired ARC ATAC-to-GEX whitelist translation. Preview and full-stream synchronization checks are recorded. Chromap 0.3.2-r518 uses the ATAC preset, MAPQ ≥30, maximum insert length 2,000, barcode Hamming distance ≤1 and correction probability ≥0.9. Tn5 shift is applied once by that preset (+4 start / −5 end); cell-level duplicates are collapsed. Original fragment files include nuclear contigs 1–5, Mt and Pt for competitive alignment.

Reference: Ensembl Plants release 62 TAIR10 toplevel FASTA, including Mt/Pt. FASTA/reference hashes, contig lengths, read checks and tool versions are in [upstream provenance](../../results/core2/diagnostics/upstream/). This is an ATAC-only Chromap reconstruction, **not an identical Cell Ranger ARC rerun**.

## Common regions and count units

Original Peaks intervals are interpreted as zero-based half-open BED. Import to GRanges adds 1 to start only. `GenomicRanges::reduce` merges the union, followed by 20 < width <10,000 bp. All 27,677 reduced regions pass; widths range 154–2,134 bp, median 884 bp. The common BED is [stored here](../../data/core2_reference/common_peaks.bed).

One deduplicated fragment row contributes once to every overlapping peak. Column 5 PCR support is not a weight. A fragment overlapping multiple peaks can contribute to multiple matrix rows. FRiP instead counts each fragment at most once in its numerator. Our fragment-region matrix is not required to equal original Cell Ranger cut-site counts.

## QC and pairing

| Quantity | Definition / use |
| --- | --- |
| nCount_ATAC | Sum of common-peak counts per nucleus |
| nFeature_ATAC | Number of common peaks with nonzero count |
| nuclear_fragments | Unique deduplicated rows on chr1–5; selection ≥1,000 |
| nucleosome_signal | Signac mononucleosomal / short-fragment ratio; selection <2 |
| TSS.enrichment | Shared nuclear protein-coding gene TSS set from TAIR10.62 GTF; selection >1 |
| FRiP_nuclear | Unique nuclear fragments hitting ≥1 common peak / all nuclear fragments; diagnostic |
| percent_organelle | 100 × (Mt + Pt fragment rows) / all fragment rows in the original file; diagnostic |

Nuclear-only fragment copies are indexed for nuclear QC; original fragments retain organelle information. TSS uses unique strand-aware gene starts, `n = NULL`, `fast = FALSE`, region extension 1,000 bp. The recorded Signac NucleosomeSignal implementation uses lengths <147 for short fragments and 147<length<294 for mono-nucleosomal fragments; exact boundary values are excluded. NucleosomeSignal is deprecated in this installed Signac version; the recorded run is retained rather than silently changed to another QC implementation.

FRiP and organelle percentage have no additional hard filter. The union-defined peak universe affects FRiP interpretation; a high ratio is not a test that peak calling is correct. No arbitrary upper count threshold was used. Corrected, gene-only Core 1 scDblFinder calls are matched by capture + barcode; retain RNA QC-pass predicted singlets, and exclude not-tested nuclei from paired inference. This is not an independent ATAC doublet call.

The paper's nCount_ATAC threshold and our nuclear-fragment threshold are different quantities. The working TSS and nucleosome cutoffs are informed by the original study, not a literal reproduction of all its QC rules. The original study used Cell Ranger ARC 2.0.2; our reconstruction, count unit and downstream choices differ. [Liu et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11199978/)

## Embedding, matching and marker checks

- TF-IDF method 1, scale factor 10,000; FindTopFeatures min.cutoff=0.
- SVD requested 50 LSI components. LSI 2–30 is the working set (29 components); LSI 1 was strongly associated with depth. Checking the first 30 does not mean only 30 were calculated. Moderate depth associations in other dimensions remain a limitation.
- UMAP: uwot, cosine, neighbors 30, min.dist 0.3, seed 42; coordinates named umap.atac.
- Neighbors: k=20 on the same LSI dimensions. Clustering: graph atac_snn, algorithm 1, resolution 0.5, seed 42. Five clusters were observed, not required.
- No Harmony, anchor transfer, WNN or reclustering to force RNA agreement was used in this ATAC branch.
- RNA population/confidence comes from corrected Core 1's primary unintegrated PCA branch. The 59 unique canonical marker genes have row-wise provenance in `data/core1_corrected_reference/canonical_markers_provenance.csv`, extracted from Shahan 2022 Data S2. They are not independent ground truth.
- GeneActivity uses nuclear protein-coding transcripts with 500 bp upstream and 0 downstream extension, gene.id=TRUE. GENE_ACTIVITY is LogNormalize-normalized at 10,000. The dot plot uses scale=FALSE: size is percentage with nonzero region count; color is average normalized ATAC signal, not RNA expression.
- Coverage uses transcript + exon annotation (including tx_id), avoiding the earlier transcript-only plotting error. SCR transcripts span chr3:20,069,987–20,072,878, with ±2 kb displayed. The x axis is genomic position in base pairs. Group-normalized coverage cannot eliminate low-nucleus sampling noise.

## Pseudobulk and differences

Sum raw ATAC counts by sample × RNA population. Forty-six observed combinations exist; missing stress groups are not fabricated as zero columns. Metadata are explicitly aligned by unique cell identifiers; internal name-cache attributes are removed from identifiers without dropping cells or changing barcode text.

The five working populations are Atrichoblast, Cortex, Endodermis, Pericycle and Root_cap. All have ≥50 nuclei in every sample. This is a limited panel of clearer lineage labels selected before inspecting differential results, not every eligible RNA label. Dividing_meristematic is a state, initial/unresolved labels remain uncertain, Trichoblast has only 45 Control 1 nuclei, and Phloem has 5. The panel does not establish absence of effects in omitted populations.

For each population, retain peaks with raw count ≥10 in ≥2 samples. DESeq2 uses design ~condition, control baseline and sfType=poscounts. VST blind=TRUE is for diagnostics only, not the input to inference. Sample PCA uses the 500 highest-variance retained peaks; sample correlation uses all retained peaks. Automatic local dispersion fitting where parametric fitting failed is recorded, not treated as a fatal error.

`DESeq2::results` uses contrast stress/control, alpha=.05, BH and independent filtering. FDR <.05 is the statistical criterion within each population. |Estimated log2FC|≥1 is additional prioritization, not testing a null at a twofold boundary. Five separate result sets do not imply experiment-wide FDR control. Descriptive replicate direction is recorded; it is not another significance test. Two biological replicates per condition remain the inference unit.

## P1 positions and scope

All 137 interaction candidates are taken from `data/project1/FULL_interaction_DEGs_padj0.05_log2FC1.tsv`. Match common peaks to gene body OR strand-aware TSS upstream 2,000 / downstream 500 bp. Keep both overlap flags. A peak may match multiple genes. This mapping window differs from the 500 bp GeneActivity extension. Existing differential statistics are reused; there is no new gene-level test or enrichment test.

No new GO analysis, motif analysis, co-accessibility, peak-to-gene links, expression-accessibility correlation model or causal regulatory network was performed. Same-nucleus identity matching and ATAC checks are the completed joint comparison; deeper joint modeling would be additional work.

## Local supplemental methods, 2026-10-15

See [supplemental_review.md](supplemental_review.md) for the predeclared counts>=5 filter, ten model-diagnostic figures, numerical initial-dispersion boundary proportions, and RNA-only selection of 15 response examples. Primary counts>=10 inference, sample replication and QC are unchanged. No new RNA model or Control 1 attribution experiment was run. Publication of these reviewed additions was authorized on 2026-10-18.
