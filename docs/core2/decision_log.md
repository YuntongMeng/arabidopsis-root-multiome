# Core 2 decision trail

| Decision/problem | Resolution | Limit retained |
| --- | --- | --- |
| RNA and peaks coexist in processed matrices | Use corrected gene-only Core 1 and Peaks-only original ATAC intervals | Archived mixed-feature RNA analysis is not current evidence |
| Exact intersection of four peak IDs is zero | Union/reduce overlapping intervals; reconstruct fragments and recount | Common intervals do not require four-way accessibility |
| No fragments in restored processed files | Recover technical SRA streams and reconstruct ATAC with Chromap | Tool substitution is disclosed; not bitwise ARC reproduction |
| Organelle contigs | Include Mt/Pt during alignment; use nuclear copies for nuclear QC and original files for organelle % | Removing contigs from a copy does not remove historical QC metadata |
| QC thresholds | Nuclear fragments ≥1000, TSS >1, nucleosome <2 across captures | No hard FRiP/organelle or arbitrary upper-count cutoff |
| Control 1 low retention | Independent full-fragment/peak-matrix audit and population-level diagnostics | TSS/remapping not independently rerun; laboratory cause unresolved |
| Doublets | Reuse exact-nucleus gene-only RNA scDblFinder calls | Not-tested nuclei excluded from paired set; not called doublets |
| LSI depth | Exclude strongly depth-associated LSI1; use 2–30 as working dimensions | Dimension optimum not demonstrated; some moderate correlations remain |
| Five ATAC clusters vs more RNA populations | Describe RNA compositions and mixed clusters | No forced one-to-one annotation |
| Transcript-only coverage error | Include transcript and exon records with tx_id | Coverage remains descriptive, not target validation |
| Crowded 59-marker plot | Keep as internal diagnostic; use smaller explanatory panel | No requirement to make all reference markers cluster-specific |
| Name alignment error | Strip name-cache attributes, then assert exact unique IDs and order | Never silently intersect or remove capture prefixes |
| Biological inference | Sample × RNA population raw-count pseudobulk, two replicates per condition | Nuclei cannot replace biological replication |
| Different sample PCA geometry | Independent population PCAs, corrected equal x/y unit scale | Panels are not a shared coordinate system |
| Very few significant peaks | Retain FDR .05 and estimated twofold priority rule | No discovery-count-driven threshold relaxation |
| Root_cap weak repeatability | Preserve as explicitly exploratory | Zero results are particularly weak evidence of absence |
| One Atrichoblast candidate | Review all four sample counts and nearby genes | Proximity and direction consistency do not establish regulation |
| P1 cross-project comparison | Fixed gene-body / TSS-nearby rule for all137 candidates | Different contrasts; no new gene-level or enrichment test |
| GO and deeper regulatory analysis | Not performed in Core 2 | Not necessary to force a larger story from one candidate |
| Core 3 | Initially deferred; user cancelled the plan on 2026-10-15 and concluded Project 3 at Core 2 | Joint identity comparison already belongs to Core 2; no separate Core 3 deliverable |
| Submission | Local preparation and user downstream validation completed; user explicitly authorized GitHub publication on 2026-10-15 | Upload code/results record only; raw inputs remain local; fresh FASTQ/QC rerun not claimed |

The recorded models/results come from the user's RStudio execution. Repository preparation reorganizes those outputs, reconstructs previously interactive steps into scripts, and redraws saved PCA coordinates. New code is not retrospective proof that a complete fresh run reproduces every result.


## Subsequent user validation

The user manually reran the assembled downstream pipeline in RStudio from the four existing unfiltered QC objects. All ten comparisons matched (statistics tolerance 1e-7). This supersedes the earlier preparation-time statement that no fresh downstream run had occurred. It does not assert fresh upstream/QC reconstruction or scientific ground truth. Original results remain unchanged.

## Final closure, 2026-10-15

The user cancelled the separate Core 3 plan and authorized GitHub publication. Project 3 now concludes with Core 1 and Core 2. No further biological analysis or new correction algorithm was run during closure. See [final scope](../project_scope.md).
