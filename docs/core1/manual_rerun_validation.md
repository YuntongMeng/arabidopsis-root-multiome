# Core 1 manual rerun validation — 2026-10-06

The user restarted the R session and invoked stages 01–07 from the original processed input matrices in RStudio, inspecting outputs with the assistant. This record summarizes the observed console tables; it is not a transcript, a claim of row-by-row identity, or a proof of biological truth.

| Checkpoint | Observed result |
| --- | --- |
| Input | Four samples; 32,833 genes each; 0 coordinate IDs after filtering |
| QC | 24,474 input nuclei; 430 removed by QC |
| Doublets | 2,494 predicted doublets removed; 21,550 singlets retained |
| Embeddings | Seven dimension choices each for PCA/Harmony; both branches retained |
| Evidence | PCA/Harmony same-condition neighbor fractions 0.8589404/0.5997309; atlas consistency 0.4852189/0.4846682 |
| Pseudobulk | Nine PCA and ten Harmony eligible populations at >=50 nuclei in every sample |
| Structural checks | All eight checks TRUE, including gene-only DE universe and no peak markers |

The user-pasted input, QC, branch and DE summaries agree with the preceding corrected run. This verifies those summaries and the workflow execution, not every gene-level numerical value independently. The eight checks are recorded in `results/core1/tables/structural_invariants.csv`; their criteria concern data structure, not fixed biological counts.

Observed messages were nonfatal: Harmony's BLAS message indicates a single-thread fallback; Seurat's UMAP message describes its R UWOT default; AnnotationDbi's 1:many message flags annotation mapping multiplicity. Marker provenance and low-confidence groups remain reviewable. Completion does not establish annotation correctness or an optimal PC count.

After this manual rerun, portability and presentation were cleaned up. Figure code was executed again using the saved analysis objects, generating 22 main and two supplementary figures. This final rendering did not rerun clustering, marker tests or DE. Personal paths now reside only in untracked local settings; archived text retains its historical paths.
