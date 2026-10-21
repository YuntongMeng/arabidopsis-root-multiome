# Core 1 version comparison and conclusion status

| Item | Initial version | Current corrected version | Note |
|---|---|---|---|
| Input | Genes + peaks mixed into RNA | Gene Expression only | The initial input error is confirmed |
| QC | Mixed-count thresholds | Thresholds from gene-only distributions | The two sets of metrics are not directly equivalent |
| Doublets | Not included | scDblFinder per capture | Predictions, not ground truth |
| Retained cells | 24,277 | 21,550 | 21,488 shared cells; several steps changed at once |
| Clustering | Old partition | Unintegrated primary + Harmony sensitivity | Old/new PCA ARI on shared cells 0.639 |
| Annotation | Some fine labels lacked evidence | Sourced broad lineages + uncertainty | Old cluster-ID labels are not inherited |
| Pseudobulk | Old groups | New groups × biological sample | DE counts are not a reproduction target |
| Project 1 overlap | Pericycle 23/118 | PCA 19/120; Harmony 19/117 fitted | Denominator is fitted candidates, not the Fisher background denominator |
| Project 1 enrichment | Mainly overlap | Fisher + BH, plus a drought-DEG control | No evidence of extra enrichment over the control set |
| GO | Seven existing candidate gene sets | Complete zeros, overlap and interpretation limits | Not a new Project 3 GO enrichment |

Still valid: the raw-count pseudobulk approach, the exploratory value of the
Project 1/Project 3 stress context, and continuity for some broader cell
identities.

Must be updated: old QC metrics, cluster membership, annotation, and DE and
context values computed on the old groups.

Not supported as conclusions: the old LRC-specific interpretation, "Pericycle
is uniquely strongest", and any cell-specific validation of the BRL3
interaction.

None of the old HVGs were coordinate peaks, so not all embedding changes can be
attributed to peaks entering the PCA directly.

Full numbers, branch concordance for each population and literature sources
are in the [detailed comparison](comparison_report.md).
