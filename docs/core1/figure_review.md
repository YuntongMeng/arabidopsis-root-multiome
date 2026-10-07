# Core 1 figure review and rendering

`results/core1/figures/` contains 22 main review figures. Two additional ARI plots live under `results/core1/diagnostics/figures/`: the full pairwise heatmap and historical 25-dimension-reference curve. Main review figures do not all need to appear in a report's main text.

The main stability figure compares adjacent dimension choices, avoiding a privileged 25-PC reference. The PCA/Harmony comparison uses coordinate standard deviation. Harmony corrects PCA coordinates and does not re-rank dimensions by explained variance; the comparison therefore is not two conventional explained-variance elbow plots. Keep the PCA elbow, adjacent ARI, marker evidence and mixing/condition diagnostics together when explaining the provisional 25-PC choice. The 40–50 range also has stability and has not been disproven.

Pairwise ARI heatmaps are established diagnostics: [Dune plotARIs](https://hectorrdb.github.io/Dune/reference/plotARIs.html) provides this exact type of view. [OSCA clustering stability](https://bioconductor.org/books/3.20/OSCA.advanced/clustering-redux.html) discusses related ARI-based diagnostics. This supports retaining the full heatmap as supplementary evidence; it does not make ARI a biological accuracy score. Harmony's coordinate SD behavior is documented in its [implementation](https://github.com/immunogenomics/harmony/blob/master/R/RunHarmony.R).

From the project root:

```r
source("scripts/core1/00_common.R")
source("exploration/core1/09_render_figures.R")
run_render_figures()
```

Or run `Rscript scripts/run_core1_corrected.R figures`. This actually executes drawing code from saved tables/checkpoints, without recalculating clustering, UMAP, markers or DE. Full `all` recomputes stages 01–07 and then draws these figures. RStudio Plots need not display every image saved by ggsave; inspect the PNG files directly.

Canonical DotPlots require zoom. Low-confidence and state labels remain visible in population counts. P1 ordinary-background enrichment must be read alongside the drought-DEG control; the GO plot is contextual overlap, not new enrichment.
