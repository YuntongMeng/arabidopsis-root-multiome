# Core 1 exploration and supporting diagnostics

Exploration records alternatives and uncertainty rather than defining another primary version.

| File or evidence | Role |
|---|---|
| `09_render_figures.R` | Re-renders all 22 figures from saved outputs; adds ARI, sample-count and P1 control views |
| `08_extra_diagnostics.R` | Old/new HVG audit, stress-marker detection, sample-wise atlas scores and raw-count preservation |
| `extract_atlas.py` | Re-extracts published panels from the locally stored Shahan workbook; no download |
| `../../scripts/core1/03_embeddings.R` | Executable paired embedding and PC-grid sensitivity, required by the current workflow |
| `../../scripts/core1/04_branch_evidence.R` | Comparative marker and neighborhood evidence before choosing the primary branch |
| `../../results/core1/tables/` | Original quantitative diagnostic results, kept with the current analysis |
| `../../results/core1/figures/` | QC distributions, elbow plots, UMAP and canonical marker evidence |
| `../../results/core1/diagnostics/` | Feature universes, session records and supporting audit files |

Run additional checks through the main runner's `diagnostics` stage. Source files
define functions; they do not independently start biological analysis.
The embedding/evidence stages stay in scripts/core1 because the current workflow
needs them. Results are not duplicated here; a single result location avoids divergence.
Method selection and limitations are documented in docs/core1/decision_log.md.
