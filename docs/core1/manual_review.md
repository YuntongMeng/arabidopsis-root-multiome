# Core 1 manual review order

1. Use the README and `analysis_story.md` to confirm the current version, the research question and the decision process.
2. Check `results/core1/tables/input_matrix_audit.csv`: feature types and coordinate IDs.
3. Look at the QC figures, thresholds, doublet calls and the assumptions in `decision_log.md`. Similar retention numbers are not a substitute for checking.
4. Look at PCA/UMAP, PC ARI, replicate mixing and condition preservation, keeping the primary analysis separate from the Harmony sensitivity branch.
5. Compare the canonical DotPlot, strict markers, atlas scores and `cluster_annotation_review.csv`, especially for low-confidence groups.
6. Check population × sample counts, eligibility, raw-count preservation and the DE tables.
7. Check the Project 1 Fisher background definition, the BH test family, the drought control and the GO gene-set overlap.
8. Use `comparison_report.md` to check which old conclusions are kept or withdrawn.

I reran stages 01–07; the observed summaries matched the previous corrected
run, and all eight structural checks passed. A short record is in
[manual_rerun_validation.md](manual_rerun_validation.md). Cell identities and
biological interpretation still have clear limits; a completed run does not
confirm the scientific conclusions.
