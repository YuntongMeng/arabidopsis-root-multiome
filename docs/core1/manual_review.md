# Core 1 人工审阅顺序

1. 从 README 和 analysis_story.md 确认当前版本、研究问题和决策过程。
2. 核对 results/core1/tables/input_matrix_audit.csv：feature type 与 coordinate IDs。
3. 看 QC figures、阈值、doublet calls 和 decision_log.md 的假设；不要以保留数量相似代替检查。
4. 看 PCA/UMAP、PC ARI、replicate mixing 和 condition preservation；区分主分析与 Harmony sensitivity。
5. 对照 canonical DotPlot、strict markers、atlas scores 和 cluster_annotation_review.csv，尤其低置信度 groups。
6. 查看 population × sample counts、eligibility、raw-count preservation 和 DE tables。
7. 查看 P1 Fisher 背景定义、BH test family、drought control 和 GO gene-set overlap。
8. 用 comparison_report.md 检查旧结论的保留/撤回边界。

用户已实际运行 01–07，观察到的 summary 与此前修正版一致，八项结构检查全通过。简短记录见 [manual_rerun_validation.md](manual_rerun_validation.md)。细胞身份与生物学解释仍有明确限制；运行完成不等于科学结论已经证实。
