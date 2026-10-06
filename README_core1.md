# Project 3 — Core 1: Arabidopsis root-tip snRNA-seq under osmotic stress

状态：Core 1 全流程已由用户从原始矩阵重跑完成，关键结果与原归档一致，用户已授权提交 GitHub。最初整理仅导出已有对象；随后修复 Seurat dimensions 返回 double 的 vapply 类型问题，用户运行了完整流程。当前正式重跑图表位于 `results/core1/`，最初归档位于 `results/core1_frozen/`。

## 文件入口

| 文件 | 内容 |
| --- | --- |
| `scripts/01_core1_rna_preprocessing_clustering.R` | 矩阵审计、分布驱动 QC、PC robustness、最终聚类、样本组成检查 |
| `scripts/02_core1_cell_population_annotation.R` | JoinLayers 修复、统一 strict markers、signature 分数、人工标签与验证 |
| `scripts/03_core1_pseudobulk_stress_DE.R` | population×sample 汇总、eligibility、共有 feature universe、正式 DE |
| `scripts/04_core1_P1_contextualization.R` | 137 个 P1 candidates 和已有 7 个 GO terms 的 P3 context |
| `scripts/00_core1_common.R` | 共享文件处理与输出辅助函数 |
| `docs/core1_analysis_summary.md` | 思考过程、方法调整、完整数值汇总和解释边界 |
| `docs/core1_manual_review.md` | 人工检查清单与已知证据缺口 |
| `results/core1_frozen/` | 本次从已完成会话导出的原结果；整理脚本不覆盖此目录 |
| `data/core1_reference/` | 固定的 annotation/signature 和 GO symbol→TAIR 输入 |

原项目目录为 `/Users/tt/Documents/arabidopsis_root_multiome`。该目录已经包含四个 RNA 输入样本、P1 输入、2026-09-30 Seurat checkpoint、PC 比较和 strict marker 表。整理包不复制大体积原始矩阵和原 Seurat checkpoint；完整复现需要原项目中的这些文件。

## 人工检查时先阅读文件

所有脚本只定义函数；`source()` 不会启动生物学分析。函数调用才会执行相应阶段。非 base R 函数均使用 `package::function`，包括 `dplyr::select`、`dplyr::filter`、`tidyr::pivot_longer`、`AnnotationDbi::mapIds` 和 `DESeq2::results`。

此次归档的 `core1_existing_results.rds` 是轻量结果 checkpoint，含已有 pseudobulk counts、metadata、feature universe、13 个 DE 结果和 annotation 输入；它不是完整 Seurat 对象。已有细胞 metadata 包括 signature 分数和最终标签。原始完整 Seurat checkpoint 仍位于原项目的 `results/reproducibility/p3_core1_checkpoint_2026-09-30.rds`。

## 将来明确要求复现时的运行方式

从原项目根目录运行，先载入辅助函数及四份主脚本：

```r
source("scripts/00_core1_common.R")
source("scripts/01_core1_rna_preprocessing_clustering.R")
source("scripts/02_core1_cell_population_annotation.R")
source("scripts/03_core1_pseudobulk_stress_DE.R")
source("scripts/04_core1_P1_contextualization.R")
```

读取归档结果的两个入口（不做新的 DE 检验）：

```r
# 导出已完成的正式 DE 到独立的 results/core1_reproduction 目录。
# run_core1_03(recompute = FALSE)
# 基于原正式 DE 表重新汇总 contextualization；不会重新检验。
# run_core1_04()
```

下面的完整运行会重新分析，只在用户明确要求复现时执行；本次没有执行：

```r
# run_core1_01()
# run_core1_02(
#   checkpoint = "results/core1_reproduction/checkpoints/01_clustered.rds",
#   recompute_evidence = TRUE
# )
# run_core1_03(recompute = TRUE)
# run_core1_04(de_file = "results/core1_reproduction/tables/pseudobulk_DE_all.csv")
```

也可以从 2026-09-30 checkpoint 恢复 annotation：`run_core1_02()` 默认复用保存的 strict markers 和每个细胞已有的 signature 分数，再附加 reviewed labels，不重新算 markers/scores。该入口会生成验证图和新的 annotated checkpoint。它会先按细胞 ID 核对历史 cluster；如果重跑后 partition 或 cluster 编号改变，脚本停止，要求重新检查标签，防止把历史 annotation 贴到不对应的 cluster。

运行时默认将新输出写入 `results/core1_reproduction/`。如指定其他 `output`，后续阶段须同时传入相应的 `checkpoint`/`de_file` 路径。输出目录可包含图、表和 stage checkpoints；完整重跑生成的关键图已收录在 `results/core1/figures/`。

## 环境与证据

本次导出会话：R 4.6.1、Seurat 5.5.1、SeuratObject 5.4.0、DESeq2 1.52.0、org.At.tair.db 3.22.0、dplyr 1.2.1、tidyr 1.3.2。完整会话记录在 `results/core1_frozen/sessionInfo_export.txt`。这是 2026-10-05 导出环境；2026-09-30 最初运行的完整 sessionInfo 未保存，不能声称已经锁定了当时所有依赖。脚本显式写出了历史使用的参数及可核查的默认随机种子。

方法来源是原项目 `.Rhistory`、2026-10-05 RStudio history、已有磁盘结果和直接导出的最终 RStudio 对象。只有原聊天最近五轮可读取；不足的早期上下文由本地记录补充。没有填造不在证据中的结果或参考文献。

`docs/core1_integrity_report.md` 记录了本次语法与结果一致性检查。另见 `docs/core1_rerun_verification.md`：完整重跑已通过，所列关键汇总与原归档一致；未声称每一行统计量均已逐项比对。
