# Core 1 analysis summary / decision trail

## 研究问题与结果定位

Project 3 使用 WT Arabidopsis root-tip snRNA-seq，对照与 sorbitol osmotic stress 各两个 biological samples。Core 1 从 RNA QC 和聚类到保守 cell-population annotation，再以 biological sample 为统计单位做 pseudobulk DE，并将 Project 1 的 BRL3×drought interaction candidates 和 enriched GO annotations 放入 P3 的细胞群体 context。P1/P3 的基因型、组织采样和处理条件不同，结果属于 contextualization，不能称为 P1 的 validation。

最终保留 **24,277 nuclei、18 transcriptional clusters、17 annotation populations**；其中 **13 populations** 满足每个样本至少 25 nuclei 的正式 DE 标准。18 clusters 与 17 populations 的差异来自 clusters 0/5 同标为 LRC_like。

## 1. 输入矩阵审计和分布驱动的 QC

分别读取四个 RNA matrix/barcode/feature 文件，使用 feature 第一列的 ID。审计矩阵维度、feature/barcode 唯一性与 count 合法性，保持 sample、condition、replicate 信息，并在 merge 时给 barcode 添加 sample 前缀。

历史记录先查看 nCount_RNA、nFeature_RNA、percent.mt（ATMG）、percent.cp（ATCG）的分布摘要和两端 quantiles：0.5/1/5%、95/99/99.5%。随后采用统一阈值：**800 ≤ nFeature_RNA ≤ 15,000，nCount_RNA ≤ 50,000**。没有机械应用通用 mt/cp cutoff；它们仍保留用于 QC 和 cluster-level diagnostics。历史未保存所有最初 QC 图与 quantile 输出，脚本保留查看步骤，不能补写未记录的具体分位数。

以下 before 来自原 barcode 文件行数，after 来自已保存的最终 cell metadata；本次未重新过滤矩阵：

| sample | features | before | after | removed |
| --- | --- | --- | --- | --- |
| control_rep1 | 52691 | 3431 | 3388 | 43 |
| control_rep2 | 57150 | 6631 | 6562 | 69 |
| stress_rep1 | 57323 | 7827 | 7774 | 53 |
| stress_rep2 | 57834 | 6585 | 6553 | 32 |

## 2. PCA 选择经历了 robustness comparison

使用 LogNormalize（scale.factor=10,000）、vst（2,000 variable features）、ScaleData 和 PCA。未添加未记录的 regression covariates。ElbowPlot 使 PC20 成为候选，但没有直接确定 PC20；以相同 clustering resolution=0.5 测试 PC10/15/20/25/30，查看 cluster assignments 和全部 pairwise ARI，最终选择 **PC1:25**。

保存的关键 ARI：PC10/25=0.461993，PC15/25=0.653200，PC20/25=0.754838，PC25/30=0.854386。PC25 是结合比较作出的人工决定，不是“找到最大 ARI 即自动选 PC25”的算法。完整矩阵见 `results/reproducibility/pc_ari_matrix.csv`。最终得到 18 clusters。

研究对象包含 stress signal，因此没有自动做 Harmony 或其他 batch correction。先检查 sample/condition UMAP、cluster×sample counts、fractions 和 cluster QC。Cluster 6 counts 为 control_rep1/2、stress_rep1/2 = **223/1096/8/22**；cluster 14 为 **83/587/4/4**，均明显偏 control_rep2。它们没有因此被机械删除；后续以 conservative label 保留，并由每个样本 nuclei 数量的 eligibility 规则限制正式 DE。仅凭这些检查不能证明完全不存在 batch effects。

## 3. Marker 修复、annotation 主证据和验证证据

Seurat v5 分层 assay 初次 FindAllMarkers 提示 data layers 未 join。修复为先 `SeuratObject::JoinLayers(..., assay="RNA")`，再运行 marker test。统一正式标准为 only.pos=TRUE、min.pct=0.25、logfc.threshold=0.25、Wilcoxon；未另设 return.thresh，当前包可核查默认值为 0.01，整理脚本已显式写出。正式表共 **26,059 returned marker rows**，cluster 12 在这一标准下没有 returned marker。

Cluster 12 的 relaxed FindMarkers（only.pos=FALSE、min.pct=0.10、logfc.threshold=0）只是 diagnostic，不能混入正式 marker 表或用于和其他 clusters 作同标准比较。其原 diagnostic 数值未归档；脚本仅在明确重新计算证据时输出到独立 diagnostics 文件。

主方法是 **reference/canonical signature scoring 加人工解释**，不是自动 atlas transfer。最初小 canonical panel 扩展为 v2（8 signatures），查看 cluster mean scores、top two signatures 和 score margin；再以 v3（12 signatures）帮助查看尚未明确的 clusters 0/1/5/8/12/13。Dataset-specific strict markers 和 canonical DotPlot 是 validation evidence。低/混合证据保留 *_like、Phloem_related、Unresolved_transition。

需要人工复核：扩展 v3 的原记录中，AT1G62300、AT1G33280、AT1G62360 等被写作 candidate/placeholder，其他 QC-enriched 项也未附参考来源。整理保留当时 gene lists 和最终标签，未将其改称已验证的 reference atlas。应核查所有 signature 基因及其来源，尤其这些扩展项；无法仅凭当前记录给出完整文献 provenance。

最终 reviewed annotation：

| cluster | annotation | confidence |
| --- | --- | --- |
| 0 | LRC_like | Moderate |
| 1 | Initial_like | Moderate |
| 2 | Trichoblast | Strong |
| 3 | Pericycle | Strong |
| 4 | Endodermis_like | Moderate |
| 5 | LRC_like | Moderate |
| 6 | Phloem_related | Moderate |
| 7 | Atrichoblast | Strong |
| 8 | Epidermal_like | Low |
| 9 | Dividing_meristematic | Strong |
| 10 | Cortex | Moderate |
| 11 | Dividing_S_phase | Strong |
| 12 | Unresolved_transition | Low |
| 13 | Columella | Strong |
| 14 | Atrichoblast_like | Moderate |
| 15 | Endodermis | Strong |
| 16 | Phloem | Strong |
| 17 | Xylem | Very strong |

曾直接把以 cluster 编号为 names 的向量赋给 Seurat metadata，造成 names/barcodes 对齐错误；整理脚本按 cell barcode 显式构建 metadata。它还核对历史 partition，防止重跑后 cluster ID 改变造成错误标签传递。

## 4. Pseudobulk 的统计单位和两次重要修复

将原始 RNA counts 按 **cell_population×biological sample** 汇总；不是把 nuclei 当作独立 replicate。只有四个样本各自至少 **25 nuclei** 的 population 进入正式 DE，共 13 populations、52 pseudobulk libraries。Endodermis 在 control_rep1 恰为 25，满足 ≥25 的边界。排除 Atrichoblast_like（min=4）、Phloem（8）、Phloem_related（8）、Xylem（8）；标签与描述性证据仍保留。

`AggregateExpression()` 将 group values 内的 `_` 转成 `-`，例如 `Unresolved_transition`→`Unresolved-transition`、`control_rep1`→`control-rep1`。原 eligible labels 与 pseudobulk metadata 字符串不一致，使 population matching 失败。原过程把 metadata population 的短横线恢复为下划线；整理脚本直接用原 population/sample metadata 构建转换后的 key 并核对完整 column alignment，避免按分隔符错误拆解复杂名称。

最初 Cortex trial 使用 merge 的 union feature universe：126,394 features，row-total≥10 后仍为 106,442。这个 trial 不作为正式结果。原四个样本 feature universe 不相同，merge 后缺失 feature 的结构性 0 会误导 DE；正式版本限制在原四个 feature lists 的 intersection，**32,833 features**，然后在各 population 内用四样本 raw-count row total≥10 预过滤。

正式模型为 **DESeq2 design=~condition**，contrast=stress/control，每个条件两个 biological samples。rep1/rep2 没有被假定为跨条件的配对 block。阈值为 **padj<0.05 且 |log2FoldChange|≥1**；up/down 指 P3 stress 相对于 control。保留 NA padj 的 all-result rows，不将它们计为 significant。

`genes_tested` 是 prefilter 后进入 DESeq2 result table 的行数，包含 NA padj 的行，并非全部通过 independent filtering 的 gene 数。正式结果：

| cell_population | genes_tested | DEGs | up | down |
| --- | --- | --- | --- | --- |
| Trichoblast | 19301 | 739 | 251 | 488 |
| Atrichoblast | 17862 | 687 | 488 | 199 |
| Pericycle | 17523 | 663 | 136 | 527 |
| Epidermal_like | 17650 | 535 | 275 | 260 |
| LRC_like | 19704 | 524 | 155 | 369 |
| Endodermis_like | 17186 | 504 | 82 | 422 |
| Initial_like | 19369 | 443 | 284 | 159 |
| Dividing_meristematic | 18184 | 410 | 323 | 87 |
| Cortex | 17150 | 259 | 86 | 173 |
| Dividing_S_phase | 17545 | 214 | 115 | 99 |
| Columella | 15649 | 151 | 55 | 96 |
| Endodermis | 13562 | 104 | 16 | 88 |
| Unresolved_transition | 13819 | 96 | 50 | 46 |

## 5. P1 candidate contextualization

正式输入是 **137 个 BRL3×drought interaction genes**，使用已存在的 P1 interaction DEG 文件，而非另选一个候选集合。共有 **129 个 distinct P1 genes** 出现在至少一个 population 的正式 P3 result table；总计 1,366 gene×population result rows，其中 155 rows 满足 P3 stress DE 阈值。129 是跨 populations 的去重数，不能与每个 population 分母混用。

每个 population 的 P1_tested 是 candidate 在该 population prefilter 后 result table 中的行数，包括 NA padj；因此比例分母因 population 而异。Pericycle 的比例最高，**23/118=19.49%（6 up、17 down）**。完整汇总如下：

| cell_population | P1_tested | P1_stress_DEG | up | down | percent_stress_DEG |
| --- | --- | --- | --- | --- | --- |
| Pericycle | 118 | 23 | 6 | 17 | 19.49 |
| Atrichoblast | 108 | 20 | 19 | 1 | 18.52 |
| Trichoblast | 119 | 19 | 9 | 10 | 15.97 |
| Epidermal_like | 107 | 16 | 11 | 5 | 14.95 |
| Dividing_meristematic | 100 | 14 | 14 | 0 | 14.00 |
| Initial_like | 112 | 15 | 14 | 1 | 13.39 |
| LRC_like | 121 | 13 | 6 | 7 | 10.74 |
| Endodermis_like | 111 | 10 | 2 | 8 | 9.01 |
| Unresolved_transition | 85 | 5 | 1 | 4 | 5.88 |
| Cortex | 106 | 6 | 4 | 2 | 5.66 |
| Columella | 94 | 5 | 2 | 3 | 5.32 |
| Dividing_S_phase | 98 | 5 | 5 | 0 | 5.10 |
| Endodermis | 87 | 4 | 0 | 4 | 4.60 |

这说明 P1 candidates 在 P3 WT osmotic stress 中具有 cell-population-specific context；不构成对 BRL3×drought interaction 效应或方向的重复验证。未在某 population result table 中出现的候选不能直接解释为不表达；它可能不在共同 feature universe 或未通过该 population 的 count prefilter。

## 6. 已有 P1 GO annotations 的 contextualization

直接读取已有 7 个 P1 enriched GO terms，不重新做 GO enrichment。geneID 中的 slash-separated symbols 拆成 **36 term-symbol memberships，13 distinct symbols**。最初用 P1 SYMBOL 字符串 reverse mapping 留下 23/36 NA；随后通过 `AnnotationDbi::select(org.At.tair.db, keytype="SYMBOL")` 全部修复。36/36 是 membership rows，不是 36 个不同 genes。最终映射固定保存在 reference 目录。

把每个 GO term 的 mapped genes 与 **P1 candidates 中的 P3 significant rows** 连接；这是预期的 many-to-many relationship（gene 可属于多个 terms、可在多个 populations 显著），各 term×population 内用 distinct TAIR 计数。20 个非零 term×population combinations，全部 up=0、down=stress_DE_genes。Pericycle 覆盖所有 7 annotations：cold acclimation、JA metabolism、long-chain fatty acid metabolism、acid response、cold response、water response、water deprivation；其中 **cold response 3 genes，全部 down**。

原 pivot 仅出现 8 个非零 populations，得到 7×9 表（含 Description）；补齐全部 13 eligible populations 后为 **7×14 表，数值矩阵为7×13**。零表示“没有满足当前 candidate overlap 和 DE 阈值的 gene”，不能解释为该 biological program 不存在。

JA metabolism / long-chain fatty acid metabolism 共享同一组四个 symbols；water / water-deprivation / acid response 共享同一组六个 symbols。七个 enriched annotations 有明显 gene overlap，不能称作七个独立 programs 被重复验证。完整非零汇总及补齐矩阵已归档。

## 归档与复现边界

本次完成的是组织代码、归档已有结果和一致性检查，没有重跑分析。历史决策、诊断和正式结果均在脚本中区分。完整结果轻量 checkpoint 与 cell metadata 已保存；完整 Seurat checkpoint 保留在原项目。当前导出环境记录可用，但最初九月运行的完整 sessionInfo、所有 QC 图、relaxed cluster12 原结果以及 signature 文献来源不完整。随后用户从原始矩阵完成全流程重跑；QC retained nuclei、PC ARI 矩阵、13 populations DE 汇总、P1 汇总及 GO 矩阵均与原归档一致。用户已授权提交 Core 1 到 GitHub。Core 2 尚未启动。
