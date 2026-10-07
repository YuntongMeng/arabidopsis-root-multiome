# Core 1 修正版：实际重跑与旧版比较

2026-10-06 · 本地 Project 3 · 已运行完成；用户人工运行 01–07 的 summary 与八项结构检查已核对。最新发布整理见 decision_log.md 与 manual_rerun_validation.md。

## 结论先读

**这次结果不是旧版的原样复现。无需你再从零重复跑一次才能知道差别，但旧版的生物学结论必须据新结果修改。** 本轮从原始 multiome 矩阵重新筛选 Gene Expression、QC、scDblFinder、PCA/Harmony、marker/annotation、DESeq2、P1 Fisher 和 GO，未把 frozen 结果当作新结果复用。

- **仍有支持：** population × biological sample 的 pseudobulk 设计；P1 候选具有 P3 osmotic-stress context；Pericycle 的候选相对同群体可检验基因背景显著富集；若干群体的 DE 方向和相当部分基因与旧版一致。
- **需要撤回或重写：** 旧 cluster 编号/标签；“Pericycle overlap 比例最高”；“Pericycle 涉及全部 7 个 GO”；“GO-linked genes 全部下调”。
- **新增解释边界：** P1 interaction 候选相对普通基因背景显著富集，但相对 **P1 非 interaction drought DEGs**，两个分支都未检出额外富集。不能据此宣称 P3 提供了 BRL3 interaction 特异性验证。
- **一个纠正先前 review 的实证：** 旧 checkpoint 的 2,000 个 HVGs 中，坐标型 peak 为 **0**。ATAC 确实污染输入、QC、归一化分母和 marker 检验，但“peaks 大量进入 HVG、直接主导 PCA”没有获得这个 checkpoint 的支持。

本次同时改变了输入处理、QC、doublets 和 annotation，另加 Harmony 分支；各结果差异不能全部单独归因于 ATAC 过滤。没有做逐项因素消融实验。

## 1. 输入修正与 QC

原每个样本混有 19,858 / 24,317 / 24,490 / 25,001 个 Peaks。严格按 feature 文件 **第 3 列等于 Gene Expression** 保留行，逐一校验矩阵行顺序、ID、维度、非负整数计数及唯一 barcode。每个样本均为 32,833 genes，共同 gene universe 也为 32,833；筛选后坐标型 peak ID = 0。不是因为数字接近基因总数就认定它们是基因。

先查看 gene-only quantiles、density 和 scatter，再制定样本内 log1p-MAD 阈值：低 features/counts 用 3 MAD、最低 500；极端高 counts 用 4 MAD；高 mt 用 3 MAD。实际最低 features 为 500/500/583/602，最高 mt 为 10.41/4.76/4.53/4.57%。cp 仅诊断。没有沿用旧 mixed-feature 的 800–15,000 features / 50,000 counts。

| 样本 | 输入 nuclei | QC 去除 | 预测 doublets 去除 | 保留 singlets |
| --- | --- | --- | --- | --- |
| control_rep1 | 3431 | 85 | 236 | 3110 |
| control_rep2 | 6631 | 180 | 666 | 5785 |
| stress_rep1 | 7827 | 94 | 917 | 6816 |
| stress_rep2 | 6585 | 71 | 675 | 5839 |

从 24,474 输入 nuclei 中，QC 去除 430、scDblFinder 再去除 2,494，最终 **21,550**。每个 capture 独立检测；rate prior=0.008/1,000 captured cells，版本及随机种子已保存。这是预测 doublets，不是已确认的 doublets；实际 loading rate 未知，homotypic doublets 可能残留。

旧版保留 24,277，与新版共享 21,488；旧独有 2,789、新独有 62。用原 mixed-count 规则重建的旧 barcode 集合与冻结 metadata 完全匹配。阈值的 800-gene floor 和 mt 4-MAD 备选只做了保留数量诊断，没有声称已完成整套 QC 参数敏感性重跑。

**Ambient RNA：** 本地 GEO RAW.tar 只有同样的 12 个 processed matrix/features/barcodes 文件，没有已验证的 raw/empty droplets。因此没有强行运行 SoupX，也没有把跨谱系表达直接视为污染比例。canonical DotPlot 可供检查，但不能排除环境 RNA。[GEO 数据来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495)

## 2. PC robustness 与 Harmony 的取舍

共同使用 gene-only singlets 和 RNA assay，LogNormalize 10,000、初选 2,000 HVGs；移除其中 1 个 organelle HVG 后，embedding 实际用 1,999 个。与旧版 HVG 共有 1,524 个基因。原始 counts 保留全部 gene features，逐样本确认到最终对象前后计数总和完全一致。

两分支均计算 10/15/20/25/30/40/50 dimensions、resolution=0.5。选 25 PCs 是依据 20–30 区域的局部稳定性与维度简约性：未整合 ARI(20,25)=0.809、ARI(25,30)=0.880；Harmony 为 0.886、0.903。40–50 PCs 形成另一稳定区域，和 25 PCs 的 partition 仍有明显差异，不能称整个 PC 范围完全稳定。

25 PCs 下未整合恰好仍为 **18 clusters**，Harmony 为 **17**；这并不意味着复现旧 partition。旧/新未整合在共享 nuclei 上的 **ARI=0.639**；新未整合/Harmony 为 0.622。

| 诊断 | 未整合 | Harmony |
| --- | --- | --- |
| 同 condition、同较可信 atlas identity 内的 replicate mixing ratio，分组平均 | 0.876 | 0.984 |
| broad atlas identity 邻域一致比例 | 0.485 | 0.485 |
| 同 condition 邻居比例 | 0.859 | 0.600 |
| 稀疏 stress-marker 原表达与邻居平滑值 Spearman | 0.229 | 0.184 |

选择 **未整合作主分析，Harmony 作敏感性分析**。Harmony 改善 control 重复混合，但没有相应提高这项 identity 一致性指标，并削弱 condition-associated geometry。这里 sample 嵌套于 condition，没有独立真值能证明被去除的部分全是 batch 或全是生物信号。四个 stress diagnostic genes 很稀疏，这个相关指标只是辅助证据。低维 condition 分离本身也不能证明生物信号保留正确。两条分支都实际跑了 raw-count DE，以下并列报告。

## 3. Annotation：有来源，也保留不确定性

主面板来自 [Shahan et al., 2022, Developmental Cell](https://doi.org/10.1016/j.devcel.2022.01.008)：原始 **Data S2K** 的 14 类、每类 top 50 signature genes；**S2F/G** 的 canonical markers 和其原始来源、空间表达范围；**S2J** 的 cell-cycle markers。原 xlsx、抽取后的面板、引用和每个新 cluster 的理由均在交付包内。不是把旧 placeholder 列表换个名字继续使用。

AddModuleScore 仅作支持，与新数据 strict markers 和 canonical DotPlot 联合判断。严格 marker 标准统一为 Wilcoxon、min.pct=.25、logFC=.25、返回 p<.01；marker 表另含调整后 p 值。marker 的细胞级 p 值不用于取代生物重复 DE。canonical panel 与 atlas panel 来自相关文献，不能声称完全独立验证。

- 用 broad lineage 合并有共同身份支持、但 developmental/condition 状态不同的子群。例如 Pericycle 的 control-enriched 与 stress-enriched 未整合子群合并后做 DE。
- Atrichoblast/Cortex/Endodermis/Trichoblast 的 `_lineage` 包含强弱不同的子群，不能把弱子群升级为确定的成熟细胞身份。
- **Root_cap** 有 SMB/BRN1/BRN2 等证据，但共享 marker 不足以可靠区分 LRC 与 columella。
- **Meristematic_unresolved、Epidermal_initial_like** 保留低置信度；**Dividing_meristematic** 是有 cell-cycle 证据的状态，不是单一谱系。
- Phloem、小型 unresolved、样本特异状态保留描述；不满足每个 replicate 的门槛就不进入正式 DE。

旧 LRC_like 的细胞大量分散到新 Atrichoblast_lineage 和 Cortex_lineage。因此旧 LRC 特异性解释不能直接继承。新 Root_cap 中 88.6% 细胞来自旧 Columella，支持根冠层面的连续性，但不据此强行恢复旧细标签。

## 4. Pseudobulk 与 DE

每个 population × biological sample 汇总**原始整数 gene counts**；共同基因 universe 内，四样本 row-total≥10 预过滤；DESeq2 `~condition`，stress/control，padj<.05 且 |log2FC|≥1。没有把 rep1/rep2 伪装成配对设计，也没有把细胞当 replicate。

主分析有 **9 populations**、Harmony 有 **10** 达标；入选群体最小样本数均≥100，所以本轮 ≥25 / ≥50 / ≥100 三种 eligibility 给出同一集合。两条件各 n=2，结论仍为探索性。低置信度 population 的 DE 只能按该保守标签解释。

| Population | 主分析 DEGs | Harmony DEGs |
| --- | --- | --- |
| Atrichoblast_lineage | 1.56e+03 | 1.55e+03 |
| Cortex_lineage | 1.46e+03 | 846 |
| Dividing_meristematic | 395 | 388 |
| Endodermis_lineage | 524 | 439 |
| Epidermal_initial_like | 403 | 440 |
| Meristematic_unresolved | 504 | 316 |
| Pericycle | 637 | 630 |
| Procambium_like | NA | 239 |
| Root_cap | 148 | 155 |
| Trichoblast_lineage | 1.75e+03 | 1.53e+03 |

表内 Procambium_like 的主分析为 NA，表示该标签只在 Harmony 分支定义，不是零个 DEGs。

这些是重新定义的 population，不能把某行数字与旧版同名或近名群体视为严格一一对应。全部 fitted/valid-pvalue/finite-padj 分母、原始 pseudobulk、样本 correlation、细胞数、size factors 和全量 DE 均已保存。

**新两分支的稳定性：**

| Population | log2FC Spearman | DEG Jaccard | 细胞成员 Jaccard |
| --- | --- | --- | --- |
| Atrichoblast_lineage | 0.946 | 0.717 | 0.827 |
| Cortex_lineage | 0.851 | 0.339 | 0.713 |
| Dividing_meristematic | 0.99 | 0.882 | 0.949 |
| Endodermis_lineage | 0.947 | 0.787 | 0.908 |
| Epidermal_initial_like | 0.968 | 0.809 | 0.819 |
| Meristematic_unresolved | 0.896 | 0.48 | 0.528 |
| Pericycle | 0.901 | 0.705 | 0.725 |
| Root_cap | 0.963 | 0.836 | 0.904 |
| Trichoblast_lineage | 0.988 | 0.825 | 0.858 |

所有共同显著基因的符号在两分支一致；但 Cortex 的 DEG Jaccard 仅 0.339，Meristematic_unresolved 为 0.480，成员和检出数量明显依赖分支。Pericycle 两分支有 524 个共同 DEG、log2FC Spearman=0.901；分裂状态和 Root_cap 更稳定。不能把“方向一致”解释为所有群体完全稳健。

与旧版成员接近的群体中，Dividing_meristematic 的 log2FC Spearman=0.945、DEG Jaccard=0.709；Pericycle 为 0.873 / 0.631；新 Root_cap 对旧 Columella 为 0.900 / 0.738。详细比较带成员交集，避免只比较名字。

## 5. P1 137 候选与 Fisher enrichment

Fisher 主背景为每个 population 内 **padj 非缺失、确实可进入显著性判定**的 genes。表格为 P1 / 非 P1 × stress DEG / 非 DEG，单侧富集检验，BH 在每个分支的各 population 内校正。另保存含 independent-filtered genes 的全部 fitted 背景敏感性版本，不混用分母。

| Population | P1 stress DEG | 可检验 P1 | OR | Fisher BH |
| --- | --- | --- | --- | --- |
| Atrichoblast_lineage | 26 | 108 | 3.28 | 3.83e-06 |
| Cortex_lineage | 19 | 87 | 2.38 | 0.00177 |
| Dividing_meristematic | 14 | 73 | 9.13 | 2.25e-08 |
| Endodermis_lineage | 11 | 69 | 3.7 | 0.000667 |
| Epidermal_initial_like | 11 | 70 | 6.59 | 5.58e-06 |
| Meristematic_unresolved | 17 | 87 | 7.74 | 1.2e-08 |
| Pericycle | 19 | 97 | 5.51 | 8.48e-08 |
| Root_cap | 5 | 56 | 6.34 | 0.00177 |
| Trichoblast_lineage | 30 | 102 | 3.47 | 4.8e-07 |

主分析所有 9 个 population 相对普通背景均有富集，Harmony 的 10 个也均有；这并不证明精确 cell-type specificity，尤其低置信度标签不可精确定位。每个 Fisher family 分别校正；没有把多种备选检验混成一次全局 BH。

**Pericycle：** 主分析 19 个 P1 stress DEGs（6 up、13 down），19/120 fitted=15.83%；真正 Fisher 背景中为 19/97，OR=5.51，BH=8.48×10⁻⁸。Harmony 也有 19 个（6 up、13 down），19/117 fitted，19/78 finite-padj，OR=5.93，BH=1.45×10⁻⁷。旧版是 23/118，不应继续引用。主分析最高 fitted overlap 比例现在是 Trichoblast_lineage，30/115=26.09%，不是 Pericycle；比例排名本身仍不是富集特异性检验。

**更合适的对照：** P1 WT drought DEGs 共 1,480；去掉 137 interaction 候选中的重叠者后，得到 1,442 个对照；再要求 interaction padj 有值且≥.05，得到 1,389 个严格对照。相对两种对照，两个分支**没有任何 population 达到 BH<.05 的 interaction 候选额外富集**。Pericycle 主分析相对第一种对照 OR=1.39、BH=0.437；Harmony OR=1.36、BH=0.688。

所以保留“候选与独立 osmotic-stress response 收敛”的描述，撤去“interaction-specific validation”或“Pericycle 唯一/最强”的表述。未显著不证明两集合等效；没有做 expression-matched control，表达量及检出功效混杂仍在。

## 6. 七个 GO annotations

重新映射后仍为 7 terms、36 term–gene memberships、13 distinct genes。使用的是 **P1 原 enrichment 表中列出的 candidate members**，不是完整 GO 注释基因库；此次没有新跑 P3 GO enrichment。

两个分支只有 **3/7 terms** 有至少一个满足当前 P1×P3 DEG 条件的 overlap：JA metabolism、long-chain fatty acid metabolism、response to cold。cold acclimation / water / water deprivation / acid response 均为 0。这是阈值下没有候选 overlap，不能解释成程序不存在。

Pericycle 主分析 cold response 从旧版 3 genes 降至 **1**，Harmony 为 **2**。JA 和长链脂肪酸条目各有 1 个下调基因，但它们共享成员，不能当两个独立发现。Atrichoblast_lineage 中这两个条目的 overlap 现在为上调，因此旧版“GO-linked 全部 down”不再成立。

全 7×population 矩阵包含零列/零行，pairwise gene-set overlap/Jaccard 也已导出。原七个 term 本来就强烈重叠，不能称为七个独立程序被重复验证。

## 7. 哪些旧结论可以继续用

| 旧结论或方法 | 本次处理 |
| --- | --- |
| P1 → population → stress pseudobulk → contextualization 的研究框架 | 保留 |
| biological sample 为统计重复；gene-only common universe | 保留并直接验证 |
| PC25、18 clusters | PC25 有新稳定性支持；18 只是未整合分支的巧合，partition 不同 |
| 旧 13 populations 和标签 | 重做；正式 DE 为 9/10，低置信度标签保守解释 |
| Pericycle 有较强 stress/P1 context | 两分支继续支持，但不能称唯一、最高或 interaction 特异 |
| Pericycle 23/118 | 替换为新分支的结果及清楚区分的分母 |
| Pericycle 全部 7 GO、全部 GO overlap 下调 | 撤回，改用上述 3-term 结果及方向 |
| 完全不需要考虑 Harmony | 改为实测敏感性比较；目前未整合作主分析 |
| peak 直接大量进入 HVG 驱动 PCA | 旧 checkpoint 不支持；不要继续当事实陈述 |

## 8. 文件、复现和审阅

本地完整 checkpoint 在 `Project 3/results/core1/checkpoints/`，包括 gene-only counts、singlets、两分支对象和 pseudobulk；约 5.6 GB，本审阅包不重复复制它们。包内提供 scripts、全部 CSV 表格、图、面板来源、方法记录、运行日志与 sessionInfo。原始 matrix 文件保留原位。

**入口：** 项目根目录 `Rscript scripts/run_core1_corrected.R help`。可逐阶段运行 audit / qc / embeddings / evidence / de / context / compare / diagnostics，或显式 all。R 函数冲突处均写 `package::function()`。默认使用标准 R libraries；可用未跟踪的本机 P3_R_LIB 设置覆盖。common helper 不再含个人绝对路径。

手工 annotation 是本次 evidence review 的产物，绑定 partition MD5。未来 partition 变化会要求重新审查，不能跳过检查后套标签。旧混合输入的 stage-1 入口已加阻止误运行提示，旧代码和 frozen archive 仍保留。

已检查：10 份新 R 文件语法；8 项结构性 invariants；152 个 Fisher p 值用独立 hypergeometric 求和复核、BH 用独立排序算法复核；27 个 frozen 文件与既有 SHA256 记录一致。图已检查可读性。**没有 GitHub commit/push。**

建议你先审阅每个 cluster 的 annotation 理由、canonical DotPlot，以及 Cortex/Meristematic 的分支敏感性，再决定是否接受这些 population 定义；不必为了核对本报告而重跑所有步骤。

主要方法与数据来源：[Shahan atlas 及补充表](https://pmc.ncbi.nlm.nih.gov/articles/PMC9014886/)、[scDblFinder 作者说明](https://plger.github.io/scDblFinder/articles/scDblFinder.html)、[P3 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495)、[root osmotic-stress marker 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10573044/)。
