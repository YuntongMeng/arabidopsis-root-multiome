# Core 1：从最初判断到当前分析

这是分析过程的阅读主线。当前流程已实际运行；这次目录整理没有重新计算生物学结果。
“当前”表示目前证据支持的选择，不表示所有问题都已有唯一答案。

## 1. 我们最初怎样理解数据

最初把 multiome feature matrix 作为 RNA 输入，按 feature ID 读入，却没有按第三列区分 Gene Expression 与 Peaks。
当时做出了 QC、聚类、注释及 pseudobulk，也讨论过坐标形式的 markers，但没有及时追查其输入来源。
这些判断和代码留在 archive/core1_initial，原结果留在 results/core1_frozen；旧结论不能自动作为正式结果继续引用。

## 2. 什么证据改变了判断

坐标形式的 marker 和异常 feature 数促使我们重新审计原始 features 文件。
直接核对证明每份样本有 32,833 个 gene features，另外含 19,858–25,001 个 ATAC peaks。
新读入先对齐矩阵与 feature 行，再只保留第三列严格等于 Gene Expression 的行；RNA 中坐标 peak IDs 为零。

这也修正了我们对错误后果的解释：旧版 2,000 个 HVGs 中坐标 peaks 为零。
混合输入确定影响 QC、归一化分母和 marker 测试，但不能声称已经观察到 peaks 直接主导旧 PCA。
新流程同时改变了多项处理，没有单因素消融实验，所以不能把所有前后差异都只归因于 ATAC。

## 3. 为什么重新选择 QC，并加入 doublet 检查

去除 peaks 后，旧 nCount/nFeature 阈值不再代表同一量。我们重新检查每个样本的分布、分位数和散点图。
分布呈连续尾部，缺少唯一明确的截断点；control 的 RNA complexity 较低，control_rep1 的 mt 基线较高。
因此采用样本内 MAD 阈值及最低 500 的保护下限，不硬套旧上限。具体规则和敏感性检查见 decision_log.md。

24,474 个输入 nuclei 经 QC 去除 430 个，再由按 sample/capture 的 scDblFinder 标记 2,494 个预测 doublets，保留 21,550 个。
这些是算法预测；加载率先验和同型 doublet 漏检仍有不确定性。

GEO RAW.tar 名称并不证明包含 raw droplets。核查本地文件后，没有确认存在空液滴矩阵，所以没有强行做 ambient correction。

## 4. 为什么同时保留 unintegrated 与 Harmony

原先“研究 stress 所以不能整合”的理由过于绝对，但“必须 Harmony”也不是已证实结论。
我们让同一批 gene-only singlets 进入两条 embedding 分支，以原始 RNA counts 做后续 DE。
Harmony 提高 replicate mixing，但 canonical identity consistency 没有相应改善，同时削弱 condition-associated geometry。
因此当前以 unintegrated 为主，Harmony 作为敏感性比较；这不能证明所有变化都是过度校正。

PC 网格显示 20–30 PC 局部较稳定，我们选择 25 PC 做后续比较。
40–50 PC 还有另一稳定区域，与 25 PC 的划分有明显差异，所以没有宣称全范围稳健。
新的 18/17 个 clusters 是运行输出，不是为了复现旧数字而设定的目标。

## 5. 为什么重新注释，而不继承旧 cluster 名称

新 partition 重新使用 Shahan 2022 Data S2 signatures、canonical markers 和 dataset-specific markers 检查。
原始文献引用字符串和 spatial domains 保留在 reference 表中；同来源证据不算完全独立验证。
我们优先保留较宽的 lineage；Root_cap 不强拆 LRC/columella，meristematic 和部分 initial-like groups 保留低置信度。
Dividing_meristematic 是细胞状态，不能简单等同一个 lineage。

旧 LRC_like 的许多细胞在新结果中分到 Atrichoblast/Cortex，所以旧 LRC-specific 解读需要撤回。
Root_cap 与旧 Columella 有较强细胞连续性，但这不足以恢复全部精细标签。
人工 mapping 与 reviewed partition 的 MD5 绑定，未来 clustering 变化必须重新检查注释证据。

## 6. 怎样让 DE 对应生物学重复

用 population × biological sample 汇总 raw gene counts，要求每个样本至少 50 个细胞。
正式基因背景是四样本共有的 32,833 个 genes；DESeq2 使用 ~condition，不臆造配对。
两条分支分别计算 DE，报告 membership 变化导致的敏感性。

当前 eligible populations 在 25/50/100 细胞阈值下相同，这只说明本次 eligibility 稳定。
两组各两个 biological replicates 仍限制推断能力，较低置信度注释不能支持精细 cell-type 结论。

## 7. P1 contextualization 为什么增加背景检验

137 个 interaction candidates 与 stress DEGs 的 overlap 本身不能证明 cell-specific interaction。
Fisher + BH 显示，相对于普通 stress-DEG 背景，候选集在各 eligible populations 中富集。
但改用 P1 非interaction drought DEGs 作对照后，两条分支均没有 population 达到额外富集的 BH 显著性。
因此保留 P1/P3 stress context 的联系，不能据此宣布 Pericycle 特异或 BRL3 interaction 已获验证。

P1 TSV 中引号曾导致解析警告。修正 quote/comment 选项后重新运行 contextualization，旧警告运行被替代。
这个实现修正也保留在 decision log，而不是删去过程记录。

## 8. 七个 GO terms 能说明什么

七行 GO 只有 13 个独特候选 genes，gene sets 互相重叠。
当前两条分支仅三个 terms 有非零 DEG overlap；这是已有候选集的 context，不是 P3 全基因 GO enrichment。
零 overlap 不能证明对应生物程序不存在，共享 genes 也不能当成多条独立证据。

## 9. 怎样继续阅读

先看本文，再看 decision_log.md 的参数依据，最后用 comparison_report.md 和当前 figures/tables 核对数值。
version_comparison.md 给出旧结论能保留到什么程度；manual_review.md 提供审阅顺序。
新证据可更新当前选择，但应记录理由和版本差异，而不是抹去旧判断。
