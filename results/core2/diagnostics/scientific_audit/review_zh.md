# Core 2 科学性复核结论

日期：2026-10-13。核查对象是用户已经完成并保存的真实分析，以及生成这些结果的代码与记录。本轮不重新比对、选择新阈值、拟合模型或生成新差异结果。

## 判断

**在本次核实范围内，没有发现 RNA/ATAC 输入混用、坐标错置、跨样本 barcode 错配、错误计数汇总、以细胞数替代生物学重复、错误模型输入或差异结果导出等原则性科学错误。现有主线方法可以保留。**

这与“所有假设都已被验证”“细胞类型标签完全准确”或“所有生物学结论已确证”不同。以下明确区分实际核查与仍未验证的部分。

## 本轮实际核查的证据

| 核查项目 | 结果 |
|---|---|
| 四个原始 feature 表的 RNA/Peaks 类型 | 重新读取来源表，分别识别 Gene Expression 与 Peaks；RNA 与 ATAC 没有作为同类 features 进入本次 ATAC 计数 |
| 共同区域与坐标 | 从四份来源 Peaks 重新构建 union/reduce，逐个区间核对保存的 27,677 个区域；BED start+1 与 GRanges 一基闭区间转换一致 |
| 最终 ATAC 输入 | 27,677 ×14,146；全部是共同核区域的非负整数原始计数；与四份上游矩阵相应子矩阵逐值一致，未以 TF-IDF、GeneActivity 或 RNA counts 替代 |
| 实际保留核集合 | 独立从原 QC 表和 corrected Core 1 RNA 调用重建筛选逻辑，最终核 ID 集合完全一致；没有发现额外隐藏的样本或核删除 |
| 同核 RNA 标签 | sample＋barcode 唯一对应；14,146 个核的 condition 和 RNA population 与 corrected Core 1 主分支一致 |
| fragments 引用与 barcode 翻译 | 四个对象内的文件引用对应各自核 fragments；每个最终核恰有一个对应文件，文件中的 barcode 与该核原 barcode 相同，无跨 capture 混用 |
| QC 元数据沿袭 | 最终核的 nuclear fragments、TSS、核小体信号、nCount/nFeature、FRiP 和细胞器百分比与原 QC 表一致 |
| 参考注释 | TAIR10.62 注释坐标均在所用参考序列长度内；保存的 shared TSS 位置与重新导入 GTF 后按链方向生成的位置完全一致 |
| pseudobulk | 使用另一种直接逐组求和方式复算全部 46 个实际存在的 sample×RNA population 组合，所有 peak counts 及核数一致；没有补造缺失组的零样本 |
| 五个 DESeq2 模型 | 输入与原始 pseudobulk 过滤子矩阵逐值一致；每个模型是四个生物学样本、两个 control＋两个 stress；设计矩阵满秩 |
| 已导出的统计量 | stress/control 的系数和标准误、Wald 统计量及 p 值与保存模型一致；从保存 p 值复算的群体内 BH 值一致；显著候选拟合收敛 |
| 显著性与倍数条件 | FDR<0.05 和估计两倍变化是分开记录的；没有把倍数筛选误写成检验真实效应超过两倍 |
| P1 匹配 | 用独立的正负链区间算术复核基因体 OR TSS 上游2000/下游500的全部匹配，区间集合和两个重叠标记一致；137个参考基因均存在，120个有匹配 |
| 四个 SRA 与处理记录 | accession/sample 配对及记录的四路读长、同步检查和核心参数一致；本轮没有重读全量 FASTQ |

本轮 counts/model/P1 主核查共111项，通过。具体项目见 `checks.csv`；fragment 映射见 `fragment_mapping_checks.csv`；注释检查见 `annotation_checks.csv`。数量不是科学正确性的证明，重要的是上表实际检查的数据关系。

## 方法本身的依据

先建立共同 peaks，再从 fragments 对每个样本重计数，是 [Signac 官方合并流程](https://stuartlab.org/signac/articles/merging)所描述的路线。

使用生物学样本内的细胞群汇总计数，进行 replicate-aware pseudobulk 差异分析，有[单细胞 ATAC 方法评估研究](https://www.nature.com/articles/s41467-024-53089-5)支持。DESeq2 使用原始计数进行推断，VST 仅作诊断，与[官方说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)一致。这些依据支持方法选择，不能自动保证本项目每个具体发现都为真。

## 保留的局限与解释要求

- 每条件只有两个生物学重复，部分群体重复一致性弱。Root_cap 仍为探索性结果。未显著不等于没有胁迫响应。
- RNA 标签和 ATAC cluster 描述有不确定性，部分为混合谱系。RNA singlet 判定不是独立的 ATAC doublet 检验，也不能证明完全无 doublets。
- 主要 QC 阈值是合理工作选择，不是已证明唯一最优的标准；没有单独硬筛 FRiP 或细胞器比例，不等于这些指标未计算。
- DESeq2 比较的是经规模标准化后的相对 peak 开放性，不支持断言全局绝对开放性没有变化。
- BH 控制范围是每个 population 内，不能宣称整个五群体实验统一 FDR 控制。
- 附近 peak 不等于靶基因；GeneActivity 不等于 RNA 表达；P1 interaction contrast 不等于 P3 stress/control contrast。
- Control 1 低信号群的存在和已检查计数的正确性已有证据，具体实验原因仍未定位。

## 尚未做的验证

本轮没有从 FASTQ 重比对，也没有从 fragments 独立重新计算所有 TSS enrichment 数值；核对的是其元数据沿袭、共有 TSS 集和坐标。前次 Control 1 独立核查已覆盖完整 fragments 计数及完整矩阵，但不代表四样本全部步骤均由第二套实现全量重跑。

补充记录：这轮只读核查之后，用户已在 RStudio 完成一次新整理脚本的真实数据下游重跑，十项新旧对照全部一致。范围是从四个现成未筛选 QC 对象开始，不包含 FASTQ 重比对、fragment 重计数或 QC 指标重算；也不等于所有生物学解释获得独立实验验证。记录见 [手动重跑对照](../manual_rerun/validation_comparison.csv)。

**因此，当前没有发现需要推翻整个 Core 2 的原则性错误。合适的处理是保留有限结果及以上解释边界，而不是因为显著发现少而改变方法或阈值。**
