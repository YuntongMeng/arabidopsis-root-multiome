# Control 1 独立核查

核查日期：2026-10-12。未筛选细胞核，未覆盖原矩阵、fragments 或 QC 对象。

## 结论

**Control 1 的低 ATAC 信号群在原论文发布的矩阵中已经存在；本次独立复算没有发现当前核片段计数、细胞器比例或共同 peaks 矩阵的计算错误。** 这支持将差异视为输入数据中部分 barcode 的 ATAC 质量差异，而不是 Control 1 沿用旧 QC 结果所造成的假象。

这不表示整个 Control 1 都不可用，也没有定位到具体的实验原因。

## 独立复算做了什么

新脚本没有调用原批量脚本的计数函数，也没有用 Signac 重计数。它用 Python 直接读取原始 fragments 的全部 **49,108,250 行**，按 barcode 和染色体分别计数，并用二分查找计算片段与共同 peaks 的重叠。

- fragments 每行代表一个去重片段，第五列 PCR 支持数不作为权重。
- 核染色体为 1–5，线粒体为 Mt，叶绿体为 Pt。
- BED 坐标采用左闭右开；仅接触边界不算重叠。
- FRiP 分子中，每个片段最多计一次；矩阵中，跨越多个 peak 的片段分别贡献给这些 peak。
- 核小体信号复现当前 Signac 的实现：长度 <147 bp 为短片段，147<长度<294 bp 为单核小体片段；恰好 147 或 294 bp 不进入这两个分类。参见 [Signac 源码](https://github.com/stuart-lab/signac/blob/master/src/group.cpp)。

R 只用于读取已有 RDS/原始稀疏矩阵，以及汇总、画图。读取原论文 mixed-feature matrix 时，明确仅保留第三列为 `Peaks` 的行；RNA 行不参与比较。

## 数值核对结果

| 检查 | 结果 |
|---|---|
| 3,431 个 barcode 的总片段数、核片段数、命中共同 peaks 的独立片段数 | 全部完全一致 |
| 每个 barcode 的矩阵总计数与检测到的 peak 数 | 全部完全一致 |
| 每个 barcode 的 FRiP、细胞器比例、核小体信号 | 全部一致，仅存在浮点存储误差，最大约 5.7×10⁻¹⁴ |
| 27,677 × 3,431 矩阵，共 94,959,787 个位置 | 全部一致；核对了 3,884,627 个非零条目及其余零位置 |
| 核染色体 fragments 副本 | 解压内容与原 fragments 中仅保留 1–5 的内容完全一致 |
| 原 fragments 的校验值、各染色体行数 | 与上游保存的验证记录一致 |
| 核查前后原 fragments、核 fragments、原 Control 1 QC 表校验值 | 未改变 |

详细数值见同目录 `existing_vs_independent.csv` 与 `audit_checks.json`。

## 低信号是否在我们处理之前就存在

以相同 barcode 对照原论文发布的 ATAC 矩阵：

| 样本 | 原发布 ATAC 计数中位数 | 原发布计数与本次共同 peaks 计数的 Spearman 相关 |
|---|---:|---:|
| Control 1 | 797 | 0.9920 |
| Control 2 | 3,437 | 0.9984 |
| Stress 1 | 3,820 | 0.9987 |
| Stress 2 | 3,979 | 0.9989 |

不同流程的计数单位和 peak 集合不同，不能要求原计数等于新计数。这里检查的是：是否同一批 barcode 始终偏低或偏高。

Control 1 按核片段数做描述性分组，**500 仅用于诊断，不是本次确定的筛选阈值**：

| Control 1 分组 | barcode 数 | 原发布 ATAC 计数中位数 | 核片段数中位数 | 细胞器比例中位数 |
|---|---:|---:|---:|---:|
| 核片段数 <500 | 1,822 | 233 | 70 | 90.69% |
| 核片段数 ≥500 | 1,609 | 5,722 | 3,570 | 6.98% |

因此，Control 1 整体中位数被占比超过一半的低核信号 barcode 拉低。这些 barcode 在原发布矩阵中也偏低；有较多核片段的另一群，其细胞器比例明显较低。细胞器比例来自保留 Mt/Pt 的原 fragments，核 fragments 副本剔除 Mt/Pt 后仍可追溯这个筛选前指标。

## 核查范围与下一步

本次已独立复算 8 个指标并完整核对共同 peaks 矩阵。**没有重新计算 TSS enrichment，也没有从 FASTQ 重新比对。** 因此不能将本次结果解释为对所有上游假设和全部 QC 指标的绝对证明；也不能仅凭这些结果区分样本制备、捕获效率、测序等具体原因。

当前可以继续比较四个样本的联合 QC 分布，确定筛选规则，并记录每个样本的保留率；无需仅因这项差异重新跑 Control 1 的整套上游。Doublet 识别仍是后续独立步骤，本次没有执行。

脚本：`outputs/core2_control1_independent_audit.py`。运行方式：`python3 /absolute/path/to/outputs/core2_control1_independent_audit.py`。脚本会拒绝覆盖已存在的独立核查目录；再次运行前应先归档该目录。中间 R 读取脚本保存在 `work/control1_independent_audit/`，审计结果保存在本目录。
