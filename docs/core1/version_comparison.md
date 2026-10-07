# Core 1 版本比较与结论状态

| 项目 | 初版 | 当前修正版 | 解释 |
|---|---|---|---|
| 输入 | genes + peaks 混入 RNA | Gene Expression only | 初版输入错误确定存在 |
| QC | mixed-count 阈值 | gene-only 分布驱动阈值 | 两者指标不能直接等同 |
| doublets | 未加入本次方法 | scDblFinder per capture | 预测结果，非真值 |
| 保留细胞 | 24,277 | 21,550 | 21,488 shared cells；多项处理同时改变 |
| clustering | 旧 partition | unintegrated 主分析 + Harmony sensitivity | 共享细胞 old/new PCA ARI 0.639 |
| annotation | 部分精细标签证据不足 | sourced broad lineages + uncertainty | 不继承旧 cluster ID 标签 |
| pseudobulk | 旧分组 | 新分组 × biological sample | DE 数量不作为复现目标 |
| P1 overlap | Pericycle 23/118 | PCA 19/120；Harmony 19/117 fitted | 分母为 fitted candidates，不等于 Fisher 背景分母 |
| P1 enrichment | overlap 为主 | Fisher + BH，另有 drought-DEG control | 无额外 control-set 富集证据 |
| GO | 已有七个候选 gene sets | 完整零值、overlap 和解释边界 | 非新 P3 GO enrichment |

可以保留：raw-count pseudobulk 的分析思路、P1/P3 stress context 的探索意义，以及部分较宽细胞身份的连续性。
必须更新：旧 QC 指标、cluster membership、annotation、按旧分组计算的 DE 与 context 数值。
不能保留为定论：旧 LRC-specific 解释、Pericycle uniquely strongest 或 cell-specific BRL3 interaction validation。

旧 HVGs 中坐标 peaks 为零；不能把全部 embedding 变化归因于 peaks 直接进入 PCA。
完整数值、各群分支 concordance 和文献来源见 [详细比较](comparison_report.md)。
