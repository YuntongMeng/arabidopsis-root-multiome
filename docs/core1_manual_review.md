# Core 1 人工检查清单

- [ ] 逐份阅读四份主脚本，确认顺序、注释和中途方法调整符合实际过程。
- [ ] 核对四个 sample 的 sample/condition/replicate、RNA 路径与 feature ID 列。
- [ ] 核对 QC 阈值与 24,277 retained nuclei；mt/cp 作为诊断指标。
- [ ] 查看完整 PC ARI 矩阵及 PC25 的决定；确认 resolution=0.5 和18 clusters。
- [ ] 查看 clusters 6/14 的 sample enrichment 与 QC，不把 replicate enrichment 自动解释成独立生物学身份。
- [ ] 确认 JoinLayers、统一 strict marker 表和 cluster12 relaxed diagnostic 的区别。
- [ ] 核对全部 signature genes、gene ID 和文献来源，重点检查原 v3 candidate/placeholder 项；保留证据不足的 conservative labels。
- [ ] 核对18 clusters→17 populations 的最终 annotation/confidence，以及0/5合并为LRC_like。
- [ ] 核对13 eligible populations与每样本≥25 nuclei；Endodermis正好满足边界。
- [ ] 核对 pseudobulk sample keys、四样本共有32,833 features和row-total≥10。
- [ ] 核对 DESeq2 ~condition、stress/control方向、padj<0.05且|log2FC|≥1，以及全部DE汇总。
- [ ] 核对 P1 137-gene输入、129 genes出现在至少一个P3结果表、各population分母，以及Pericycle23/118。
- [ ] 核对36 GO memberships→13 symbols的mapping；检查完整7×13数值矩阵。
- [ ] 保持P1/P3 contextualization解释边界及GO gene overlap限制。
- [ ] 如将来实际重跑，先核对环境，并审查cluster partition/IDs；脚本不会对改变后的cluster编号直接套用旧annotation。
- [ ] 完整重跑及关键汇总核对已通过；用户已授权GitHub提交。signature来源仍需进一步科学审阅。

已有证据不足的部分：最初完整sessionInfo、全部原QC图与quantile输出、cluster12 relaxed原表、signature文献provenance。这些缺口被记录，未用新计算或编造信息补齐。
