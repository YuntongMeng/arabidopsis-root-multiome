# Core 1 integrity check

2026-10-05：仅检查静态脚本及已经完成的结果；没有重新运行生物学分析。

- 五份R文件（四主脚本＋helper）全部parse通过。
- 所有显式package::function引用均可在当前安装环境中解析。
- 在独立R会话中source全部脚本，仅定义函数，没有触发生物学分析。
- 未发现裸select/Select/filter/rename/mapIds/pivot_longer/pivot_wider调用。
- 原导出metadata：24,277 cells、18 clusters，annotation表18行。
- 已有DE表224,504行、13 populations；gene×population无重复。
- 依据保存的padj/log2FC列复核各population的DEGs/up/down，全部与de_summary一致。
- P1输入137 unique TAIR；129出现于至少一个population DE表。
- 逐population复核P1_tested、P1_stress_DEG、up/down及比例，与归档summary一致。
- 由原population×sample count表复核≥25 nuclei eligibility，恰为归档的13 populations。
- GO mapping36 memberships、13 symbols、全部有TAIR；20个非零combinations均全部down。
- 完整GO表7行14列（Description+13 populations）；Pericycle包含所有7 terms。
- Common feature list32,833项。
- 原strict marker表26,059行；cluster12没有returned marker。
- QC before只读取原barcode文件行数，after来自保存的metadata，没有重新过滤counts。

这些检查确认整理所引用的数字与已有表一致，不等同于验证完整分析重跑后结果完全一致。
