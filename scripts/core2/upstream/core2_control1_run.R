# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Control 1：只完成 QC 之前的 ATAC upstream，不进行细胞核筛选。
# 此入口调用保存在同一目录的 Python 流程和 R 计数脚本。
# 输入：已下载的 SRR24921616.sra、配对 barcode 参考表、统一 peak 坐标。
# 处理：完整 FASTQ 导出 → barcode 转换 → TAIR10 比对 → 去重/Tn5 shift
#       → fragments 排序、压缩及索引 → 统一 peak 重新计数。
# 完成标记允许复用已完成的导出、比对和 fragments，避免重复长步骤。

task_dir <- Sys.getenv("P3_CORE2_WORK", "")
stopifnot(nzchar(task_dir))
python <- file.path(task_dir, "work", "atac_env", "bin", "python")
pipeline <- file.path(Sys.getenv("P3_ROOT", getwd()), "scripts", "core2", "upstream", "core2_control1_upstream.py")
stopifnot(file.exists(python), file.exists(pipeline))

status <- system2(python, args = c(shQuote(pipeline), "all"))
stopifnot(status == 0L)

# 读取检查清单；这里不创建 Seurat 对象，也不执行 QC。
manifest <- utils::read.csv(file.path(
  task_dir, "outputs", "control1_upstream", "count_matrix_manifest.csv"
))
print(manifest, row.names = FALSE)

# 后续学习 QC 时再读取矩阵，不在这里筛选：
# atac_counts <- readRDS(file.path(
#   task_dir, "outputs", "control1_upstream", "control1_common_peak_counts.rds"
# ))
