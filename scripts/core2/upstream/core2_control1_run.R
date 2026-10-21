# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
# Control 1: ATAC upstream up to (not including) QC; no nucleus filtering.
# This entry point calls the Python pipeline and R counting script in the same directory.
# Inputs: downloaded SRR24921616.sra, paired barcode reference table, common peak coordinates.
# Steps: full FASTQ export -> barcode translation -> TAIR10 alignment -> dedup/Tn5 shift
#        -> sort, compress and index fragments -> recount on common peaks.
# Completion markers let finished export, alignment and fragment steps be reused, avoiding long reruns.

task_dir <- Sys.getenv("P3_CORE2_WORK", "")
stopifnot(nzchar(task_dir))
python <- file.path(task_dir, "work", "atac_env", "bin", "python")
pipeline <- file.path(Sys.getenv("P3_ROOT", getwd()), "scripts", "core2", "upstream", "core2_control1_upstream.py")
stopifnot(file.exists(python), file.exists(pipeline))

status <- system2(python, args = c(shQuote(pipeline), "all"))
stopifnot(status == 0L)

# Read the check manifest; no Seurat object is created and no QC is run here.
manifest <- utils::read.csv(file.path(
  task_dir, "outputs", "control1_upstream", "count_matrix_manifest.csv"
))
print(manifest, row.names = FALSE)

# Read the matrix later, during QC; no filtering here:
# atac_counts <- readRDS(file.path(
#   task_dir, "outputs", "control1_upstream", "control1_common_peak_counts.rds"
# ))
