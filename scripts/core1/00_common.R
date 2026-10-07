# Corrected Core 1. Source defines helpers only; run stages explicitly.
# All inference is rebuilt from gene-only counts. Frozen results are comparison only.
root <- Sys.getenv('P3_ROOT', getwd())
out <- file.path(root, 'results/core1')
# Use standard R libraries unless an optional local override is supplied.
lib <- Sys.getenv('P3_R_LIB', '')
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
if (!file.exists(file.path(root, 'data/rna/sample_metadata.tsv')))
  stop('Start R in the Project 3 root, or set P3_ROOT to that directory.')
for (d in c('tables','figures','checkpoints','diagnostics')) dir.create(file.path(out,d), recursive=TRUE, showWarnings=FALSE)
write_table <- function(x,name) utils::write.csv(x,file.path(out,'tables',paste0(name,'.csv')),row.names=FALSE)
checkpoint <- function(name) file.path(out,'checkpoints',paste0(name,'.rds'))
plot_save <- function(p,name,w=10,h=7) ggplot2::ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=160,bg="white")
samples <- utils::read.delim(file.path(root,'data/rna/sample_metadata.tsv'))
stopifnot(!anyDuplicated(samples$sample),all(table(samples$condition)>=2))
session_save <- function(stage) writeLines(capture.output(utils::sessionInfo()),file.path(out,'diagnostics',paste0(stage,'_sessionInfo.txt')))
set.seed(42)
options(future.globals.maxSize=8*1024^3)
