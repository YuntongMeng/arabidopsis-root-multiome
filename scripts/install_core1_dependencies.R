# Run once in the project root. Package versions of the reviewed run are recorded
# in results/core1/diagnostics; installing current packages is not version locking.
lib <- Sys.getenv('P3_R_LIB', '')
if (nzchar(lib)) {
  dir.create(lib,recursive=TRUE,showWarnings=FALSE)
  .libPaths(c(lib,.libPaths()))
}
options(repos=c(CRAN='https://cloud.r-project.org'))
cran <- c('Seurat','SeuratObject','Matrix','harmony','mclust','ggplot2','dplyr','tidyr','readxl','RANN','BiocManager')
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) utils::install.packages(missing)
bioc <- c('BiocParallel','scDblFinder','SingleCellExperiment','SummarizedExperiment','DESeq2','AnnotationDbi','org.At.tair.db')
missing <- bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
