# Decision 2 (after inspecting 01 quantiles/density/scatter): distributions are
# continuous; no uniquely justified knee. Controls have lower complexity and
# control_rep1 has higher mt baseline. Use sample-wise robust log1p outliers,
# not old mixed-feature thresholds. A 500-gene/count floor removes the sparse
# extreme tail, without imposing stress-sample complexity on controls.
# High count cutoff is deliberately lenient (4 MAD); scDblFinder handles doublets.
# Chloroplast percentages are low (sample medians <0.24%); diagnostic only.
run_qc_doublets <- function() {
  counts <- readRDS(checkpoint('01_gene_counts')); qc <- readRDS(checkpoint('01_qc'))
  thresholds <- list(); calls <- list(); retained <- list()
  bound <- function(v,n,side) exp(stats::median(log1p(v))+side*n*stats::mad(log1p(v)))-1
  for (i in seq_len(nrow(samples))) {
    s <- samples$sample[i]; md <- qc[qc$sample==s,]
    th <- data.frame(sample=s,min_features=max(500,ceiling(bound(md$nFeature_RNA,3,-1))),
       min_counts=max(500,ceiling(bound(md$nCount_RNA,3,-1))),
       max_counts=floor(bound(md$nCount_RNA,4,1)),max_mt=bound(md$percent_mt,3,1))
    md$fail_low_features <- md$nFeature_RNA<th$min_features
    md$fail_low_counts <- md$nCount_RNA<th$min_counts
    md$fail_high_counts <- md$nCount_RNA>th$max_counts
    md$fail_high_mt <- md$percent_mt>th$max_mt
    md$qc_pass <- !Reduce(`|`,md[,grep('^fail_',names(md)),drop=FALSE])
    # Alternative QC definitions are diagnostic counts, not alternate formal results.
    th$retained_mt4MAD <- sum(md$nFeature_RNA>=th$min_features & md$nCount_RNA>=th$min_counts & md$nCount_RNA<=th$max_counts & md$percent_mt<=bound(md$percent_mt,4,1))
    th$retained_gene_floor800 <- sum(md$qc_pass & md$nFeature_RNA>=800)
    th$retained_primary <- sum(md$qc_pass)
    thresholds[[s]] <- th
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts[[s]][,md$cell[md$qc_pass],drop=FALSE]))
    # Each GEO biological sample is one separately processed capture here.
    # Random artificial doublets suit continuous developmental trajectories.
    # Explicit standard 10x rate prior; no ground-truth loading-rate estimate.
    message('scDblFinder ',s,' n=',ncol(sce)); set.seed(4200+i)
    sce <- scDblFinder::scDblFinder(sce,clusters=FALSE,dbr.per1k=.008,
      BPPARAM=BiocParallel::SerialParam(RNGseed=4200+i))
    md$scDblFinder_score <- NA_real_; md$scDblFinder_class <- 'not_tested_QC_fail'
    j <- match(colnames(sce),md$cell)
    md$scDblFinder_score[j] <- sce$scDblFinder.score
    md$scDblFinder_class[j] <- as.character(sce$scDblFinder.class)
    md$retained <- md$qc_pass & md$scDblFinder_class=='singlet'
    calls[[s]] <- md; retained[[s]] <- counts[[s]][,md$cell[md$retained],drop=FALSE]
    write_table(md,paste0('qc_doublet_calls_',s)); rm(sce); gc()
  }
  md <- dplyr::bind_rows(calls)
  write_table(dplyr::bind_rows(thresholds),'qc_thresholds_and_sensitivity')
  write_table(md,'qc_doublet_calls_all')
  sm <- md |> dplyr::group_by(sample) |> dplyr::summarise(input=dplyr::n(),qc_removed=sum(!qc_pass),qc_pass=sum(qc_pass),doublets=sum(scDblFinder_class=='doublet'),singlets=sum(retained),old_qc_pass=sum(old_qc_pass),.groups='drop')
  write_table(sm,'qc_doublet_retention')
  saveRDS(retained,checkpoint('02_singlet_counts'),compress=FALSE)
  saveRDS(md[md$retained,],checkpoint('02_singlet_metadata'),compress=FALSE)
  plot_save(ggplot2::ggplot(md[md$qc_pass,],ggplot2::aes(scDblFinder_score,fill=scDblFinder_class))+ggplot2::geom_histogram(bins=50)+ggplot2::facet_wrap(~sample)+ggplot2::theme_bw(),'02_doublet_scores')
  session_save('02_qc_doublets')
}
