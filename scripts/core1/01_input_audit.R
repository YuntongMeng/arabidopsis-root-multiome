# Decision 1: Read the matrix and its row metadata together, BEFORE constructing
# a Seurat object or calculating library sizes. ReadMtx alone does not select RNA.
run_input_audit <- function() {
  counts <- list(); qcs <- list(); audits <- list(); feature_types <- list()
  for (i in seq_len(nrow(samples))) {
    s <- samples$sample[i]; message('Reading ',s)
    dir <- file.path(root,'data/rna',s)
    ff <- list.files(dir,pattern='_features.tsv.gz$',full.names=TRUE)
    mf <- list.files(dir,pattern='_matrix.mtx.gz$',full.names=TRUE)
    bf <- list.files(dir,pattern='_barcodes.tsv.gz$',full.names=TRUE)
    stopifnot(length(ff)==1,length(mf)==1,length(bf)==1)
    f <- utils::read.delim(gzfile(ff),header=FALSE,stringsAsFactors=FALSE)
    b <- utils::read.delim(gzfile(bf),header=FALSE,stringsAsFactors=FALSE)[[1]]
    x <- methods::as(Matrix::readMM(gzfile(mf)), 'CsparseMatrix')
    stopifnot(nrow(x)==nrow(f),ncol(x)==length(b),ncol(f)>=3,
              !anyNA(f$V3),!anyDuplicated(f$V1),!anyDuplicated(b),
              all(is.finite(x@x)),all(x@x>=0),all(x@x==round(x@x)))
    # Assignment precedes positional subsetting; preserves exact row correspondence.
    dimnames(x) <- list(f$V1,paste(s,b,sep='_'))
    keep <- f$V3=='Gene Expression'; g <- x[keep,,drop=FALSE]
    stopifnot(nrow(g)>0,identical(rownames(g),f$V1[keep]),
              !any(grepl('^[^:]+:[0-9]+-[0-9]+$',rownames(g))))
    nc <- Matrix::colSums(g); nf <- Matrix::colSums(g>0)
    md <- data.frame(cell=colnames(g),sample=s,condition=samples$condition[i],
      nCount_RNA=nc,nFeature_RNA=nf,
      percent_mt=100*Matrix::colSums(g[grepl('^ATMG',rownames(g)),,drop=FALSE])/pmax(nc,1),
      percent_cp=100*Matrix::colSums(g[grepl('^ATCG',rownames(g)),,drop=FALSE])/pmax(nc,1),
      mixed_nCount=Matrix::colSums(x),mixed_nFeature=Matrix::colSums(x>0))
    md$old_qc_pass <- md$mixed_nFeature>=800 & md$mixed_nFeature<=15000 & md$mixed_nCount<=50000
    counts[[s]] <- g; qcs[[s]] <- md
    audits[[s]] <- data.frame(sample=s,total_features=nrow(x),gene_features=nrow(g),
      excluded_features=sum(!keep),peak_features=sum(f$V3=='Peaks'),nuclei=ncol(g),
      coordinate_ids_after=sum(grepl('^[^:]+:[0-9]+-[0-9]+$',rownames(g))),
      gene_UMI=sum(g),excluded_counts=sum(x)-sum(g))
    feature_types[[s]] <- data.frame(sample=s,as.data.frame(table(f$V3)))
    rm(x); gc()
  }
  qc <- dplyr::bind_rows(qcs)
  write_table(dplyr::bind_rows(audits),'input_matrix_audit')
  write_table(dplyr::bind_rows(feature_types),'input_feature_types')
  write_table(qc,'qc_all_barcodes')
  qq <- dplyr::bind_rows(lapply(names(qcs),function(s) dplyr::bind_rows(lapply(c('nCount_RNA','nFeature_RNA','percent_mt','percent_cp'),function(v) {
    q <- stats::quantile(qcs[[s]][[v]],c(0,.005,.01,.025,.05,.25,.5,.75,.95,.975,.99,.995,1))
    data.frame(sample=s,metric=v,quantile=names(q),value=as.numeric(q))
  }))))
  write_table(qq,'gene_only_qc_quantiles')
  common <- Reduce(base::intersect,lapply(counts,rownames))
  writeLines(common,file.path(out,'diagnostics/common_gene_features.txt'))
  saveRDS(counts,checkpoint('01_gene_counts'),compress=FALSE)
  saveRDS(qc,checkpoint('01_qc'),compress=FALSE)
  long <- tidyr::pivot_longer(qc,cols=dplyr::all_of(c('nCount_RNA','nFeature_RNA')),names_to='metric',values_to='value')
  plot_save(ggplot2::ggplot(long,ggplot2::aes(value,colour=sample))+ggplot2::geom_density()+ggplot2::scale_x_log10()+ggplot2::facet_wrap(~metric,scales='free')+ggplot2::theme_bw(),'01_gene_only_qc_density')
  plot_save(ggplot2::ggplot(qc,ggplot2::aes(nCount_RNA,nFeature_RNA,colour=condition))+ggplot2::geom_point(size=.3,alpha=.2)+ggplot2::scale_x_log10()+ggplot2::scale_y_log10()+ggplot2::facet_wrap(~sample)+ggplot2::theme_bw(),'01_gene_only_qc_scatter')
  session_save('01_audit')
}
