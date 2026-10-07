# Explicit reviewed mapping is specific to saved corrected partitions. It is
# evidence-derived, not a carry-over of old cluster IDs. Freeze mapping only
# after 04 marker/score/DotPlot inspection; changed partitions require review.
run_annotation_pseudobulk <- function() {
  stopifnot(requireNamespace('Seurat',quietly=TRUE),requireNamespace('DESeq2',quietly=TRUE))
  reviewed_md5<-readLines(file.path(root,'data/core1_corrected_reference/reviewed_partition.md5'))
  stopifnot(identical(unname(tools::md5sum(file.path(out,'tables/pc_cluster_assignments.csv'))),reviewed_md5))
  obj<-readRDS(checkpoint('04_branch_evidence'))
  ann<-utils::read.csv(file.path(root,'data/core1_corrected_reference/cluster_annotation_review.csv'))
  genes<-readLines(file.path(out,'diagnostics/common_gene_features.txt'))
  x<-SeuratObject::LayerData(obj,layer='counts')[genes,,drop=FALSE]
  stopifnot(!any(grepl('^[^:]+:[0-9]+-[0-9]+$',rownames(x))),all(x@x==round(x@x)))
  allres<-list(); pbout<-list(); samplediag<-list()
  for(b in c('pca','harmony')) {
    a<-ann[ann$branch==b,]; cl<-as.character(obj[[paste0(b,'_cluster')]][,1])
    stopifnot(!anyDuplicated(a$cluster),setequal(cl,as.character(a$cluster)))
    pop<-a$cell_population[match(cl,a$cluster)]
    obj[[paste0(b,'_population')]]<-pop
    # Mappings are validated against saved cell-level partitions as well as IDs.
    expected<-utils::read.csv(file.path(out,'tables/pc_cluster_assignments.csv'))
    stopifnot(identical(colnames(obj),expected$cell),identical(cl,as.character(expected[[paste0(b,'_25')]])))
    tb<-table(pop,factor(obj$sample,levels=samples$sample))
    eligibility<-data.frame(cell_population=rownames(tb),min_cells=apply(tb,1,min),
      eligible25=apply(tb,1,min)>=25,eligible50=apply(tb,1,min)>=50,eligible100=apply(tb,1,min)>=100)
    write_table(eligibility,paste0(b,'_pseudobulk_eligibility'))
    tc<-as.data.frame.matrix(tb);tc$cell_population<-rownames(tc);write_table(tc,paste0(b,'_population_sample_counts'))
    plot_save(Seurat::DimPlot(obj,reduction=paste0('umap_',b),group.by=paste0(b,'_population'),label=TRUE,repel=TRUE),paste0('05_',b,'_population_umap'),w=13,h=8)
    for(p in eligibility$cell_population[eligibility$eligible50]) {
      message('DESeq2 ',b,' ',p)
      cells<-which(pop==p)
      # Sparse indicator multiplication avoids parsing Seurat-generated names.
      indicator<-Matrix::sparseMatrix(i=seq_along(cells),j=match(obj$sample[cells],samples$sample),x=1,dims=c(length(cells),nrow(samples)))
      pb<-x[,cells,drop=FALSE] %*% indicator
      colnames(pb)<-samples$sample
      stopifnot(all(Matrix::colSums(pb)>0),all(pb@x==round(pb@x)),all(tb[p,]>=50))
      meta<-samples;rownames(meta)<-meta$sample
      meta$condition<-factor(meta$condition,levels=c('control','stress'))
      keep<-Matrix::rowSums(pb)>=10
      dds<-DESeq2::DESeqDataSetFromMatrix(as.matrix(pb[keep,,drop=FALSE]),meta,~condition)
      dds<-DESeq2::DESeq(dds,quiet=TRUE)
      res<-as.data.frame(DESeq2::results(dds,contrast=c('condition','stress','control'),alpha=.05))
      res$gene<-rownames(res);res$cell_population<-p;res$branch<-b
      res$significant<-!is.na(res$padj)&res$padj<.05&abs(res$log2FoldChange)>=1
      key<-paste(b,p,sep='__');allres[[key]]<-res;pbout[[key]]<-pb
      norm<-DESeq2::counts(dds,normalized=TRUE)
      cors<-stats::cor(log2(norm+1),method='spearman')
      samplediag[[key]]<-data.frame(branch=b,population=p,sample=samples$sample,n_cells=as.numeric(tb[p,]),
        library_size=as.numeric(Matrix::colSums(pb)),size_factor=DESeq2::sizeFactors(dds))
      write_table(data.frame(sample=rownames(cors),cors),paste0('pb_sample_correlation_',key))
      write_table(res,paste0('DE_',key))
    }
  }
  de<-dplyr::bind_rows(allres)
  write_table(de,'pseudobulk_DE_all_branches')
  sm<-de |> dplyr::group_by(branch,cell_population) |> dplyr::summarise(
    fitted_genes=dplyr::n(),valid_pvalue=sum(!is.na(pvalue)),finite_padj=sum(!is.na(padj)),
    DEGs=sum(significant),up=sum(significant & log2FoldChange>=1),down=sum(significant & log2FoldChange<=-1),.groups='drop')
  write_table(sm,'pseudobulk_DE_summary')
  write_table(dplyr::bind_rows(samplediag),'pseudobulk_sample_diagnostics')
  md<-obj[[]];md$cell<-rownames(md);write_table(md,'annotated_cell_metadata')
  saveRDS(list(counts=pbout,results=allres,samples=samples,common_features=genes),checkpoint('05_pseudobulk'),compress=FALSE)
  saveRDS(obj,checkpoint('05_annotated'),compress=FALSE)
  write_table(ann,'cluster_annotation_review')
  session_save('05_pseudobulk')
}
