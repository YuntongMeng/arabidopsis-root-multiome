# Evidence stage, BEFORE final annotation and DE. Atlas signatures are published
# top-50 genes (S2K); canonical validation uses separately curated S2F/G markers.
# Both derive from Shahan 2022 and are supportive, not independent ground truth.
run_branch_evidence <- function(npcs) {
  stopifnot(requireNamespace('Seurat',quietly=TRUE),requireNamespace('Matrix',quietly=TRUE))
  obj <- readRDS(checkpoint('03_embeddings'))
  asg <- utils::read.csv(file.path(out,'tables/pc_cluster_assignments.csv'),check.names=FALSE)
  stopifnot(identical(asg$cell,colnames(obj)))
  atlas <- utils::read.csv(file.path(root,'data/core1_corrected_reference/atlas_top50_signatures.csv'),check.names=FALSE)
  panels <- split(atlas$gene,atlas[['cell type group']]); panels <- lapply(panels,base::intersect,rownames(obj))
  obj <- Seurat::AddModuleScore(obj,features=panels,name='Atlas',seed=42,ctrl=100)
  score_cols <- paste0('Atlas',seq_along(panels))
  write_table(data.frame(score_column=score_cols,signature=names(panels),n_genes=lengths(panels)),'atlas_score_map')
  scores <- as.matrix(obj[[]][,score_cols]);colnames(scores)<-names(panels)
  top <- max.col(scores,ties.method='first'); second <- apply(scores,1,function(v) sort(v,decreasing=TRUE)[2])
  obj$atlas_top <- colnames(scores)[top];obj$atlas_score<-scores[cbind(seq_len(nrow(scores)),top)]
  obj$atlas_margin<-obj$atlas_score-second
  broad <- function(v) {v[grepl('Pericycle',v)]<-'Pericycle';v[grepl('phloem|Phloem',v)]<-'Phloem';v[grepl('xylem',v)]<-'Xylem';v}
  obj$atlas_broad <- broad(obj$atlas_top)
  write_table(data.frame(cell=colnames(obj),scores,atlas_top=obj$atlas_top,atlas_broad=obj$atlas_broad,score=obj$atlas_score,margin=obj$atlas_margin),'cell_signature_scores')
  canon <- utils::read.csv(file.path(root,'data/core1_corrected_reference/canonical_markers_provenance.csv'))
  genes <- unique(canon$Locus);genes<-genes[genes %in% rownames(obj)]
  symbols <- canon$Gene[match(genes,canon$Locus)]
  stress_genes <- c('AT5G52310','AT5G52300','AT5G66400','AT5G05410')
  rawdata <- SeuratObject::LayerData(obj,layer='data')
  obj$stress_marker_mean <- Matrix::colMeans(rawdata[stress_genes,,drop=FALSE])
  cycles<-utils::read.csv(file.path(root,'data/core1_corrected_reference/cell_cycle_markers.csv'))
  cc <- unique(stats::na.omit(cycles[[1]]));cc<-intersect(cc,rownames(obj))
  obj<-Seurat::AddModuleScore(obj,features=list(cc),name='CellCycle',seed=42)
  mix <- list();global <- list();clusters <- list();canontabs<-list()
  for(red in c('pca','harmony')) {
    cl<-asg[[paste0(red,'_',npcs)]]; obj[[paste0(red,'_cluster')]]<-cl
    SeuratObject::Idents(obj)<-cl
    obj<-Seurat::RunUMAP(obj,reduction=red,dims=seq_len(npcs),seed.use=42,
      reduction.name=paste0('umap_',red),reduction.key=paste0('UMAP',red,'_'),verbose=FALSE)
    for(g in c(paste0(red,'_cluster'),'sample','condition','atlas_broad')) plot_save(Seurat::DimPlot(obj,reduction=paste0('umap_',red),group.by=g,label=g==paste0(red,'_cluster'),raster=FALSE),paste0('04_',red,'_',g))
    plot_save(Seurat::DotPlot(obj,features=genes,group.by=paste0(red,'_cluster'))+Seurat::RotatedAxis()+ggplot2::scale_x_discrete(labels=setNames(symbols,genes)),paste0('04_',red,'_canonical_DotPlot'),w=18,h=9)
    means<-stats::aggregate(scores,list(cluster=cl),mean)
    long<-tidyr::pivot_longer(means,-cluster,names_to='signature',values_to='mean_score')
    write_table(long,paste0(red,'_cluster_signature_scores'))
    ranks<-long |> dplyr::group_by(cluster) |> dplyr::arrange(dplyr::desc(mean_score),.by_group=TRUE) |> dplyr::slice_head(n=3) |> dplyr::ungroup()
    write_table(ranks,paste0(red,'_cluster_top3_signatures'))
    ctab<-as.data.frame.matrix(table(cluster=cl,sample=obj$sample));ctab$cluster<-rownames(ctab)
    write_table(ctab,paste0(red,'_cluster_sample_counts'))
    can<-Seurat::DotPlot(obj,features=genes,group.by=paste0(red,'_cluster'))$data
    can$symbol<-symbols[match(can$features.plot,genes)];write_table(can,paste0(red,'_canonical_expression'))
    # Local metrics are descriptive. No cell-level inferential p-values.
    e<-SeuratObject::Embeddings(obj,red)[,seq_len(npcs),drop=FALSE]
    nn<-RANN::nn2(e,k=31)$nn.idx[,-1]
    same_cond<-rowMeans(matrix(obj$condition[nn],nrow=nrow(nn))==obj$condition)
    smooth<-rowMeans(matrix(obj$stress_marker_mean[nn],nrow=nrow(nn)))
    type_purity<-rowMeans(matrix(obj$atlas_broad[nn],nrow=nrow(nn))==obj$atlas_broad)
    global[[red]]<-data.frame(branch=red,PCs=npcs,clusters=length(unique(cl)),
      mean_neighbor_same_condition=mean(same_cond),atlas_broad_neighbor_consistency=mean(type_purity),
      stress_marker_neighbor_spearman=stats::cor(obj$stress_marker_mean,smooth,method='spearman'))
    md<-obj[[]]
    for(cond in unique(md$condition)) for(ty in unique(md$atlas_broad)) {
      ix<-which(md$condition==cond & md$atlas_broad==ty & md$atlas_score>.1 & md$atlas_margin>.1)
      ss<-table(factor(md$sample[ix],levels=samples$sample[samples$condition==cond]))
      if(length(ix)<60 || any(ss<30)) next
      ni<-RANN::nn2(e[ix,,drop=FALSE],k=31)$nn.idx[,-1]
      obs<-mean(matrix(md$sample[ix][ni],nrow=nrow(ni))!=md$sample[ix])
      expected<-1-sum((as.numeric(ss)/sum(ss))^2)
      mix[[paste(red,cond,ty)]]<-data.frame(branch=red,condition=cond,signature=ty,n=length(ix),
        other_replicate_fraction=obs,expected_random=expected,mixing_ratio=obs/expected)
    }
    message('Strict markers ',red)
    markers<-Seurat::FindAllMarkers(obj,only.pos=TRUE,min.pct=.25,logfc.threshold=.25,
      test.use='wilcox',return.thresh=.01,verbose=FALSE)
    stopifnot(!any(grepl('^[^:]+:[0-9]+-[0-9]+$',markers$gene)))
    markers$SYMBOL<-unname(AnnotationDbi::mapIds(org.At.tair.db::org.At.tair.db,keys=unique(markers$gene),keytype='TAIR',column='SYMBOL',multiVals='first')[markers$gene])
    write_table(markers,paste0(red,'_strict_markers_all'))
    tops<-markers |> dplyr::group_by(cluster) |> dplyr::slice_max(avg_log2FC,n=25,with_ties=FALSE) |> dplyr::ungroup()
    write_table(tops,paste0(red,'_top25_markers'))
  }
  write_table(dplyr::bind_rows(global),'branch_global_metrics')
  write_table(dplyr::bind_rows(mix),'replicate_mixing_within_condition_signature')
  md<-obj[[]];md$cell<-rownames(md);write_table(md,'branch_cell_metadata')
  # Counts unchanged by integration; same assay used by both branches.
  saveRDS(obj,checkpoint('04_branch_evidence'),compress=FALSE)
  session_save('04_evidence')
}
