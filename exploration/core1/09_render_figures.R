# Re-render all Core 1 figures from saved analysis outputs. No clustering, UMAP,
# marker tests or DE are recalculated here. Source defines the function only.
# Use the main runner's figures stage, or source 00_common and call this function.
run_render_figures <- function() {
  stopifnot(requireNamespace('Matrix',quietly=TRUE),requireNamespace('Seurat',quietly=TRUE))
  message('Rendering input QC and doublet figures')
  qc <- readRDS(checkpoint('01_qc'))
  long <- tidyr::pivot_longer(qc,cols=dplyr::all_of(c('nCount_RNA','nFeature_RNA')),names_to='metric',values_to='value')
  plot_save(ggplot2::ggplot(long,ggplot2::aes(value,colour=sample))+ggplot2::geom_density()+ggplot2::scale_x_log10()+ggplot2::facet_wrap(~metric,scales='free')+ggplot2::theme_bw(),'01_gene_only_qc_density')
  plot_save(ggplot2::ggplot(qc,ggplot2::aes(nCount_RNA,nFeature_RNA,colour=condition))+ggplot2::geom_point(size=.3,alpha=.2)+ggplot2::scale_x_log10()+ggplot2::scale_y_log10()+ggplot2::facet_wrap(~sample)+ggplot2::theme_bw(),'01_gene_only_qc_scatter')
  md <- utils::read.csv(file.path(out,'tables/qc_doublet_calls_all.csv'))
  plot_save(ggplot2::ggplot(md[md$qc_pass,],ggplot2::aes(scDblFinder_score,fill=scDblFinder_class))+ggplot2::geom_histogram(bins=50)+ggplot2::facet_wrap(~sample)+ggplot2::theme_bw(),'02_doublet_scores')
  rm(qc,long,md);gc()
  message('Rendering dimension spread and clustering stability')
  obj <- readRDS(checkpoint('03_embeddings'))
  plot_save(Seurat::ElbowPlot(obj,ndims=50),'03_pca_elbow')
  spread <- dplyr::bind_rows(lapply(c('pca','harmony'),function(b) {
    scores <- SeuratObject::Embeddings(obj,reduction=b)
    data.frame(branch=b,dimension=seq_len(ncol(scores)),coordinate_sd=apply(scores,2,stats::sd))
  }))
  write_table(spread,'embedding_coordinate_sd')
  p <- ggplot2::ggplot(spread,ggplot2::aes(dimension,coordinate_sd,color=branch))+ggplot2::geom_line(linewidth=.8)+
    ggplot2::scale_color_manual(values=c(pca='#31759B',harmony='#BE6B37'))+
    ggplot2::labs(title='PCA and Harmony coordinate spread',x='Dimension index',y='Standard deviation of cell coordinates',color=NULL,
      caption='Harmony corrects PCA coordinates; its dimensions are not re-ranked by variance.
This is coordinate spread, not two explained-variance curves or a stand-alone PC-selection rule.')+
    ggplot2::theme_bw()+ggplot2::theme(legend.position='top',plot.caption=ggplot2::element_text(hjust=0))
  plot_save(p,'03_embedding_coordinate_sd',w=10,h=6)
  rm(obj);gc()
  supplementary <- file.path(out,'diagnostics','figures')
  dir.create(supplementary,recursive=TRUE,showWarnings=FALSE)
  diagnostic_plot_save <- function(p,name,w=10,h=6) ggplot2::ggsave(file.path(supplementary,paste0(name,'.png')),p,width=w,height=h,dpi=160,bg='white')
  ari <- utils::read.csv(file.path(out,'tables/pc_ARI_long.csv'))
  dims <- c(10,15,20,25,30,40,50)
  adjacent <- dplyr::bind_rows(lapply(c('pca','harmony'),function(b) {
    dplyr::bind_rows(lapply(seq_len(length(dims)-1),function(i) {
      a <- paste0(b,'_',dims[i]);bb <- paste0(b,'_',dims[i+1])
      data.frame(branch=b,pair=paste(dims[i],dims[i+1],sep='–'),ARI=ari$ARI[ari$a==a & ari$b==bb])
    }))
  }))
  adjacent$pair <- factor(adjacent$pair,levels=paste(head(dims,-1),tail(dims,-1),sep='–'))
  write_table(adjacent,'pc_ARI_adjacent')
  p <- ggplot2::ggplot(adjacent,ggplot2::aes(pair,ARI,color=branch,group=branch))+
    ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=3)+
    ggplot2::scale_color_manual(values=c(pca='#31759B',harmony='#BE6B37'))+ggplot2::scale_y_continuous(limits=c(0,1))+
    ggplot2::labs(title='Clustering stability between adjacent dimension choices',x='Dimensions compared',y='Adjusted Rand Index',color=NULL,
      caption='Same retained nuclei, clustering resolution and seed; only dimensions change.
ARI measures partition agreement, not biological accuracy. Both 20–30 and 40–50 deserve inspection.')+
    ggplot2::theme_bw()+ggplot2::theme(legend.position='top',plot.caption=ggplot2::element_text(hjust=0))
  plot_save(p,'03_ARI_adjacent',w=10,h=6)
  z <- ari[(ari$a=='pca_25' & grepl('^pca_',ari$b)) | (ari$a=='harmony_25' & grepl('^harmony_',ari$b)),]
  z$branch <- ifelse(grepl('^pca',z$a),'Unintegrated PCA','Harmony')
  z$dimensions <- as.integer(sub('.*_','',z$b))
  z <- z[z$dimensions!=25,];z<-z[order(z$branch,z$dimensions),]
  p <- ggplot2::ggplot(z,ggplot2::aes(dimensions,ARI,color=branch,group=branch))+
    ggplot2::annotate('rect',xmin=20,xmax=30,ymin=0,ymax=1,fill='#EAF0F5',alpha=.8)+
    ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=3)+
    ggplot2::scale_color_manual(values=c('Unintegrated PCA'='#31759B','Harmony'='#BE6B37'))+
    ggplot2::scale_x_continuous(breaks=c(10,15,20,25,30,40,50))+
    ggplot2::scale_y_continuous(limits=c(0,1))+
    ggplot2::labs(title='Clustering agreement relative to 25 dimensions',x='Dimensions used for clustering',y='Adjusted Rand Index',color=NULL,
      caption='Each branch uses its own 25-dimension reference. Self-comparison (ARI = 1) omitted.
Agreement is not biological accuracy; this plot alone cannot select an optimal dimension count.')+
    ggplot2::theme_bw()+ggplot2::theme(legend.position='top',plot.caption=ggplot2::element_text(hjust=0))
  diagnostic_plot_save(p,'03_ARI_reference25',w=10,h=6)
  keys <- c(paste0('pca_',c(10,15,20,25,30,40,50)),paste0('harmony_',c(10,15,20,25,30,40,50)))
  ari$a <- factor(ari$a,levels=rev(keys));ari$b<-factor(ari$b,levels=keys)
  p <- ggplot2::ggplot(ari,ggplot2::aes(b,a,fill=ARI))+ggplot2::geom_tile()+
    ggplot2::geom_text(ggplot2::aes(label=sprintf('%.2f',ARI),color=ifelse(ARI>.65,'white','#222222')),size=2.8)+ggplot2::scale_color_identity()+
    ggplot2::geom_vline(xintercept=7.5,color='#666666',linewidth=.4)+ggplot2::geom_hline(yintercept=7.5,color='#666666',linewidth=.4)+
    ggplot2::scale_fill_gradient(low='white',high='#28668B',limits=c(0,1))+
    ggplot2::labs(title='All-pairs clustering agreement',subtitle='Includes within-branch and between-branch comparisons',x=NULL,y=NULL)+
    ggplot2::theme_minimal()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=50,hjust=1),panel.grid=ggplot2::element_blank())
  diagnostic_plot_save(p,'03_ARI_all_pairs',w=11,h=10)
  message('Rendering both embedding branches and canonical markers')
  obj <- readRDS(checkpoint('04_branch_evidence'))
  canon <- utils::read.csv(file.path(root,'data/core1_corrected_reference/canonical_markers_provenance.csv'))
  genes <- unique(canon$Locus);genes<-genes[genes%in%rownames(obj)]
  symbols <- canon$Gene[match(genes,canon$Locus)]
  for(red in c('pca','harmony')) {
    for(g in c(paste0(red,'_cluster'),'sample','condition','atlas_broad'))
      plot_save(Seurat::DimPlot(obj,reduction=paste0('umap_',red),group.by=g,label=g==paste0(red,'_cluster'),raster=FALSE),paste0('04_',red,'_',g))
    plot_save(Seurat::DotPlot(obj,features=genes,group.by=paste0(red,'_cluster'))+Seurat::RotatedAxis()+ggplot2::scale_x_discrete(labels=setNames(symbols,genes)),paste0('04_',red,'_canonical_DotPlot'),w=18,h=9)
  }
  rm(obj);gc()
  message('Rendering population labels and biological-sample composition')
  obj <- readRDS(checkpoint('05_annotated'))
  for(b in c('pca','harmony'))
    plot_save(Seurat::DimPlot(obj,reduction=paste0('umap_',b),group.by=paste0(b,'_population'),label=TRUE,repel=TRUE),paste0('05_',b,'_population_umap'),w=13,h=8)
  rm(obj);gc()
  cells <- utils::read.csv(file.path(out,'tables/annotated_cell_metadata.csv'))
  ann <- utils::read.csv(file.path(root,'data/core1_corrected_reference/cluster_annotation_review.csv'))
  counts <- dplyr::bind_rows(lapply(c('pca','harmony'),function(b) {
    tb <- as.data.frame(table(population=cells[[paste0(b,'_population')]],sample=factor(cells$sample,levels=samples$sample)))
    tb$branch <- b
    low <- unique(ann$cell_population[ann$branch==b & ann$confidence=='low'])
    states <- 'Dividing_meristematic'
    tb$label <- paste0(tb$population,ifelse(tb$population%in%low,' [low confidence]',ifelse(tb$population%in%states,' [state]','')))
    tb$eligible <- tb$population %in% unique(tb$population[vapply(tb$population,function(p)min(tb$Freq[tb$population==p])>=50,logical(1))])
    tb
  }))
  p <- ggplot2::ggplot(counts,ggplot2::aes(sample,label,fill=log10(Freq+1)))+
    ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=Freq),size=3)+
    ggplot2::facet_wrap(~branch,ncol=1,scales='free_y')+
    ggplot2::scale_fill_gradient(low='white',high='#9ABDCF',name='log10(n + 1)')+
    ggplot2::labs(title='Population counts across biological samples',x=NULL,y=NULL,
      caption='All populations shown, including exclusions. DE requires >=50 nuclei in every sample.
Identity labels are reviewed broad annotations, not automatic atlas predictions; low-confidence and state labels are marked.')+
    ggplot2::theme_bw()+ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0),axis.text.x=ggplot2::element_text(angle=20,hjust=1))
  plot_save(p,'05_population_sample_counts',w=12,h=12)
  message('Rendering P1 and GO context')
  en <- utils::read.csv(file.path(out,'tables/P1_Fisher_enrichment_all.csv'))
  complete <- utils::read.csv(file.path(out,'tables/GO_population_complete.csv'))
  f <- dplyr::filter(en,universe=='finite_padj_primary',test=='interaction_vs_background')
  plot_save(ggplot2::ggplot(f,ggplot2::aes(cell_population,-log10(BH),fill=branch))+ggplot2::geom_col(position='dodge')+ggplot2::geom_hline(yintercept=-log10(.05),linetype=2)+ggplot2::coord_flip()+ggplot2::theme_bw()+ggplot2::labs(y='-log10 BH-adjusted Fisher p',x=NULL),'06_P1_enrichment')
  plot_save(ggplot2::ggplot(complete,ggplot2::aes(cell_population,Description,fill=DE_genes))+ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=DE_genes),size=3)+ggplot2::facet_wrap(~branch,ncol=1,scales='free_x')+ggplot2::scale_fill_gradient(low='white',high='#416a9b')+ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=40,hjust=1))+ggplot2::labs(x=NULL,y=NULL),'06_GO_context',w=13,h=9)
  session_save('09_figure_render')
  render_P1_control_figure(en)
  message('22 main figures and 2 supplementary ARI figures rendered. No analysis stages recomputed.')
}

# Show the alternative background explicitly: significance relative to ordinary
# genes does not establish extra enrichment over drought-responsive controls.
render_P1_control_figure <- function(en) {
  d <- dplyr::filter(en,universe=='finite_padj_primary',test%in%c('interaction_vs_background','interaction_vs_drought_only'))
  d$background <- ifelse(d$test=='interaction_vs_background','Other testable genes','P1 drought DEGs excluding candidates')
  d$significance <- ifelse(d$BH<.05,'BH < 0.05','BH >= 0.05')
  pos <- ggplot2::position_dodge(width=.6)
  p <- ggplot2::ggplot(d,ggplot2::aes(odds_ratio,cell_population,color=background,shape=significance,group=background))+
    ggplot2::geom_vline(xintercept=1,linetype=2,color='#777777')+
    ggplot2::geom_errorbar(ggplot2::aes(xmin=CI_low,xmax=CI_high),orientation='y',position=pos,width=.15)+
    ggplot2::geom_point(position=pos,size=2.8)+ggplot2::scale_x_log10()+
    ggplot2::scale_color_manual(values=c('Other testable genes'='#31759B','P1 drought DEGs excluding candidates'='#BE6B37'))+
    ggplot2::scale_shape_manual(values=c('BH < 0.05'=16,'BH >= 0.05'=1))+
    ggplot2::facet_wrap(~branch,ncol=1,scales='free_y')+
    ggplot2::labs(title='P1 candidate enrichment depends on the comparison background',x='Fisher odds ratio (log scale; two-sided 95% confidence interval)',y=NULL,color='Background',shape=NULL,
      caption='Enrichment tests are one-sided; BH adjustment is within branch and test family.
No population passes BH < 0.05 versus the drought-DEG control; this is not proof of equal effects.
Low-confidence populations remain exploratory. Gene sets are not matched for expression or detection power.')+
    ggplot2::theme_bw()+ggplot2::theme(legend.position='bottom',legend.box='vertical',plot.caption=ggplot2::element_text(hjust=0),plot.title=ggplot2::element_text(size=14))
  plot_save(p,'06_P1_background_control',w=12,h=12)
}
