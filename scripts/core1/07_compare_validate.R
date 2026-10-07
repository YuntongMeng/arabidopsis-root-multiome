read_frozen_table <- function(name) {
  p <- file.path(root,'results/core1_frozen',paste0(name,'.csv'))
  if(file.exists(p)) return(utils::read.csv(p))
  if(file.exists(paste0(p,'.gz'))) return(utils::read.csv(gzfile(paste0(p,'.gz'))))
  stop('Missing archived comparison table: ',name)
}

run_compare_validate <- function() {
  old<-read_frozen_table('cell_metadata')
  new<-utils::read.csv(file.path(out,'tables/annotated_cell_metadata.csv'))
  de<-utils::read.csv(file.path(out,'tables/pseudobulk_DE_all_branches.csv'))
  oldde<-read_frozen_table('de_all')
  shared<-intersect(old$cell,new$cell);oi<-old[match(shared,old$cell),];ni<-new[match(shared,new$cell),]
  qc<-utils::read.csv(file.path(out,'tables/qc_doublet_calls_all.csv'))
  retention<-data.frame(old_retained=nrow(old),corrected_retained=nrow(new),both=length(shared),
    old_only=length(setdiff(old$cell,new$cell)),new_only=length(setdiff(new$cell,old$cell)))
  write_table(retention,'old_new_cell_retention')
  ari<-data.frame(comparison=c('old_vs_corrected_pca','old_vs_corrected_harmony','corrected_pca_vs_harmony'),
     n_cells=c(length(shared),length(shared),nrow(new)),
     ARI=c(mclust::adjustedRandIndex(oi$seurat_clusters,ni$pca_cluster),mclust::adjustedRandIndex(oi$seurat_clusters,ni$harmony_cluster),mclust::adjustedRandIndex(new$pca_cluster,new$harmony_cluster)))
  write_table(ari,'old_new_partition_ARI')
  for(b in c('pca','harmony')) {
    tb<-as.data.frame(table(old_population=oi$cell_population,new_population=ni[[paste0(b,'_population')]]))
    tb<-tb[tb$Freq>0,];write_table(tb,paste0('old_to_',b,'_annotation_transition'))
  }
  # Branch DE concordance compares the same named broad populations, but cell
  # membership may differ. Shared significance is not an independent validation.
  comparisons<-list()
  for(p in intersect(de$cell_population[de$branch=='pca'],de$cell_population[de$branch=='harmony'])) {
    a<-de[de$branch=='pca'&de$cell_population==p,];b<-de[de$branch=='harmony'&de$cell_population==p,]
    z<-merge(a,b,by='gene');sig1<-a$gene[a$significant];sig2<-b$gene[b$significant]
    both<-z$significant.x&z$significant.y
    ca<-new$cell[new$pca_population==p];cb<-new$cell[new$harmony_population==p]
    comparisons[[p]]<-data.frame(population=p,common_fitted=nrow(z),logFC_spearman=stats::cor(z$log2FoldChange.x,z$log2FoldChange.y,method='spearman'),
      pca_DEGs=length(sig1),harmony_DEGs=length(sig2),shared_DEGs=length(intersect(sig1,sig2)),DEG_Jaccard=length(intersect(sig1,sig2))/length(union(sig1,sig2)),
      shared_DEG_sign_agreement=mean(sign(z$log2FoldChange.x[both])==sign(z$log2FoldChange.y[both])),cell_Jaccard=length(intersect(ca,cb))/length(union(ca,cb)))
  }
  write_table(dplyr::bind_rows(comparisons),'branch_DE_concordance')
  oldcomp<-list()
  for(p in unique(de$cell_population[de$branch=='pca'])) {
    ca<-new$cell[new$pca_population==p]
    for(op in unique(oldde$cell_population)) {
      cb<-old$cell[old$cell_population==op];inter<-intersect(ca,cb)
      if(length(inter)<.1*length(ca))next
      a<-de[de$branch=='pca'&de$cell_population==p,];b<-oldde[oldde$cell_population==op,];
      b$significant<-!is.na(b$padj)&b$padj<.05&abs(b$log2FoldChange)>=1
      z<-merge(a,b,by='gene');sig1<-a$gene[a$significant];sig2<-b$gene[b$significant]
      oldcomp[[paste(p,op)]]<-data.frame(new_population=p,old_population=op,shared_cells=length(inter),fraction_of_new=length(inter)/length(ca),
        cell_Jaccard=length(inter)/length(union(ca,cb)),old_DEGs=length(sig2),new_DEGs=length(sig1),
        common_fitted=nrow(z),logFC_spearman=stats::cor(z$log2FoldChange.x,z$log2FoldChange.y,method='spearman'),
        shared_DEGs=length(intersect(sig1,sig2)),DEG_Jaccard=length(intersect(sig1,sig2))/length(union(sig1,sig2)))
    }
  }
  write_table(dplyr::bind_rows(oldcomp),'old_new_DE_comparison_membership_qualified')
  # Biological invariants only; no fixed cell/cluster/population/DE answer.
  aud<-utils::read.csv(file.path(out,'tables/input_matrix_audit.csv'))
  checks<-data.frame(check=c('no_peak_IDs_after_input','nonempty_retained','unique_retained_cell_IDs','retained_are_QC_singlets','unique_gene_population_branch','DE_only_common_gene_universe','no_peak_markers','primary_eligibility_satisfied'),pass=TRUE)
  checks$pass[1]<-all(aud$coordinate_ids_after==0)
  checks$pass[2]<-nrow(new)>0
  checks$pass[3]<-!anyDuplicated(new$cell)
  checks$pass[4]<-all(qc$retained[match(new$cell,qc$cell)])
  checks$pass[5]<-!anyDuplicated(de[c('branch','cell_population','gene')])
  checks$pass[6]<-all(de$gene%in%readLines(file.path(out,'diagnostics/common_gene_features.txt')))
  checks$pass[7]<-all(vapply(c('pca','harmony'),function(b){m<-utils::read.csv(file.path(out,'tables',paste0(b,'_strict_markers_all.csv')));!any(grepl('^[^:]+:[0-9]+-[0-9]+$',m$gene))},logical(1)))
  d<-utils::read.csv(file.path(out,'tables/pseudobulk_sample_diagnostics.csv'));checks$pass[8]<-all(d$n_cells>=50)
  stopifnot(all(checks$pass));write_table(checks,'structural_invariants')
  # New input audit also reproduces the old barcode filter solely as a check.
  stopifnot(setequal(qc$cell[qc$old_qc_pass],old$cell))
  writeLines('Old mixed-input QC selection matches archived barcodes exactly; used for comparison only.',file.path(out,'diagnostics/old_QC_reconstruction.txt'))
  session_save('07_compare')
}
