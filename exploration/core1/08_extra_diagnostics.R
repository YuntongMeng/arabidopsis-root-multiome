# Additional evidence checks; no analysis starts on source.
run_extra_diagnostics <- function() {
  stopifnot(requireNamespace('Seurat',quietly=TRUE))
  old<-readRDS(file.path(root,'archive/core1_initial/results/reproducibility/p3_core1_checkpoint_2026-09-30.rds'))
  hv<-SeuratObject::VariableFeatures(old);coordinate<-grepl('^[^:]+:[0-9]+-[0-9]+$',hv)
  newhv<-readLines(file.path(out,'diagnostics/embedding_HVGs.txt'))
  write_table(data.frame(old_HVGs=length(hv),old_peak_HVGs=sum(coordinate),new_HVGs=length(newhv),shared_gene_HVGs=length(intersect(hv,newhv))),'old_new_HVG_audit')
  rm(old);gc()
  obj<-readRDS(checkpoint('05_annotated')); md<-obj[[]]; x<-SeuratObject::LayerData(obj,layer='counts')
  sg<-c('AT5G52310','AT5G52300','AT5G66400','AT5G05410')
  d<-dplyr::bind_rows(lapply(sg,function(g) dplyr::bind_rows(lapply(samples$sample,function(s) data.frame(gene=g,sample=s,n_cells=sum(md$sample==s),fraction_detected=mean(x[g,md$sample==s]>0),raw_sum=sum(x[g,md$sample==s]))))))
  write_table(d,'stress_diagnostic_marker_detection')
  # Cell lineage evidence consistency by sample: exported for uncertain calls.
  for(b in c('pca','harmony')) {
    groups<-paste(md[[paste0(b,'_population')]],md$sample,sep=' | ')
    scores<-md[,grep('^Atlas[0-9]+$',names(md)),drop=FALSE]
    s<-stats::aggregate(scores,list(population_sample=groups),mean)
    write_table(s,paste0(b,'_population_sample_atlas_scores'))
  }
  # Raw assay integrity: sample totals over all retained nuclei exactly match
  # pre-embedding singlet input totals, despite two reductions and annotations.
  xs<-readRDS(checkpoint('02_singlet_counts'))
  z<-data.frame(sample=samples$sample,pre_embedding_sum=vapply(xs,sum,numeric(1)),
    final_raw_sum=vapply(samples$sample,function(s)sum(x[,md$sample==s]),numeric(1)))
  z$identical<-z$pre_embedding_sum==z$final_raw_sum;stopifnot(all(z$identical));write_table(z,'raw_count_preservation_audit')
}
