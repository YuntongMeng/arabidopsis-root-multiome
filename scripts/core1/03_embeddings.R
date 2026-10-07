# Decision 3: identical gene-only singlets and RNA expression in both branches.
# Harmony by sample can remove condition-related structure because each sample
# belongs to one condition. Compare it as a sensitivity branch before selection.
run_embeddings <- function() {
  stopifnot(requireNamespace("Matrix",quietly=TRUE),requireNamespace("Seurat",quietly=TRUE))
  xs <- readRDS(checkpoint('02_singlet_counts')); md <- readRDS(checkpoint('02_singlet_metadata'))
  # Confirm common gene sets and align by TAIR ID before cbind.
  # Load Matrix before reading S4 matrices to ensure dimension methods dispatch.
  universe <- rownames(xs[[1]])
  stopifnot(all(vapply(xs,function(x) setequal(rownames(x),universe),logical(1))))
  xs <- lapply(xs,function(x) x[universe,,drop=FALSE])
  stopifnot(all(vapply(xs,function(x) identical(rownames(x),universe),logical(1))))
  x <- do.call(cbind,xs); md <- md[match(colnames(x),md$cell),]; rownames(md)<-md$cell
  obj <- SeuratObject::CreateSeuratObject(x,meta.data=md,project='P3_gene_only')
  rm(xs,x);gc()
  obj <- Seurat::NormalizeData(obj,verbose=FALSE)
  obj <- Seurat::FindVariableFeatures(obj,nfeatures=2000,verbose=FALSE)
  # Organelle genes retained in raw counts/DE, excluded from HVG embedding only.
  hv <- SeuratObject::VariableFeatures(obj)
  hv <- hv[!grepl('^AT[MC]G',hv)];SeuratObject::VariableFeatures(obj)<-hv
  writeLines(hv,file.path(out,'diagnostics/embedding_HVGs.txt'))
  stopifnot(!any(grepl('^[^:]+:[0-9]+-[0-9]+$',hv)))
  obj <- Seurat::ScaleData(obj,features=hv,verbose=FALSE)
  obj <- Seurat::RunPCA(obj,features=hv,npcs=50,seed.use=42,verbose=FALSE)
  plot_save(Seurat::ElbowPlot(obj,ndims=50),'03_pca_elbow')
  write_table(data.frame(PC=1:50,stdev=SeuratObject::Stdev(obj,'pca')),'pca_stdev')
  set.seed(42)
  obj <- harmony::RunHarmony(obj,group.by.vars='sample',reduction.use='pca',dims.use=1:50,
    reduction.save='harmony',theta=2,lambda=1,max_iter=20,ncores=2,project.dim=FALSE)
  assignments <- data.frame(cell=colnames(obj)); grid <- list()
  for (red in c('pca','harmony')) for (np in c(10,15,20,25,30,40,50)) {
    message('Clustering ',red,' PC ',np)
    cand <- Seurat::FindNeighbors(obj,reduction=red,dims=seq_len(np),verbose=FALSE)
    cand <- Seurat::FindClusters(cand,resolution=.5,random.seed=42,verbose=FALSE)
    key <- paste0(red,'_',np); assignments[[key]]<-as.character(SeuratObject::Idents(cand))
    grid[[key]] <- data.frame(reduction=red,PCs=np,clusters=length(unique(assignments[[key]])))
  }
  write_table(assignments,'pc_cluster_assignments')
  write_table(dplyr::bind_rows(grid),'pc_cluster_numbers')
  ari <- dplyr::bind_rows(lapply(names(assignments)[-1],function(a) dplyr::bind_rows(lapply(names(assignments)[-1],function(b) data.frame(a=a,b=b,ARI=mclust::adjustedRandIndex(assignments[[a]],assignments[[b]]))))))
  write_table(ari,'pc_ARI_long')
  # Save pre-choice object; selection is documented after inspecting ARI/elbow.
  saveRDS(obj,checkpoint('03_embeddings'),compress=FALSE)
  session_save('03_embeddings')
}
