# P1 set enrichment is conditional on each population's testable background.
# Do not rank biological specificity by overlap percentage alone.
run_context <- function() {
  de<-utils::read.csv(file.path(out,'tables/pseudobulk_DE_all_branches.csv'))
  p1<-utils::read.delim(file.path(root,'data/project1/FULL_interaction_DEGs_padj0.05_log2FC1.tsv'),quote='',comment.char='')
  candidates<-unique(p1$TAIR)
  stopifnot(!anyNA(candidates),!anyDuplicated(de[c('branch','cell_population','gene')]))
  wt<-utils::read.delim(file.path(root,'data/core1_corrected_reference/P1_WT_drought_DEGs.tsv'),quote='',comment.char='')
  ia<-utils::read.delim(file.path(root,'data/core1_corrected_reference/P1_interaction_all_results.tsv'),quote='',comment.char='')
  controls<-setdiff(wt$TAIR,candidates)
  controls_strict<-intersect(controls,ia$TAIR[!is.na(ia$padj)&ia$padj>=.05])
  write_table(data.frame(gene=unique(c(candidates,controls)),P1_interaction=unique(c(candidates,controls))%in%candidates,
    WT_drought_non_candidate=unique(c(candidates,controls))%in%controls,
    WT_drought_interaction_nonsignificant=unique(c(candidates,controls))%in%controls_strict),'P1_control_gene_sets')
  ft <- function(d,set,background=NULL) {
    if(!is.null(background))d<-d[d$gene %in% c(set,background),]
    z<-d$gene %in% set; sig<-d$significant
    a<-sum(z&sig);bb<-sum(z&!sig);cc<-sum(!z&sig);dd<-sum(!z&!sig)
    m<-matrix(c(a,bb,cc,dd),nrow=2,byrow=TRUE,
      dimnames=list(c('set','background'),c('DEG','not_DEG')))
    if(any(rowSums(m)==0))return(data.frame(set_DEG=a,set_not_DEG=bb,bg_DEG=cc,bg_not_DEG=dd,odds_ratio=NA,CI_low=NA,CI_high=NA,pvalue=NA))
    f<-stats::fisher.test(m,alternative='greater'); ci<-stats::fisher.test(m)
    data.frame(set_DEG=a,set_not_DEG=bb,bg_DEG=cc,bg_not_DEG=dd,odds_ratio=unname(f$estimate),CI_low=ci$conf.int[1],CI_high=ci$conf.int[2],pvalue=f$p.value)
  }
  enrich<-list(); sums<-list()
  groups<-split(de,paste(de$branch,de$cell_population,sep='__'))
  for(key in names(groups)) {
    d<-groups[[key]]; b<-d$branch[1];p<-d$cell_population[1]
    z<-d$gene%in%candidates
    sums[[key]]<-data.frame(branch=b,cell_population=p,P1_fitted=sum(z),P1_finite_padj=sum(z&!is.na(d$padj)),P1_DEG=sum(z&d$significant),up=sum(z&d$significant&d$log2FoldChange>=1),down=sum(z&d$significant&d$log2FoldChange<=-1))
    for(u in c('finite_padj_primary','all_fitted_sensitivity')) {
      q<-if(u=='finite_padj_primary')d[!is.na(d$padj),] else d
      tests<-list(interaction_vs_background=ft(q,candidates),drought_only_vs_background=ft(q,controls),
        interaction_vs_drought_only=ft(q,candidates,controls),interaction_vs_drought_strict=ft(q,candidates,controls_strict))
      for(n in names(tests))enrich[[paste(key,u,n)]]<-cbind(data.frame(branch=b,cell_population=p,universe=u,test=n),tests[[n]])
    }
  }
  en<-dplyr::bind_rows(enrich) |> dplyr::group_by(branch,universe,test) |> dplyr::mutate(BH=stats::p.adjust(pvalue,method='BH')) |> dplyr::ungroup()
  sm<-dplyr::bind_rows(sums);sm$percent_fitted<-100*sm$P1_DEG/sm$P1_fitted
  write_table(sm,'P1_contextualization_summary');write_table(en,'P1_Fisher_enrichment_all')
  write_table(dplyr::filter(en,universe=='finite_padj_primary',test=='interaction_vs_background'),'P1_Fisher_primary')
  overlap<-de[de$gene%in%candidates,];sig<-overlap[overlap$significant,]
  write_table(overlap,'P1_candidate_all_results');write_table(sig,'P1_candidate_significant_results')
  write_table(data.frame(TAIR=candidates,in_common=candidates%in%readLines(file.path(out,'diagnostics/common_gene_features.txt')),
    in_any_fitted=candidates%in%overlap$gene,in_any_significant=candidates%in%sig$gene),'P1_candidate_input_status')
  go<-utils::read.delim(file.path(root,'data/project1/FULL_interaction_GO_enrichment_positive.tsv'),quote='',comment.char='')
  long<-tidyr::separate_rows(dplyr::select(go,ID,Description,geneID),geneID,sep='/')
  mp<-AnnotationDbi::select(org.At.tair.db::org.At.tair.db,keys=unique(long$geneID),keytype='SYMBOL',columns=c('TAIR','SYMBOL'))
  mp<-unique(mp);stopifnot(!anyNA(mp$TAIR),!anyDuplicated(mp$SYMBOL))
  gm<-dplyr::left_join(long,mp,by=c('geneID'='SYMBOL'));stopifnot(!anyNA(gm$TAIR))
  write_table(gm,'GO_gene_mapping')
  write_table(data.frame(TAIR=gm$TAIR,in_P1_candidates=gm$TAIR%in%candidates),'GO_candidate_membership_check')
  context<-dplyr::inner_join(gm,sig,by=c('TAIR'='gene'),relationship='many-to-many')
  write_table(context,'GO_gene_population_context')
  gs<-context |> dplyr::group_by(branch,cell_population,ID,Description) |> dplyr::summarise(DE_genes=dplyr::n_distinct(TAIR),up=dplyr::n_distinct(TAIR[log2FoldChange>=1]),down=dplyr::n_distinct(TAIR[log2FoldChange<=-1]),.groups='drop')
  combos<-dplyr::bind_rows(lapply(unique(de$branch),function(b) merge(data.frame(branch=b,cell_population=unique(de$cell_population[de$branch==b])),go[,c('ID','Description')],by=NULL)))
  complete<-dplyr::left_join(combos,gs,by=c('branch','cell_population','ID','Description'))
  for(v in c('DE_genes','up','down'))complete[[v]][is.na(complete[[v]])]<-0
  write_table(complete,'GO_population_complete')
  for(b in unique(de$branch))write_table(tidyr::pivot_wider(dplyr::select(complete[complete$branch==b,],ID,Description,cell_population,DE_genes),names_from=cell_population,values_from=DE_genes),paste0(b,'_GO_matrix'))
  sets<-split(gm$TAIR,gm$ID)
  pairs<-dplyr::bind_rows(lapply(names(sets),function(a)dplyr::bind_rows(lapply(names(sets),function(b) data.frame(term1=a,term2=b,n1=length(unique(sets[[a]])),n2=length(unique(sets[[b]])),overlap=length(intersect(sets[[a]],sets[[b]])),Jaccard=length(intersect(sets[[a]],sets[[b]]))/length(union(sets[[a]],sets[[b]])))))))
  write_table(pairs,'GO_gene_set_pairwise_overlap')
  f<-dplyr::filter(en,universe=='finite_padj_primary',test=='interaction_vs_background')
  plot_save(ggplot2::ggplot(f,ggplot2::aes(cell_population,-log10(BH),fill=branch))+ggplot2::geom_col(position='dodge')+ggplot2::geom_hline(yintercept=-log10(.05),linetype=2)+ggplot2::coord_flip()+ggplot2::theme_bw()+ggplot2::labs(y='-log10 BH-adjusted Fisher p',x=NULL),'06_P1_enrichment')
  plot_save(ggplot2::ggplot(complete,ggplot2::aes(cell_population,Description,fill=DE_genes))+ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=DE_genes),size=3)+ggplot2::facet_wrap(~branch,ncol=1,scales='free_x')+ggplot2::scale_fill_gradient(low='white',high='#416a9b')+ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=40,hjust=1))+ggplot2::labs(x=NULL,y=NULL),'06_GO_context',w=13,h=9)
  session_save('06_context')
}
