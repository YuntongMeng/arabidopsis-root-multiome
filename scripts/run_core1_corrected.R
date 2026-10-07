#!/usr/bin/env Rscript
# Run from the Project 3 root, or set P3_ROOT. New dependencies may be supplied
# through optional P3_R_LIB; otherwise standard R libraries are used.
args <- commandArgs(trailingOnly=TRUE)
stage <- if(length(args))args[1] else 'help'
root_arg <- Sys.getenv('P3_ROOT',getwd())
source(file.path(root_arg,'scripts/core1/00_common.R'))
for(f in list.files(file.path(root,'scripts/core1'),pattern='^[0-9][0-9]_.*\\.R$',full.names=TRUE)) if(!grepl('00_common',f))source(f)
source(file.path(root,'exploration/core1/08_extra_diagnostics.R'))
source(file.path(root,'exploration/core1/09_render_figures.R'))
if(stage=='help') {
  cat('Stages: audit, qc, embeddings, evidence, de, context, compare, diagnostics, figures, all\n')
  cat('Historical diagnostics are optional and excluded from all.\n')
  cat('Run all explicitly to recompute. Annotation is bound to the reviewed partition MD5.\n')
  cat('If partition changes, review marker evidence and mapping before de; never bypass the check.\n')
} else {
  calls<-list(audit=run_input_audit,qc=run_qc_doublets,embeddings=run_embeddings,
    evidence=function()run_branch_evidence(npcs=25),de=run_annotation_pseudobulk,
    context=run_context,compare=run_compare_validate,diagnostics=run_extra_diagnostics,figures=run_render_figures)
  stopifnot(stage%in%c(names(calls),'all'))
  if(stage=='all')for(fun in calls[setdiff(names(calls),'diagnostics')])fun() else calls[[stage]]()
}
