# Core 2: describe ATAC peak support around the full Project 1 candidate list.
# source() defines the function; the function call performs the comparison.
# This is positional contextualization, not target-gene validation, gene-level
# differential testing, GO enrichment, or a new multiple-testing procedure.

core2_P1_contextualization <- function(
    task_dir,
    project_dir = Sys.getenv("P3_ROOT", getwd()),
    analysis_dir = file.path(task_dir, "outputs", "core2_analysis"),
    populations = c(
      "Atrichoblast_lineage", "Cortex_lineage", "Endodermis_lineage",
      "Pericycle", "Root_cap"
    ),
    tss_upstream = 2000L,
    tss_downstream = 500L) {

  for (package in c("rtracklayer", "GenomicRanges", "IRanges", "S4Vectors")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Missing R package: ", package)
    }
  }
  stopifnot(
    length(tss_upstream) == 1L, is.finite(tss_upstream), tss_upstream >= 0,
    length(tss_downstream) == 1L, is.finite(tss_downstream), tss_downstream >= 1,
    length(populations) > 0L, !anyDuplicated(populations)
  )
  p1_path <- file.path(
    project_dir, "data", "project1",
    "FULL_interaction_DEGs_padj0.05_log2FC1.tsv"
  )
  common_path <- file.path(task_dir, "outputs", "core2_common_peaks.rds")
  gtf_path <- file.path(
    task_dir, "work", "annotation", "Arabidopsis_thaliana.TAIR10.62.gtf.gz"
  )
  stopifnot(file.exists(p1_path), file.exists(common_path), file.exists(gtf_path))
  results_dir <- file.path(analysis_dir, "differential_accessibility")
  out_dir <- file.path(analysis_dir, "P1_contextualization")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p1 <- utils::read.delim(
    p1_path, quote = "", comment.char = "", check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(all(c("TAIR", "SYMBOL", "log2FoldChange", "padj") %in% names(p1)))
  p1$TAIR <- trimws(as.character(p1$TAIR))
  stopifnot(!anyNA(p1$TAIR), all(nzchar(p1$TAIR)), !anyDuplicated(p1$TAIR))

  common <- readRDS(common_path)
  stopifnot(inherits(common, "GRanges"))
  common_ids <- paste0(
    as.character(GenomicRanges::seqnames(common)), "-",
    GenomicRanges::start(common), "-", GenomicRanges::end(common)
  )
  stopifnot(
    !anyDuplicated(common_ids),
    all(as.character(GenomicRanges::seqnames(common)) %in% as.character(1:5))
  )
  full_annotation <- rtracklayer::import(gtf_path, format = "gtf")
  genes <- full_annotation[which(full_annotation$type == "gene")]
  stopifnot(!anyDuplicated(as.character(genes$gene_id)))
  gene_index <- base::match(p1$TAIR, as.character(genes$gene_id))
  known <- !is.na(gene_index)
  reference <- data.frame(
    TAIR = p1$TAIR,
    P1_symbol = p1$SYMBOL,
    P1_interaction_log2FC = p1$log2FoldChange,
    P1_interaction_padj = p1$padj,
    in_reference = known,
    gene_name = NA_character_, chromosome = NA_character_,
    gene_start = NA_integer_, gene_end = NA_integer_,
    strand = NA_character_, TSS_position = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (any(known)) {
    gg <- genes[gene_index[known]]
    reference$gene_name[known] <- as.character(gg$gene_name)
    reference$chromosome[known] <- as.character(GenomicRanges::seqnames(gg))
    reference$gene_start[known] <- GenomicRanges::start(gg)
    reference$gene_end[known] <- GenomicRanges::end(gg)
    reference$strand[known] <- as.character(GenomicRanges::strand(gg))
    reference$TSS_position[known] <- ifelse(
      reference$strand[known] == "+",
      reference$gene_start[known], reference$gene_end[known]
    )
  }
  unnamed <- is.na(reference$gene_name) | reference$gene_name == ""
  reference$gene_name[unnamed] <- reference$TAIR[unnamed]
  reference$in_nuclear_analysis <- known &
    reference$chromosome %in% as.character(1:5)
  nuclear_candidate_index <- which(reference$in_nuclear_analysis)
  candidate_genes <- genes[gene_index[nuclear_candidate_index]]
  stopifnot(all(as.character(GenomicRanges::strand(candidate_genes)) %in% c("+", "-")))

  # GTF import and common GRanges both use 1-based, closed coordinates.
  tss_windows <- GenomicRanges::promoters(
    candidate_genes, upstream = tss_upstream, downstream = tss_downstream
  )
  if (length(tss_windows)) {
    GenomicRanges::start(tss_windows) <- pmax(1L, GenomicRanges::start(tss_windows))
  }
  body_hits <- IRanges::findOverlaps(common, candidate_genes, ignore.strand = TRUE)
  tss_hits <- IRanges::findOverlaps(common, tss_windows, ignore.strand = TRUE)
  body_pairs <- data.frame(
    peak_index = S4Vectors::queryHits(body_hits),
    candidate_index = S4Vectors::subjectHits(body_hits)
  )
  tss_pairs <- data.frame(
    peak_index = S4Vectors::queryHits(tss_hits),
    candidate_index = S4Vectors::subjectHits(tss_hits)
  )
  pairs <- unique(rbind(body_pairs, tss_pairs))
  pair_keys <- function(x) paste(x$peak_index, x$candidate_index, sep = "__")
  pair_id <- pair_keys(pairs)
  gene_rows <- nuclear_candidate_index[pairs$candidate_index]
  mapping <- data.frame(
    TAIR = reference$TAIR[gene_rows],
    gene_name = reference$gene_name[gene_rows],
    peak = common_ids[pairs$peak_index],
    chromosome = as.character(GenomicRanges::seqnames(common))[pairs$peak_index],
    peak_start = GenomicRanges::start(common)[pairs$peak_index],
    peak_end = GenomicRanges::end(common)[pairs$peak_index],
    overlaps_gene_body = pair_id %in% pair_keys(body_pairs),
    overlaps_TSS_window = pair_id %in% pair_keys(tss_pairs),
    stringsAsFactors = FALSE
  )
  mapping$TSS_distance_bp <- pmax(
    mapping$peak_start - reference$TSS_position[gene_rows],
    reference$TSS_position[gene_rows] - mapping$peak_end, 0
  )
  mapping$n_P1_candidates_for_peak <- as.integer(
    table(mapping$peak)[mapping$peak]
  )
  stopifnot(!anyDuplicated(mapping[c("TAIR", "peak")]))
  utils::write.csv(mapping, file.path(out_dir, "P1_gene_peak_matches.csv"), row.names = FALSE)
  utils::write.csv(reference, file.path(out_dir, "P1_reference_status.csv"), row.names = FALSE)

  summaries <- list()
  all_status <- list()
  all_long <- list()
  sample_order <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
  requested_fields <- c(
    "baseMean", "log2FoldChange", "pvalue", "padj",
    "FDR_significant", "priority_change", "category", "replicate_pattern",
    "direction_supported", paste0("raw_", sample_order),
    paste0("normalized_", sample_order)
  )

  for (population in populations) {
    result_path <- file.path(results_dir, population, "all_peak_results.csv")
    stopifnot(file.exists(result_path))
    r <- utils::read.csv(result_path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(
      all(c("peak", "population", "analysis_status", requested_fields) %in% names(r)),
      !anyDuplicated(r$peak), all(r$peak %in% common_ids),
      all(r$population == population), length(unique(r$analysis_status)) == 1L
    )
    ri <- base::match(mapping$peak, r$peak)
    long <- cbind(mapping, r[ri, requested_fields, drop = FALSE])
    rownames(long) <- NULL
    long$in_population_model <- !is.na(ri)
    long$population <- rep(population, nrow(long))
    long$analysis_status <- rep(unique(r$analysis_status), nrow(long))
    all_long[[population]] <- long

    status <- reference
    status$population <- population
    status$analysis_status <- unique(r$analysis_status)
    for (column in c(
      "n_common_peaks", "n_modelled_peaks", "n_peaks_with_adjusted_p",
      "n_FDR_peaks", "n_priority_peaks", "n_direction_supported_priority_peaks"
    )) status[[column]] <- 0L
    status$status <- "Reference_annotation_missing"
    grouped <- split(seq_len(nrow(long)), factor(long$TAIR, levels = p1$TAIR), drop = FALSE)

    for (i in seq_len(nrow(status))) {
      d <- long[grouped[[p1$TAIR[i]]], , drop = FALSE]
      status$n_common_peaks[i] <- nrow(d)
      status$n_modelled_peaks[i] <- sum(d$in_population_model)
      status$n_peaks_with_adjusted_p[i] <- sum(d$in_population_model & !is.na(d$padj))
      status$n_FDR_peaks[i] <- sum(d$FDR_significant %in% TRUE)
      status$n_priority_peaks[i] <- sum(d$priority_change %in% TRUE)
      status$n_direction_supported_priority_peaks[i] <- sum(d$direction_supported %in% TRUE)
      status$status[i] <- if (!status$in_reference[i]) {
        "Reference_annotation_missing"
      } else if (!status$in_nuclear_analysis[i]) {
        "Outside_nuclear_analysis"
      } else if (status$n_common_peaks[i] == 0L) {
        "No_common_peak_match"
      } else if (status$n_modelled_peaks[i] == 0L) {
        "Matching_peaks_not_modelled_in_population"
      } else if (status$n_peaks_with_adjusted_p[i] == 0L) {
        "Modelled_peaks_without_adjusted_p"
      } else if (status$n_priority_peaks[i] > 0L) {
        "Nearby_priority_peak_present"
      } else if (status$n_FDR_peaks[i] > 0L) {
        "Nearby_FDR_peak_below_priority_effect_cutoff"
      } else {
        "Modelled_peaks_without_FDR_support"
      }
    }
    stopifnot(nrow(status) == nrow(p1), !anyNA(status$status))
    all_status[[population]] <- status
    summaries[[population]] <- data.frame(
      population = population,
      P1_total = nrow(p1),
      with_common_peak = sum(status$n_common_peaks > 0L),
      modelled = sum(status$n_modelled_peaks > 0L),
      with_FDR_peak = sum(status$n_FDR_peaks > 0L),
      with_priority_peak = sum(status$n_priority_peaks > 0L),
      row.names = NULL
    )
  }

  summary <- do.call(rbind, summaries)
  status_table <- do.call(rbind, all_status)
  long_table <- do.call(rbind, all_long)
  rownames(summary) <- rownames(status_table) <- rownames(long_table) <- NULL
  stopifnot(nrow(status_table) == nrow(p1) * length(populations))
  utils::write.csv(summary, file.path(out_dir, "P1_summary.csv"), row.names = FALSE)
  utils::write.csv(status_table, file.path(out_dir, "P1_gene_status_all_populations.csv"),
                   row.names = FALSE)
  utils::write.csv(long_table, file.path(out_dir, "P1_peak_results_long.csv"), row.names = FALSE)

  saveRDS(list(
    candidate_source = normalizePath(p1_path),
    common_peak_source = normalizePath(common_path),
    annotation_source = normalizePath(gtf_path),
    input_results_dir = normalizePath(results_dir), populations = populations,
    rule = "Overlap gene body OR strand-aware TSS window; both flags retained",
    tss_upstream = tss_upstream, tss_downstream = tss_downstream,
    coordinate_system = "1-based closed GRanges; no new BED conversion",
    allow_multiple_candidate_genes_per_peak = TRUE,
    new_gene_level_significance_test = FALSE
  ), file.path(out_dir, "parameters.rds"))
  writeLines(c(
    "# Project 1 candidate contextualization in Core 2", "",
    "This describes genomic proximity to the original Project 1 interaction candidates.",
    paste0("Mapping rule: overlap the gene body OR the strand-aware TSS window (upstream ",
           tss_upstream, " bp; downstream ", tss_downstream, " bp)."),
    "The two overlap flags are retained. A peak may match multiple candidate genes.",
    "This proximity window is separate from the earlier GeneActivity 500-bp upstream window.",
    "Peaks, counts, models, p values, and primary selection thresholds are unchanged.",
    "Counts in P1_summary refer to candidate genes with at least one matching peak.",
    "with_FDR_peak does not mean the gene itself has a gene-level adjusted p value.",
    "Multiple nearby peaks are not treated as independent gene-level tests.",
    "P1 interaction log2FC is not the same contrast as Core 2 stress/control log2FC.",
    "Matching signs cannot by themselves establish cross-project agreement.",
    "No_common_peak_match does not prove a gene is inactive or chromatin is closed.",
    "Modelled_peaks_without_FDR_support does not prove absence of a stress response.",
    "Reference_annotation_missing and Outside_nuclear_analysis are retained explicitly.",
    "Root_cap remains exploratory because of weak within-condition reproducibility.",
    "No GO analysis, target-gene validation, or new enrichment/significance test is performed."
  ), file.path(out_dir, "README.md"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  print(summary, row.names = FALSE)
  invisible(summary)
}
