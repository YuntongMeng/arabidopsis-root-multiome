# Project 3 Core 1: contextualize P1 BRL3 x drought candidates in P3 WT
# root-tip sorbitol/osmotic stress. Conditions/genotypes differ: contextualization
# is the intended interpretation; neither candidate nor GO overlap is validation.
# Final input 137 candidates; 129 occur in >=1 final population DE result table.
# Pericycle: 23/118 = 19.49% (6 up, 17 down), the largest fraction.
# GO: 7 annotations, 36 term-symbol rows representing 13 distinct symbols.
# All 36 rows mapped; Pericycle overlaps all 7 terms; cold response=3, all down.
# Source 00_core1_common.R first. This function only summarizes saved DE results;
# it never performs GO enrichment, expression normalization, or DE testing.
run_core1_04 <- function(root = getwd(), output = "results/core1_reproduction",
                         de_file = "results/core1_frozen/de_all.csv",
                         use_frozen_mapping = TRUE) {
  core1_require(c("dplyr", "tidyr"))
  p <- core1_paths(root, output)
  de_all <- utils::read.csv(file.path(p$root, de_file), stringsAsFactors = FALSE)
  p1_deg <- utils::read.delim(file.path(p$root,
    "data/project1/FULL_interaction_DEGs_padj0.05_log2FC1.tsv"), stringsAsFactors = FALSE)
  p1_genes <- unique(p1_deg$TAIR)
  stopifnot(length(p1_genes) == 137L, !anyNA(p1_genes),
            !anyDuplicated(de_all[c("gene", "cell_population")]))
  eligible <- readLines(file.path(p$root, "results/core1_frozen/eligible_populations.txt"))
  stopifnot(setequal(unique(de_all$cell_population), eligible))
  # Eligibility/prefilter determines each denominator: it is not always 137.
  # P1_tested counts candidate rows retained in the DESeq2 results, including
  # rows with NA padj. It is not the count passing independent filtering.
  p1_p3_overlap <- dplyr::filter(de_all, gene %in% p1_genes)
  p1_p3_sig <- p1_p3_overlap |>
    dplyr::filter(!is.na(padj), padj < .05, abs(log2FoldChange) >= 1)
  candidate_summary <- p1_p3_overlap |>
    dplyr::group_by(cell_population) |>
    dplyr::summarise(P1_tested = dplyr::n_distinct(gene),
      P1_stress_DEG = dplyr::n_distinct(gene[!is.na(padj) & padj < .05 & abs(log2FoldChange) >= 1]),
      up = dplyr::n_distinct(gene[!is.na(padj) & padj < .05 & log2FoldChange >= 1]),
      down = dplyr::n_distinct(gene[!is.na(padj) & padj < .05 & log2FoldChange <= -1]), .groups = "drop") |>
    tidyr::complete(cell_population = eligible,
                    fill = list(P1_tested = 0L, P1_stress_DEG = 0L, up = 0L, down = 0L)) |>
    dplyr::mutate(percent_stress_DEG = ifelse(P1_tested > 0, 100 * P1_stress_DEG / P1_tested, NA_real_)) |>
    dplyr::arrange(dplyr::desc(percent_stress_DEG))
  common <- readLines(file.path(p$root, "results/core1_frozen/common_features.txt"))
  input_status <- data.frame(TAIR = p1_genes, in_common_features = p1_genes %in% common,
    in_any_population_result = p1_genes %in% p1_p3_overlap$gene,
    significant_in_any_population = p1_genes %in% p1_p3_sig$gene)
  core1_write(input_status, p, "P1_candidate_input_status.csv")
  core1_write(p1_p3_overlap, p, "P1_P3_all_candidate_results.csv")
  core1_write(p1_p3_sig, p, "P1_P3_significant_candidate_results.csv")
  core1_write(candidate_summary, p, "P1_P3_contextualization_summary.csv")

  # Read the existing seven enriched GO terms; do NOT rerun enrichment in P3.
  p1_go <- utils::read.delim(file.path(p$root,
    "data/project1/FULL_interaction_GO_enrichment_positive.tsv"), stringsAsFactors = FALSE)
  stopifnot(nrow(p1_go) == 7L, !anyDuplicated(p1_go$ID))
  term_symbols <- p1_go |>
    dplyr::select(ID, Description, geneID) |>
    tidyr::separate_rows(geneID, sep = "/") |>
    dplyr::rename(SYMBOL = geneID)
  stopifnot(nrow(term_symbols) == 36L)

  # Initial reverse lookup against P1's own SYMBOL column left 23/36 rows NA.
  # NA meant symbol-string mismatch (alias/case/old names), not gene absence.
  # AnnotationDbi::select(org.At.tair.db, keytype="SYMBOL") repaired all 36 rows.
  # Preserve that resolved mapping as a versioned input for stable reproduction.
  # 36/36 refers to term-gene memberships, NOT 36 distinct genes (there are 13).
  if (use_frozen_mapping) {
    symbol_map <- utils::read.csv(file.path(p$root, "data/core1_reference/symbol_map_db.csv"),
                                  stringsAsFactors = FALSE)
  } else {
    core1_require(c("AnnotationDbi", "org.At.tair.db"))
    symbol_map <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db,
      keys = unique(term_symbols$SYMBOL), keytype = "SYMBOL", columns = c("TAIR", "SYMBOL"))
  }
  symbol_map <- dplyr::distinct(symbol_map, SYMBOL, TAIR)
  if (anyNA(symbol_map$TAIR) || anyDuplicated(symbol_map$SYMBOL)) {
    stop("Unresolved or ambiguous GO symbol mapping. Resolve explicitly before joining.")
  }
  go_map <- dplyr::left_join(term_symbols, symbol_map, by = "SYMBOL", relationship = "many-to-one")
  stopifnot(nrow(go_map) == 36L, !anyNA(go_map$TAIR))
  # Expected many-to-many join: a gene can belong to several GO terms and can
  # be stress-responsive in several populations. Count distinct TAIR per term.
  go_context <- dplyr::inner_join(go_map, p1_p3_sig, by = c("TAIR" = "gene"),
                                  relationship = "many-to-many")
  go_summary <- go_context |>
    dplyr::group_by(ID, Description, cell_population) |>
    dplyr::summarise(stress_DE_genes = dplyr::n_distinct(TAIR),
      up = dplyr::n_distinct(TAIR[log2FoldChange >= 1]),
      down = dplyr::n_distinct(TAIR[log2FoldChange <= -1]), .groups = "drop") |>
    dplyr::arrange(Description, dplyr::desc(stress_DE_genes))
  # The first pivot contained only eight populations with nonzero overlap:
  # seven rows x nine columns including Description. Zero-overlap populations
  # had vanished entirely. Complete all terms x all 13 eligible populations.
  complete_go <- go_summary |>
    tidyr::complete(tidyr::nesting(ID, Description), cell_population = eligible,
                    fill = list(stress_DE_genes = 0L, up = 0L, down = 0L))
  go_matrix <- complete_go |>
    dplyr::select(Description, cell_population, stress_DE_genes) |>
    tidyr::pivot_wider(names_from = cell_population, values_from = stress_DE_genes, values_fill = 0)
  stopifnot(nrow(go_matrix) == 7L, ncol(go_matrix) == 14L)
  core1_write(go_map, p, "P1_GO_symbol_TAIR_mapping.csv")
  core1_write(go_context, p, "P1_GO_P3_gene_context.csv")
  core1_write(go_summary, p, "P1_GO_population_nonzero_summary.csv")
  core1_write(complete_go, p, "P1_GO_population_complete_summary.csv")
  core1_write(go_matrix, p, "P1_GO_population_matrix_complete.csv")
  # JA metabolism and long-chain fatty-acid metabolism share the same four
  # symbols. Water/water-deprivation share six; acid response also shares those
  # six. These are overlapping annotations, not seven independent programs.
  core1_session(p, "04")
  invisible(list(candidate_summary = candidate_summary, go_summary = go_summary, go_matrix = go_matrix))
}
