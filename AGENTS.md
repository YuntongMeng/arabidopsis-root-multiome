# Project 3 editing conventions

- Preserve the decision trail and distinguish formal results from diagnostics.
- Use explicit R package namespaces for non-base functions, especially ambiguous
  names: dplyr::select, dplyr::filter, dplyr::rename, tidyr::pivot_longer,
  AnnotationDbi::mapIds, DESeq2::results, SeuratObject::JoinLayers.
- Never use a bare select() or Select().
- Preserve results/core1_frozen as the archived completed analysis.
- Do not rerun biological analysis or push GitHub without the user's instruction.
- Keep uncertain annotation provenance visible; do not fabricate citations.
