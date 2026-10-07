# Shared file handling only. Sourcing the scripts never starts an analysis.
# Policy: use package::function for non-base R functions; never use bare select.
core1_paths <- function(root = getwd(), output = "results/core1_reproduction") {
  root <- normalizePath(root, mustWork = TRUE)
  out <- file.path(root, output)
  for (sub in c("tables", "figures", "checkpoints", "diagnostics")) {
    dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
  }
  list(root = root, out = out)
}

core1_write <- function(x, p, name, sub = "tables", row.names = FALSE) {
  utils::write.csv(x, file.path(p$out, sub, name), row.names = row.names)
}

core1_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
}

core1_session <- function(p, stage) {
  writeLines(capture.output(utils::sessionInfo()),
             file.path(p$out, "diagnostics", paste0(stage, "_sessionInfo.txt")))
}

core1_samples <- function(root) {
  samples <- utils::read.delim(file.path(root, "data/rna/sample_metadata.tsv"),
                             stringsAsFactors = FALSE)
  stopifnot(nrow(samples) == 4L, !anyDuplicated(samples$sample),
            setequal(samples$sample, c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")))
  samples$stem <- c(GSM7504014 = "GSM7504014_Control1",
                    GSM7504015 = "GSM7504015_Control2",
                    GSM7504016 = "GSM7504016_Case1",
                    GSM7504017 = "GSM7504017_Case2")[samples$gsm]
  stopifnot(!anyNA(samples$stem))
  samples
}

core1_features <- function(root, samples) {
  lapply(seq_len(nrow(samples)), function(i) {
    path <- file.path(root, "data/rna", samples$sample[i],
                      paste0(samples$stem[i], "_features.tsv.gz"))
    utils::read.delim(gzfile(path), header = FALSE, stringsAsFactors = FALSE)[[1]]
  })
}

core1_plot <- function(plot, p, filename, width = 9, height = 6) {
  ggplot2::ggsave(file.path(p$out, "figures", filename), plot = plot,
                  width = width, height = height, dpi = 300)
}
