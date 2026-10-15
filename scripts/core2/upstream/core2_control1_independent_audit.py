# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
#!/usr/bin/env python3
"""Read-only Control 1 audit, independent of the batch QC and Signac counters.

Run: python3 /absolute/path/to/outputs/core2_control1_independent_audit.py
R is used only to read existing RDS/MatrixMarket objects and summarize/plot.
Python independently counts original BED fragment rows and every matrix entry.
No cell filtering, peak calling, mapping, or existing-output modifications.
"""
import array
import bisect
import csv
import gzip
import hashlib
import json
import math
import subprocess
import time
import os
from pathlib import Path

TASK = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
OUT = TASK / "outputs/core2_qc/control1_independent_audit"
WORK = TASK / "work/control1_independent_audit"
R = os.environ.get("P3_RSCRIPT", "Rscript")

PREPARE_R = r'''
args <- commandArgs(trailingOnly = TRUE)
task <- args[1]; work <- args[2]; out <- args[3]
.libPaths(c(file.path(task, "work", "R_library"), .libPaths()))
stopifnot(requireNamespace("Matrix", quietly = TRUE),
          requireNamespace("GenomicRanges", quietly = TRUE))
peaks <- readRDS(file.path(task, "outputs/core2_common_peaks.rds"))
upstream <- file.path(task, "outputs/control1_upstream")
mtx <- readRDS(file.path(upstream, "control1_common_peak_counts.rds"))
cells <- readLines(file.path(upstream, "matrix_barcodes.tsv"))
bed <- utils::read.table(file.path(task, "outputs/core2_common_peaks.bed"),
                         sep = "\t", colClasses = c("character", "integer", "integer"))
stopifnot(length(peaks) == 27677L, length(cells) == 3431L,
          GenomicRanges::isDisjoint(peaks), !anyDuplicated(cells),
          identical(as.vector(colnames(mtx)), cells),
          identical(as.vector(rownames(mtx)), paste(
            as.character(GenomicRanges::seqnames(peaks)),
            GenomicRanges::start(peaks), GenomicRanges::end(peaks), sep = "-")),
          identical(bed[[1]], as.character(GenomicRanges::seqnames(peaks))),
          all(bed[[2]] + 1L == GenomicRanges::start(peaks)),
          all(bed[[3]] == GenomicRanges::end(peaks)),
          all(mtx@x > 0), all(mtx@x == floor(mtx@x)))
Matrix::writeMM(mtx, file.path(work, "expected_control1.mtx"))
rm(mtx); invisible(gc())
qc <- utils::read.csv(file.path(task, "outputs/core2_qc/all_nuclei_qc.csv"),
                      check.names = FALSE)
sample_ids <- c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")
paired <- list()
for (s in sample_ids) {
  message("Read published ATAC-only counts: ", s)
  folder <- file.path(Sys.getenv("P3_ROOT", getwd()), "data", "rna", s)
  one <- function(pattern) {
    paths <- list.files(folder, pattern = pattern, full.names = TRUE)
    stopifnot(length(paths) == 1L); paths
  }
  f <- utils::read.table(one("features\\.tsv\\.gz$"), sep = "\t",
                         header = FALSE, quote = "", comment.char = "")
  b <- local({con <- gzfile(one("barcodes\\.tsv\\.gz$"), "rt")
              on.exit(close(con)); readLines(con)})
  x <- local({con <- gzfile(one("matrix\\.mtx\\.gz$"), "rt")
              on.exit(close(con)); Matrix::readMM(con)})
  peak_rows <- which(f[[3]] == "Peaks")
  stopifnot(nrow(x) == nrow(f), ncol(x) == length(b),
            length(peak_rows) > 0L, !anyDuplicated(b))
  d <- qc[qc$sample == s, , drop = FALSE]
  stopifnot(!anyDuplicated(d$barcode), setequal(d$barcode, b))
  d$published_peak_counts <- as.numeric(Matrix::colSums(
    x[peak_rows, , drop = FALSE]))[match(d$barcode, b)]
  d$published_peak_features <- length(peak_rows)
  paired[[s]] <- d
  rm(x, f); invisible(gc())
}
utils::write.csv(do.call(rbind, paired), file.path(out, "published_paired_counts.csv"),
                 row.names = FALSE)
utils::capture.output(utils::sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
'''

REPORT_R = r'''
args <- commandArgs(trailingOnly = TRUE); out <- args[1]
d <- utils::read.csv(file.path(out, "published_paired_counts.csv"))
summaries <- lapply(split(d, d$sample), function(x) data.frame(
  sample = x$sample[1], nuclei = nrow(x),
  median_published_counts = stats::median(x$published_peak_counts),
  median_rebuilt_counts = stats::median(x$nCount_ATAC),
  median_nuclear_fragments = stats::median(x$nuclear_fragments),
  median_organelle_pct = stats::median(x$percent_organelle),
  spearman_published_vs_rebuilt = stats::cor(x$published_peak_counts,
                                           x$nCount_ATAC, method = "spearman"),
  spearman_published_vs_nuclear = stats::cor(x$published_peak_counts,
                                           x$nuclear_fragments, method = "spearman")))
utils::write.csv(do.call(rbind, summaries), file.path(out, "published_concordance.csv"),
                 row.names = FALSE)
c1 <- d[d$sample == "control_rep1", ]
strata <- lapply(split(c1, ifelse(c1$nuclear_fragments < 500, "below_500", "at_least_500")),
  function(x) data.frame(group = ifelse(x$nuclear_fragments[1] < 500, "below_500", "at_least_500"),
    nuclei = nrow(x), median_published_counts = stats::median(x$published_peak_counts),
    median_rebuilt_counts = stats::median(x$nCount_ATAC),
    median_nuclear_fragments = stats::median(x$nuclear_fragments),
    median_organelle_pct = stats::median(x$percent_organelle),
    median_TSS = stats::median(x$TSS.enrichment)))
utils::write.csv(do.call(rbind, strata), file.path(out, "control1_descriptive_strata.csv"),
                 row.names = FALSE)
grDevices::png(file.path(out, "published_vs_rebuilt.png"), width = 1800, height = 1500, res = 160)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (s in c("control_rep1", "control_rep2", "stress_rep1", "stress_rep2")) {
  x <- d[d$sample == s, ]
  graphics::plot(log10(x$published_peak_counts + 1), log10(x$nCount_ATAC + 1),
    pch = 16, cex = .4, col = grDevices::adjustcolor("steelblue", alpha.f = .3),
    xlim = c(0, 5), ylim = c(0, 5), main = s,
    xlab = "log10(published ATAC counts + 1)", ylab = "log10(rebuilt common-peak counts + 1)")
}
invisible(grDevices::dev.off())
'''


def log(message):
    print(time.strftime("%H:%M:%S"), message, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def overlap_interval(starts, ends, start, end):
    # BED half-open: touching a boundary is not an overlap.
    return bisect.bisect_right(ends, start), bisect.bisect_left(starts, end)


def main():
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite prior audit: {OUT}")
    OUT.mkdir(parents=True)
    WORK.mkdir(parents=True, exist_ok=True)
    started = time.time()
    # Tests target independent boundary/PCR/nucleosome counting conventions.
    assert overlap_interval([100, 130], [110, 140], 110, 130) == (1, 1)
    assert overlap_interval([100, 130], [110, 140], 105, 135) == (0, 2)
    lengths = [146, 147, 148, 293, 294]
    assert sum(n < 147 for n in lengths) == 1
    assert sum(147 < n < 294 for n in lengths) == 2
    prepare = WORK / "read_inputs.R"
    prepare.write_text(PREPARE_R)
    log("Export existing matrix; read published ATAC rows only for four-sample comparison")
    subprocess.run([R, str(prepare), str(TASK), str(WORK), str(OUT)], check=True)

    upstream = TASK / "outputs/control1_upstream"
    original = upstream / "control1.fragments.tsv.gz"
    nuclear_copy = upstream / "control1.nuclear.fragments.tsv.gz"
    qc_path = TASK / "outputs/core2_qc/control_rep1/nucleus_qc.csv"
    preserved = {str(p): sha256(p) for p in (qc_path, original, nuclear_copy)}
    cells = (upstream / "matrix_barcodes.tsv").read_text().splitlines()
    cell_index = {b.encode(): i for i, b in enumerate(cells)}
    ncell = len(cells)
    peaks = {}
    with open(TASK / "outputs/core2_common_peaks.bed") as f:
        for row, line in enumerate(f):
            chrom, start, end = line.split()
            starts, ends, rows = peaks.setdefault(chrom.encode(), ([], [], []))
            start, end = int(start), int(end)
            assert start < end and (not ends or start >= ends[-1])
            starts.append(start); ends.append(end); rows.append(row)
    npeak = sum(len(p[0]) for p in peaks.values())
    assert npeak == 27677 and ncell == 3431 and len(cell_index) == ncell
    names = ("total_fragments", "nuclear_fragments", "mitochondrial_fragments",
             "plastid_fragments", "fragments_in_common_peaks", "nCount_ATAC",
             "nFeature_ATAC", "nucleosome_free", "mononucleosomal")
    counts = {name: [0] * ncell for name in names}
    # Dense integer array is ~380 MB; no Signac/GenomicRanges counting is invoked.
    matrix = array.array("I", [0]) * (npeak * ncell)
    assert matrix.itemsize >= 4
    subset_hash = hashlib.sha256()
    contigs = {}
    previous = None
    log(f"Independently scan original fragments; recount all {npeak * ncell:,} matrix positions")
    with gzip.open(original, "rb") as f:
        for scanned, line in enumerate(f, 1):
            fields = line.rstrip(b"\n").split(b"\t")
            if len(fields) != 5:
                raise ValueError(f"Invalid fragment row {scanned}")
            chrom, start, end, bc, support = fields
            start, end = int(start), int(end)
            assert end > start >= 0 and int(support) >= 1
            key = (chrom, start, end, bc)
            if previous is not None and key <= previous:
                raise ValueError(f"Duplicate or unsorted fragment at row {scanned}")
            previous = key
            assert chrom in peaks or chrom in (b"Mt", b"Pt")
            contigs[chrom.decode()] = contigs.get(chrom.decode(), 0) + 1
            is_nuclear = chrom in peaks
            if is_nuclear:
                subset_hash.update(line)
            cell = cell_index.get(bc)
            if cell is not None:
                # One deduplicated row = one fragment, regardless of PCR support.
                counts["total_fragments"][cell] += 1
                if is_nuclear:
                    counts["nuclear_fragments"][cell] += 1
                    length = end - start
                    counts["nucleosome_free"][cell] += length < 147
                    counts["mononucleosomal"][cell] += 147 < length < 294
                    starts, ends, rows = peaks[chrom]
                    left, right = overlap_interval(starts, ends, start, end)
                    if left < right:
                        counts["fragments_in_common_peaks"][cell] += 1
                        counts["nCount_ATAC"][cell] += right - left
                        for p in range(left, right):
                            at = rows[p] * ncell + cell
                            if matrix[at] == 0:
                                counts["nFeature_ATAC"][cell] += 1
                            matrix[at] += 1
                elif chrom == b"Mt":
                    counts["mitochondrial_fragments"][cell] += 1
                else:
                    counts["plastid_fragments"][cell] += 1
            if scanned % 5_000_000 == 0:
                log(f"Scanned {scanned:,} fragment rows")
    log("Compare complete matrix against existing RDS; verify nuclear copy")
    with open(WORK / "expected_control1.mtx") as f:
        header = next(f)
        assert "coordinate" in header and "general" in header
        for line in f:
            if not line.startswith("%"):
                nr, nc, expected_nnz = map(int, line.split()); break
        assert (nr, nc) == (npeak, ncell)
        mismatches = 0
        examples = []
        for seen, line in enumerate(f, 1):
            r, c, val = line.split()
            r, c, val = int(r) - 1, int(c) - 1, int(float(val))
            at = r * ncell + c
            if matrix[at] != val:
                mismatches += 1
                if len(examples) < 20:
                    examples.append({"row_1based": r + 1, "barcode": cells[c],
                                     "existing": val, "independent": matrix[at]})
            matrix[at] = 0
        assert seen == expected_nnz
    extra_nonzero = sum(v != 0 for v in matrix)
    del matrix
    nuclear_hash = hashlib.sha256()
    with gzip.open(nuclear_copy, "rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            nuclear_hash.update(chunk)
    with open(qc_path) as f:
        old = {r["barcode"]: r for r in csv.DictReader(f)}
    assert set(old) == set(cells)
    independent = []
    for i, bc in enumerate(cells):
        d = {"barcode": bc, **{k: v[i] for k, v in counts.items()}}
        d["FRiP_nuclear"] = d["fragments_in_common_peaks"] / d["nuclear_fragments"]
        d["percent_organelle"] = 100 * (d["mitochondrial_fragments"] + d["plastid_fragments"]) / d["total_fragments"]
        d["nucleosome_signal"] = d["mononucleosomal"] / d["nucleosome_free"] if d["nucleosome_free"] else math.nan
        independent.append(d)
    write_csv(OUT / "independent_control1_metrics.csv", independent)
    comparisons = []
    for metric in ("total_fragments", "nuclear_fragments", "fragments_in_common_peaks",
                   "nCount_ATAC", "nFeature_ATAC", "FRiP_nuclear", "percent_organelle",
                   "nucleosome_signal"):
        diffs = [abs(r[metric] - float(old[r["barcode"]][metric])) for r in independent]
        tol = 1e-8 if metric in ("FRiP_nuclear", "percent_organelle", "nucleosome_signal") else 0
        comparisons.append({"metric": metric, "nuclei_compared": ncell,
                            "nuclei_disagree": sum(not (v <= tol) for v in diffs),
                            "max_absolute_difference": max(diffs), "tolerance": tol})
    write_csv(OUT / "existing_vs_independent.csv", comparisons)
    validation = json.loads((upstream / "fragment_validation.json").read_text())
    checks = {
        "original_fragment_sha256_matches_upstream_manifest": preserved[str(original)] == validation["fragments_sha256"],
        "full_fragment_rows_match_upstream_manifest": scanned == validation["unique_fragments"],
        "contig_counts_match_upstream_manifest": contigs == validation["fragments_by_contig"],
        "nuclear_file_is_exact_original_nuclear_subset": nuclear_hash.hexdigest() == subset_hash.hexdigest(),
        "complete_count_matrix_matches_independent_recount": mismatches == 0 and extra_nonzero == 0,
        "all_eight_qc_metrics_match": all(r["nuclei_disagree"] == 0 for r in comparisons),
        "original_inputs_unmodified": all(sha256(p) == h for p, h in preserved.items()),
    }
    facts = {"checks": checks, "fragment_rows_scanned": scanned,
             "matrix_dimensions": [npeak, ncell], "matrix_positions": npeak * ncell,
             "matrix_nonzero_entries_compared": expected_nnz,
             "matrix_mismatched_existing_entries": mismatches,
             "matrix_unexpected_nonzero_entries": extra_nonzero, "mismatch_examples": examples,
             "contig_counts": contigs, "input_sha256": preserved,
             "nucleosome_definition": "free length <147; mono 147<length<294; 147/294 excluded (Signac source)",
             "no_filtering_applied": True, "elapsed_seconds": round(time.time() - started, 2)}
    (OUT / "audit_checks.json").write_text(json.dumps(facts, indent=2) + "\n")
    report_r = WORK / "summarize.R"
    report_r.write_text(REPORT_R)
    subprocess.run([R, str(report_r), str(OUT)], check=True)
    log(json.dumps(facts, indent=2))
    if not all(checks.values()):
        raise RuntimeError("At least one audit failed; inspect audit_checks.json before interpreting biology")
    log(f"AUDIT COMPLETE: {OUT}")


if __name__ == "__main__":
    main()
