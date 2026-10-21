#!/usr/bin/env python3
"""Validate the recorded Core 2 publication package; no biological refitting.

Run from the project root, or set P3_ROOT. Writes nothing by default.
Reported counts are version-specific record invariants, not required outcomes
for a new analysis. A changed run needs a new reviewed result record.
"""
import csv
import gzip
import hashlib
import math
import os
import json
from pathlib import Path
from collections import Counter


def read_csv(path):
    path = Path(path)
    if not path.exists():
        path = Path(str(path) + ".gz")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def true(value):
    return str(value).upper() == "TRUE"


def verify(root):
    root = Path(root)
    tables = root / "results/core2/tables"
    retention = read_csv(tables / "retention_summary.csv")
    samples = ["control_rep1", "control_rep2", "stress_rep1", "stress_rep2"]
    retained = dict(zip(samples, [1239, 3817, 4841, 4249]))
    assert [d["sample"] for d in retention] == samples
    for row in retention:
        assert int(row["basic_ATAC_QC"]) == sum(int(row[k]) for k in
            ["RNA_singlet_retained", "RNA_doublet_excluded", "RNA_not_tested_excluded"])
        assert int(row["RNA_singlet_retained"]) == retained[row["sample"]]
    md = read_csv(tables / "nucleus_annotations.csv.gz")
    assert len(md) == 14146 and len({r["cell"] for r in md}) == len(md)
    assert dict(Counter(r["sample"] for r in md)) == retained
    assert all(r["cell"] == r["sample"] + "_" + r["barcode"] for r in md)
    assert all(r["condition"] == ("control" if r["sample"].startswith("control") else "stress") for r in md)
    assert all(float(r["nuclear_fragments"]) >= 1000 and float(r["TSS.enrichment"]) > 1 and
               float(r["nucleosome_signal"]) < 2 for r in md)

    groups = Counter((r["sample"], r["rna_population"]) for r in md)
    pb = read_csv(tables / "atac_pseudobulk_metadata.csv")
    assert len(pb) == len(groups) == 46
    assert all(groups[(r["sample"], r["rna_population"])] == int(r["n_nuclei"]) for r in pb)
    assert not any(r["condition"] == "stress" and r["rna_population"] == "Unresolved_control_enriched" for r in pb)

    previous = {}
    intervals = 0
    with (root / "data/core2_reference/common_peaks.bed").open() as handle:
        for line in handle:
            chrom, start, end = line.rstrip().split("\t")
            start, end = int(start), int(end)
            assert chrom in list("12345") and start >= 0 and 20 < end - start < 10000
            assert start >= previous.get(chrom, 0)
            previous[chrom] = end
            intervals += 1
    assert intervals == 27677

    summary = read_csv(tables / "differential_accessibility/differential_accessibility_summary.csv")
    total_significant = 0
    for row in summary:
        population = row["population"]
        full = read_csv(tables / "differential_accessibility" / population / "all_peak_results.csv")
        assert len(full) == int(row["peaks_in_model"])
        assert len({r["peak"] for r in full}) == len(full)
        significant = [r for r in full if r["padj"] not in {"", "NA"} and float(r["padj"]) < .05]
        assert len(significant) == int(row["FDR_significant"])
        assert all(true(r["FDR_significant"]) == (r["padj"] not in {"", "NA"} and float(r["padj"]) < .05) for r in full)
        assert sum(float(r["log2FoldChange"]) >= 1 for r in significant) == int(row["higher_in_stress_2fold"])
        assert sum(float(r["log2FoldChange"]) <= -1 for r in significant) == int(row["lower_in_stress_2fold"])
        total_significant += len(significant)
        if population == "Atrichoblast_lineage":
            assert len(significant) == 1
            candidate = significant[0]
            assert candidate["peak"] == "4-5457000-5457860"
            assert math.isclose(float(candidate["padj"]), .00452605090866568, rel_tol=1e-10)
            assert min(float(candidate["normalized_"+s]) for s in samples[2:]) > max(float(candidate["normalized_"+s]) for s in samples[:2])
    assert total_significant == 1

    p1 = read_csv(tables / "P1_contextualization/P1_gene_status_all_populations.csv")
    p1_summary = read_csv(tables / "P1_contextualization/P1_summary.csv")
    assert len(p1) == 685 and len({(r["population"], r["TAIR"]) for r in p1}) == 685
    assert all(true(r["in_reference"]) and true(r["in_nuclear_analysis"]) for r in p1)
    for row in p1_summary:
        d = [r for r in p1 if r["population"] == row["population"]]
        assert len(d) == int(row["P1_total"]) == 137
        for field, detail in [("with_common_peak", "n_common_peaks"), ("modelled", "n_modelled_peaks"),
                              ("with_FDR_peak", "n_FDR_peaks"), ("with_priority_peak", "n_priority_peaks")]:
            assert sum(int(r[detail]) > 0 for r in d) == int(row[field])
        assert int(row["with_common_peak"]) == 120
        assert int(row["with_FDR_peak"]) == int(row["with_priority_peak"]) == 0

    # Check the manual rerun and the read-only audit as separate evidence.
    diagnostics = root / "results/core2/diagnostics"
    rerun = read_csv(diagnostics / "manual_rerun/validation_comparison.csv")
    assert len(rerun) == 10 and all(true(r["matches_recorded_run"]) for r in rerun)
    record = json.loads((diagnostics / "manual_rerun/rerun_record.json").read_text())
    entry = diagnostics / "manual_rerun" / record["executed_entry_point_archive"]
    assert hashlib.sha256(entry.read_bytes()).hexdigest() == record["entry_point_sha256"]
    for name, digest in record["stage_source_sha256"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest, name
    audit = read_csv(diagnostics / "scientific_audit/checks.csv")
    assert len(audit) == 111 and all(true(r["passed"]) for r in audit)
    fragment_map = read_csv(diagnostics / "scientific_audit/fragment_mapping_checks.csv")
    assert len(fragment_map) == 4
    assert {r["sample"]: int(r["selected_nuclei"]) for r in fragment_map} == retained
    assert all(true(r["raw_barcode_mapping_matches"]) and true(r["nuclear_path_matches"]) for r in fragment_map)
    annotation = read_csv(diagnostics / "scientific_audit/annotation_checks.csv")
    assert len(annotation) == 2 and all(true(r["passed"]) for r in annotation)

    manifest = diagnostics / "publication_manifest.csv"
    if manifest.exists():
        for row in read_csv(manifest):
            path = root / row["path"]
            assert path.is_file(), path
            assert path.stat().st_size == int(row["bytes"]), path
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], path
    return {"retained_nuclei": len(md), "common_peaks": intervals,
            "observed_pseudobulks": len(pb), "FDR_significant_peaks": total_significant,
            "P1_candidates": 137, "P1_with_common_peak": 120, "P1_with_FDR_peak": 0,
            "manual_rerun_comparisons_passed": len(rerun), "scientific_audit_checks_passed": len(audit)}


if __name__ == "__main__":
    root = Path(os.environ.get("P3_ROOT", Path(__file__).resolve().parent.parent))
    result = verify(root)
    print("Core 2 recorded-result verification: PASS")
    for key, value in result.items():
        print(f"  {key}: {value}")
