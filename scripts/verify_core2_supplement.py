#!/usr/bin/env python3
"""Verify supplemental results without fitting models or changing files."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


def read(path):
    path = Path(path)
    if not path.exists():
        path = Path(str(path) + '.gz')
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def truth(value):
    return str(value).upper() == 'TRUE'


def verify(root):
    root = Path(root)
    supplement = root / 'results/core2/supplement'
    diag = read(supplement / 'model_diagnostics/diagnostic_summary.csv')
    assert len(diag) == 5
    for row in diag:
        n = int(row['finite_initial_estimates'])
        floor = int(row['initial_at_boundary'])
        assert n == int(row['model_peaks']) and 0 <= floor <= n
        assert math.isclose(float(row['percent_initial_at_boundary']), 100*floor/n, abs_tol=1e-10)
        assert int(row['raw_p_below_001']) <= int(row['raw_p_below_005'])
    assert len(list((supplement / 'model_diagnostics').glob('*.png'))) == 10
    assert len(list((supplement / 'model_diagnostics').glob('*.pdf'))) == 10
    comparisons = read(supplement / 'sensitivity_count5/threshold_comparison.csv')
    assert len(comparisons) == 10
    samples = ['control_rep1', 'control_rep2', 'stress_rep1', 'stress_rep2']
    for setting in comparisons:
        pop = setting['population']
        if setting['setting'] == 'primary_count10':
            results = read(root / 'results/core2/tables/differential_accessibility' / pop / 'all_peak_results.csv')
        else:
            assert setting['setting'] == 'exploratory_count5'
            results = read(supplement / 'sensitivity_count5' / pop / 'all_peak_results.csv')
        assert len(results) == int(setting['peaks_in_model'])
        assert len({r['peak'] for r in results}) == len(results)
        significant = [r for r in results if truth(r['FDR_significant'])]
        assert len(significant) == int(setting['FDR_significant'])
        assert all(sum(int(r['raw_'+s]) >= int(setting['min_count']) for s in samples) >= 2 for r in results)
        assert all(truth(r['FDR_significant']) == (r['padj'] not in {'NA',''} and float(r['padj']) < .05) for r in results)
    shared = read(supplement / 'sensitivity_count5/shared_peak_comparison.csv')
    assert len(shared) == 5 and all(truth(r['shared_raw_counts_identical']) for r in shared)
    selected = read(supplement / 'RNA_ATAC_response/RNA_selected_genes.csv')
    assert len(selected) == 15
    for pop in {r['population'] for r in selected}:
        # Independent check that ATAC p values did not select the RNA examples.
        rna = read(root / 'results/core1/tables' / ('DE_pca__'+pop+'.csv.gz'))
        eligible = [r for r in rna if r['padj'] not in {'NA',''} and float(r['padj']) < .05 and abs(float(r['log2FoldChange'])) >= 1]
        eligible.sort(key=lambda r: (float(r['padj']), -abs(float(r['log2FoldChange'])), r['gene']))
        actual = sorted([r for r in selected if r['population'] == pop], key=lambda r: int(r['RNA_rank']))
        assert [r['gene_id'] for r in actual] == [r['gene'] for r in eligible[:3]]
    genes = read(supplement / 'RNA_ATAC_response/RNA_ATAC_gene_summary.csv')
    long = read(supplement / 'RNA_ATAC_response/RNA_ATAC_peak_context.csv')
    assert len(genes) == 15 and len({(r['population'],r['gene_id'],r['peak']) for r in long}) == len(long)
    for g in genes:
        rows = [r for r in long if (r['population'],r['gene_id']) == (g['population'],g['gene_id'])]
        assert len(rows) == int(g['common_peaks'])
        assert sum(truth(r['ATAC_in_model']) for r in rows) == int(g['modelled_peaks'])
        assert sum(truth(r['ATAC_FDR_significant']) for r in rows) == int(g['nearby_FDR_peaks'])
        eligible = [r for r in rows if truth(r['ATAC_in_model'])]
        eligible.sort(key=lambda r: (not truth(r['overlaps_TSS_window']),float(r['TSS_distance_bp']),r['peak']))
        assert g['representative_peak'] == (eligible[0]['peak'] if eligible else 'NA')
    integrity = read(supplement / 'original_input_integrity.csv')
    assert all(truth(r['unchanged']) and r['md5_before'] == r['md5_after'] for r in integrity)
    validation = read(supplement / 'primary_model_validation.csv')
    assert len(validation) == 5 and all(truth(r['original_statistics_reproduced']) and truth(r['original_counts_match_pseudobulk']) for r in validation)
    record = json.loads((supplement / 'revision_record.json').read_text())
    for path,digest in record['source_sha256'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest() == digest, path
    return {'diagnostic_figures': 10, 'threshold_comparison_rows': 10,
            'RNA_examples': 15, 'unchanged_primary_models': 5}


if __name__ == '__main__':
    import os
    root = Path(os.environ.get('P3_ROOT', Path(__file__).resolve().parent.parent))
    print(json.dumps(verify(root), indent=2))
