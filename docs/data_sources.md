# Core 1 inputs and provenance

The four processed Multiome matrices come from [GEO GSE235495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235495). The GEO archive name `RAW.tar` does not mean it contains unfiltered empty droplets.

| Sample | GEO sample | File stem | Directory |
| --- | --- | --- | --- |
| control_rep1 | GSM7504014 | GSM7504014_Control1 | data/rna/control_rep1/ |
| control_rep2 | GSM7504015 | GSM7504015_Control2 | data/rna/control_rep2/ |
| stress_rep1 | GSM7504016 | GSM7504016_Case1 | data/rna/stress_rep1/ |
| stress_rep2 | GSM7504017 | GSM7504017_Case2 | data/rna/stress_rep2/ |

Download the supplementary matrices from GEO. For each stem, place `*_matrix.mtx.gz`, `*_barcodes.tsv.gz`, and `*_features.tsv.gz` in the listed directory. The runner reads these original files; it does not download them automatically. `sample_metadata.tsv` is included.

The processed matrices contain both genes and ATAC peaks. **Only rows with feature column 3 exactly `Gene Expression` are retained.** Column 1 supplies gene IDs; column 2 is a label, not feature type. ATAC peak counts are excluded before QC, normalization and all later RNA steps.

P1 inputs under `data/project1/` contain the 137 interaction candidates and seven GO terms. `data/core1_corrected_reference/` includes the WT drought-DEG control, complete interaction results, reviewed partition-bound labels and marker provenance. The [Shahan et al. 2022 root atlas](https://doi.org/10.1016/j.devcel.2022.01.008) Data S2 workbook and extracted signature panel are retained; the extraction tool is under `exploration/core1/`. Canonical marker sources are listed row by row in `canonical_markers_provenance.csv`.

Original matrices and large RDS checkpoints are excluded from Git. After restoring the GEO inputs and installing dependencies, `all` computes new checkpoints from those matrices. The optional historical `diagnostics` stage needs a local old Seurat object and is excluded from `all`.

Current large CSV tables are published as lossless `.csv.gz` copies; smaller tables remain CSV. A rerun writes ordinary CSV locally. Read a compressed published table with `utils::read.csv(gzfile("path/to/table.csv.gz"))`. The two full frozen tables needed by script 07 are also distributed in gzip form, and its reader accepts either form. The frozen originals have not been edited.

The [core1-v1 release](https://github.com/YuntongMeng/arabidopsis-root-multiome/releases/tag/core1-v1) is a historical mixed-input snapshot, **not the corrected analysis**. Do not restore its scripts/results over the current workflow. Its numerical conclusions have been superseded; it remains available for the decision trail.
