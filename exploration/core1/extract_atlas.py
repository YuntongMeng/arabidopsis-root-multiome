"""Re-extract versioned marker panels from the unmodified Shahan 2022 Data S2.
Run with the project root as an optional argument. Requires pandas and openpyxl.
The workbook is already in data/core1_corrected_reference; no download occurs.
"""
import sys
from pathlib import Path
import pandas as pd
root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
dest = root / 'data/core1_corrected_reference'
p = dest / 'Shahan2022_DataS2.xlsx'
k = pd.read_excel(p, sheet_name='K) Atlas_markers_1', header=3)
k['source'] = 'Shahan et al. 2022 Dev Cell, Data S2K'
k['source_url'] = 'https://doi.org/10.1016/j.devcel.2022.01.008'
k.to_csv(dest / 'atlas_top50_signatures.csv', index=False)
frames = []
for sheet in ['F) Markers for SEMITONES', 'G) Markers for novoSpaRc']:
    d = pd.read_excel(p, sheet_name=sheet, header=1)
    d.columns = [x.strip() for x in d.columns]
    d['source_sheet'] = sheet
    d['source_url'] = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9014886/'
    frames.append(d)
pd.concat(frames).to_csv(dest / 'canonical_markers_provenance.csv', index=False)
j = pd.read_excel(p, sheet_name='J)Cell cycle and ploidy markers', header=3)
j.iloc[:, :3].to_csv(dest / 'cell_cycle_markers.csv', index=False)
