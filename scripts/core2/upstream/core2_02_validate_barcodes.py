# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
"""Validate the first 1,000 Control 1 SRA spots against paired ARC barcode lists.

Inputs under work/ were obtained from the public SRA and the Teichmann laboratory
mirror. This is an exact-match pilot, not cell calling or quality filtering.
"""
import os
from pathlib import Path
from collections import Counter
import csv
import gzip
import hashlib
import json

ROOT = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
REFERENCE = ROOT / "work/barcode_reference"
PROJECT = Path(os.environ.get("P3_ROOT", os.getcwd())).expanduser().resolve()

gex = gzip.decompress((REFERENCE / "gex_737K-arc-v1.txt.gz").read_bytes()).decode().splitlines()
atac = gzip.decompress((REFERENCE / "atac_737K-arc-v1.txt.gz").read_bytes()).decode().splitlines()
assert len(gex) == len(atac) == 736320
assert len(set(gex)) == len(gex) and len(set(atac)) == len(atac)
assert all(len(x) == 16 and set(x) <= set("ACGT") for x in gex + atac)
translation = dict(zip(atac, gex))
barcode_file = list((PROJECT / "data/rna/control_rep1").glob("*barcodes.tsv.gz"))
assert len(barcode_file) == 1
with gzip.open(barcode_file[0], "rt") as handle:
    existing = {line.strip().split("-")[0] for line in handle}

preview = ROOT / "work/control1_reads_preview.tsv"
spots = []
spacers = Counter()
for line in preview.read_text().splitlines():
    sequence, lengths, starts, read_types = line.split("\t")
    lengths = [int(x) for x in lengths.split(",")]
    starts = [int(x) for x in starts.split(",")]
    assert lengths == [8, 24, 50, 50] and starts == [0, 8, 32, 82]
    assert len(sequence) == 132
    index = sequence[8:32]
    spacers[index[:8]] += 1
    spots.append(index)
assert len(spots) == 1000

complement = str.maketrans("ACGTN", "TGCAN")
rows = []
for mode in ["first16", "last16", "rc_first16", "rc_last16"]:
    valid = matched = 0
    matched_nuclei = set()
    for index in spots:
        barcode = index[:16] if "first" in mode else index[-16:]
        if mode.startswith("rc_"):
            barcode = barcode.translate(complement)[::-1]
        gex_barcode = translation.get(barcode)
        if gex_barcode is not None:
            valid += 1
            if gex_barcode in existing:
                matched += 1
                matched_nuclei.add(gex_barcode)
    rows.append(dict(
        sample="control_rep1", run="SRR24921616", sra_spots_examined=len(spots),
        extraction_mode=mode, exact_atac_whitelist_matches=valid,
        spots_matching_existing_nuclei=matched,
        distinct_existing_nuclei_matched=len(matched_nuclei),
    ))

output = ROOT / "outputs/core2_control1_barcode_pilot.csv"
with output.open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
provenance = {
    "run": "SRR24921616", "spot_range": "1-1000 (not random)",
    "tool": "NCBI SRA Toolkit 3.4.1 vdb-dump",
    "source": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR24921616/SRR24921616",
    "columns": "READ,READ_LEN,READ_START,READ_TYPE", "table": "SEQUENCE",
    "preview_sha256": hashlib.sha256(preview.read_bytes()).hexdigest(),
    "barcode_references": json.loads((REFERENCE / "sources.json").read_text()),
    "existing_barcodes_in_gex_reference": len(existing & set(gex)),
    "existing_barcodes_total": len(existing), "leading_spacer_counts": dict(spacers),
    "limits": "Exact-match pilot only. No barcode error correction, alignment, fragment generation, cell calling or QC retention estimate. Reference lists are a documented laboratory mirror, not obtained directly from the 10x binary distribution.",
}
(ROOT / "outputs/core2_control1_barcode_pilot_provenance.json").write_text(json.dumps(provenance, indent=2))
for row in rows:
    print(row)
