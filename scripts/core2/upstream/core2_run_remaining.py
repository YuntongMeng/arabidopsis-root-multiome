# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
#!/usr/bin/env python3
"""Process the remaining three samples sequentially, stopping before nucleus QC."""
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
PYTHON = ROOT / "work/atac_env/bin/python"
PIPELINE = Path(__file__).resolve().parent / "core2_sample_upstream.py"
SAMPLES = [("control_rep2", "control2", "SRR24921615"),
           ("stress_rep1", "stress1", "SRR24921614"),
           ("stress_rep2", "stress2", "SRR24921613")]

for sample, tag, run in SAMPLES:
    marker = ROOT / "work/raw_sra" / run / "download.json"
    print(f"{sample}: waiting for checked download", flush=True)
    while not marker.exists():
        time.sleep(10)
    downloaded = json.loads(marker.read_text())
    assert downloaded["run"] == run
    destination = ROOT / "outputs" / (tag + "_upstream")
    destination.mkdir(parents=True, exist_ok=True)
    completion = destination / "upstream_completion.json"
    if completion.exists():
        assert json.loads(completion.read_text())["sample"] == sample
        print(f"{sample}: completed output already exists", flush=True)
        continue
    print(f"{sample}: starting validated upstream pipeline", flush=True)
    with (destination / "pipeline.log").open("a") as log:
        subprocess.run([str(PYTHON), "-u", str(PIPELINE), "all", "--sample", sample],
                       cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    assert completion.exists()
    print(f"{sample}: UPSTREAM COMPLETE BEFORE QC", flush=True)

print("ALL THREE SAMPLES COMPLETE; cross-sample verification remains", flush=True)
