# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
#!/usr/bin/env python3
"""Download three public ATAC SRA runs with resumable, checked byte ranges.

Only public downloads; no user data are uploaded. SRA validation is performed
separately by the upstream pipeline. Temporary chunk files are removed only
after their assembled SRA has passed the expected size check and been hashed.
"""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import urllib.request

ROOT = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
RUNS = {
    "SRR24921615": 5358583855,
    "SRR24921614": 5537979214,
    "SRR24921613": 5322812391,
}
CHUNKS = 4

def download_chunk(url, path, start, end):
    length = end - start + 1
    for attempt in range(30):
        existing = path.stat().st_size if path.exists() else 0
        assert existing <= length, (path, existing, length)
        if existing == length:
            return
        request = urllib.request.Request(url, headers={
            "Range": f"bytes={start + existing}-{end}",
            "User-Agent": "Core2-public-SRA-download/1.0",
        })
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                expected_range = f"bytes {start + existing}-{end}/"
                assert response.status == 206
                assert response.headers["Content-Range"].startswith(expected_range)
                with path.open("ab") as output:
                    while block := response.read(1024 * 1024):
                        output.write(block)
            assert path.stat().st_size == length
            return
        except Exception as exc:
            print(f"{path.name}: retry {attempt + 1}: {exc}", flush=True)
            time.sleep(min(5 + attempt, 30))
    raise RuntimeError(f"Download retries exhausted: {path}")

def download_run(run, expected_bytes):
    directory = ROOT / "work/raw_sra" / run
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / (run + ".sra")
    manifest = directory / "download.json"
    if destination.exists() and manifest.exists():
        assert destination.stat().st_size == expected_bytes
        print(f"{run}: reuse completed download", flush=True)
        return
    url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as response:
        remote_bytes = int(response.headers["Content-Length"])
        etag = response.headers.get("ETag")
    assert remote_bytes == expected_bytes, (run, remote_bytes, expected_bytes)
    print(f"{run}: downloading {expected_bytes:,} bytes in {CHUNKS} resumable chunks", flush=True)
    pieces = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=CHUNKS) as pool:
        pending = []
        for i in range(CHUNKS):
            start, end = expected_bytes * i // CHUNKS, expected_bytes * (i + 1) // CHUNKS - 1
            path = directory / f"chunk_{i}.part"
            pieces.append(path)
            pending.append(pool.submit(download_chunk, url, path, start, end))
        for future in pending:
            future.result()
    partial = directory / (run + ".sra.part")
    digest = hashlib.sha256()
    with partial.open("wb") as target:
        for path in pieces:
            with path.open("rb") as source:
                while block := source.read(8 * 1024 * 1024):
                    target.write(block)
                    digest.update(block)
    assert partial.stat().st_size == expected_bytes
    partial.replace(destination)
    manifest.write_text(json.dumps({"run": run, "url": url,
        "expected_bytes": expected_bytes, "bytes": destination.stat().st_size,
        "remote_etag": etag, "sha256": digest.hexdigest(),
        "note": "Byte-range checks passed; vdb-validate must also pass before processing."}, indent=2) + "\n")
    for path in pieces:
        path.unlink()
    print(f"{run}: DOWNLOAD COMPLETE, SHA-256 {digest.hexdigest()}", flush=True)

if __name__ == "__main__":
    assert shutil.disk_usage(ROOT).free > 40 * 1024**3
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(download_run, run, size) for run, size in RUNS.items()]
        for future in futures:
            future.result()
