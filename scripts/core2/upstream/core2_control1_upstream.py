# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
#!/usr/bin/env python3
"""Reconstruct Control 1 ATAC fragments and common-peak counts; no nucleus QC.

Run with the isolated Python in work/atac_env/bin/python. Each stage is logged.
The original project and downloaded SRA are read-only inputs. Output is confined
to this task's outputs/control1_upstream and work/control1_upstream directories.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone

ROOT = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
WORK = ROOT / "work/control1_upstream"
OUT = ROOT / "outputs/control1_upstream"
LOG = OUT / "logs"
TOOLS = Path(os.environ.get("P3_CORE2_TOOLS", str(ROOT / "work/atac_env/bin")))
SRA_TOOLS = Path(os.environ.get("P3_SRA_TOOLS", str(ROOT / "work/tools/sratoolkit.3.4.1-mac-arm64/bin")))
PROJECT = Path(os.environ.get("P3_ROOT", os.getcwd())).expanduser().resolve()
RUN = "SRR24921616"
SRA = PROJECT / "data/atac/raw" / RUN / (RUN + ".sra")
EXPECTED_SPOTS = 130436487
EXPECTED_BYTES = 5679168369
FASTQ = WORK / "fastq"
REF = WORK / "reference"
ENV = dict(os.environ, NCBI_HOME=str(ROOT / "work/ncbi"))

def stamp():
    return datetime.now(timezone.utc).isoformat()

def setup():
    for path in (WORK, OUT, LOG, FASTQ, REF):
        path.mkdir(parents=True, exist_ok=True)

def run(command, name):
    command = [str(x) for x in command]
    with (LOG / "commands.jsonl").open("a") as handle:
        handle.write(json.dumps({"time": stamp(), "name": name, "argv": command}) + "\n")
    print(stamp(), name, "started", flush=True)
    with (LOG / (name + ".log")).open("w") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=ENV, cwd=WORK)
    if result.returncode:
        raise RuntimeError(f"{name}: exit {result.returncode}; see {LOG / (name + '.log')}")
    print(stamp(), name, "completed", flush=True)

def save_json(path, obj):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2) + "\n")
    temp.replace(path)

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def export():
    assert SRA.stat().st_size == EXPECTED_BYTES
    done = WORK / "export.complete.json"
    if done.exists():
        assert all((FASTQ / f"{RUN}_{i}.fastq").is_file() for i in range(1, 5))
        print("Reuse completed full FASTQ export", flush=True)
        return
    assert not list(FASTQ.glob("*.fastq")), "Partial export found: inspect before retrying"
    assert shutil.disk_usage(WORK).free > 150 * 1024**3
    temp = WORK / "sra_tmp"
    temp.mkdir(exist_ok=True)
    run([SRA_TOOLS / "fasterq-dump", SRA, "--split-files", "--include-technical",
         "-M", "1", "-e", "4", "-m", "1G", "-O", FASTQ, "-t", temp,
         "--seq-defline", "@$ac.$si", "--qual-defline", "+", "--details"], "01_export")
    paths = [FASTQ / f"{RUN}_{i}.fastq" for i in range(1, 5)]
    assert all(p.exists() and p.stat().st_size > 0 for p in paths)
    save_json(done, {"time": stamp(), "input_bytes": SRA.stat().st_size,
                     "files_bytes": {p.name: p.stat().st_size for p in paths},
                     "records_expected_per_stream": EXPECTED_SPOTS,
                     "note": "Full synchronization/record validation happens in prepare stage."})

def reference():
    url = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
    archive = REF / url.rsplit("/", 1)[-1]
    fasta = REF / "TAIR10.fa"
    if not archive.exists():
        partial = archive.with_suffix(archive.suffix + ".part")
        run(["/usr/bin/curl", "--fail", "--location", "--silent", "--show-error",
             "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "60",
             "--retry", "10", "--retry-all-errors", "--retry-delay", "5",
             "--continue-at", "-", "--output", partial, url], "02_download_reference")
        partial.replace(archive)
    if not fasta.exists():
        with gzip.open(archive, "rb") as source, fasta.open("wb") as target:
            shutil.copyfileobj(source, target)
    run([TOOLS / "samtools", "faidx", fasta], "02_reference_faidx")
    lengths = dict((line.split("\t")[0], int(line.split("\t")[1]))
                   for line in Path(str(fasta) + ".fai").read_text().splitlines())
    expected = {"1":30427671, "2":19698289, "3":23459830, "4":18585056, "5":26975502}
    assert all(lengths.get(k) == v for k, v in expected.items()), lengths
    assert "Mt" in lengths and "Pt" in lengths, lengths
    for line in (ROOT / "outputs/core2_common_peaks.bed").read_text().splitlines():
        chrom, start, end = line.split("\t")[:3]
        assert 0 <= int(start) < int(end) <= lengths[chrom]
    index = REF / "TAIR10.chromap.index"
    if not index.exists():
        run([TOOLS / "chromap", "-i", "-r", fasta, "-o", index], "02_reference_index")
    save_json(OUT / "reference_provenance.json", {
        "url": url, "assembly": "TAIR10", "source_release": 62,
        "download_sha256": sha256(archive), "fasta_sha256": sha256(fasta),
        "contig_lengths": lengths, "note": "Includes Mt and Pt for competitive alignment; common peaks are nuclear chromosomes 1–5."
    })

def fastq_records(path):
    with path.open() as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
            assert header.startswith("@") and plus.startswith("+")
            assert sequence and len(sequence) == len(quality)
            yield header.split()[0], sequence, quality

def prepare():
    reference_dir = ROOT / "work/barcode_reference"
    expected_hashes = {"atac": "9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26",
                       "gex": "1279bb7ba28ddcd805e2f16ebd5d6933b7b1e7f128e18aabfe6ae8cd7b7bf710"}
    references = {}
    for kind in ("atac", "gex"):
        path = reference_dir / f"{kind}_737K-arc-v1.txt.gz"
        assert sha256(path) == expected_hashes[kind]
        with gzip.open(path, "rt") as handle:
            references[kind] = handle.read().splitlines()
        values = references[kind]
        assert len(values) == len(set(values)) == 736320
        assert all(len(s) == 16 and set(s) <= set("ACGT") for s in values)
    barcode_paths = list((PROJECT / "data/rna/control_rep1").glob("*barcodes.tsv.gz"))
    assert len(barcode_paths) == 1
    with gzip.open(barcode_paths[0], "rt") as handle:
        cells = handle.read().splitlines()
    assert len(cells) == len(set(cells)) == 3431
    gex_set = set(references["gex"])
    assert all(s.endswith("-1") and s[:-2] in gex_set for s in cells)
    (OUT / "matrix_barcodes.tsv").write_text("\n".join(cells) + "\n")
    (WORK / "atac_whitelist.txt").write_text("\n".join(references["atac"]) + "\n")
    with (WORK / "barcode_translation.tsv").open("w") as handle:
        for atac, gex in zip(references["atac"], references["gex"]):
            handle.write(f"{gex}-1\t{atac}\n")
    preview_dir = SRA.parent / "preview_10000"
    records = [list(fastq_records(preview_dir / f"{RUN}_{i}.fastq")) for i in range(1, 5)]
    assert all(len(x) == 10000 for x in records)
    for entries in zip(*records):
        assert len(set(x[0] for x in entries)) == 1
        assert [len(x[1]) for x in entries] == [8, 24, 50, 50]
    translation = dict(zip(references["atac"], references["gex"]))
    existing = set(s[:-2] for s in cells)
    complement = str.maketrans("ACGTN", "TGCAN")
    modes = []
    for mode in ("first16", "last16", "rc_first16", "rc_last16"):
        barcodes = [x[1][:16] if "first" in mode else x[1][-16:] for x in records[1]]
        if mode.startswith("rc_"):
            barcodes = [s.translate(complement)[::-1] for s in barcodes]
        translated = [translation.get(s) for s in barcodes]
        modes.append({"mode": mode, "exact_whitelist_matches": sum(x is not None for x in translated),
                      "records_matching_matrix": sum(x in existing for x in translated),
                      "distinct_matrix_nuclei": len(set(translated) & existing)})
    winner = next(x for x in modes if x["mode"] == "rc_last16")
    assert winner["exact_whitelist_matches"] > 8000
    assert winner["exact_whitelist_matches"] > 10 * max(x["exact_whitelist_matches"] for x in modes if x is not winner)
    save_json(OUT / "barcode_validation.json", {
        "sample": "control_rep1", "records": 10000, "streams_synchronized": True,
        "read_lengths": [8, 24, 50, 50], "modes": modes,
        "selected_format": "r1:0:-1,r2:0:-1,bc:8:23:-",
        "matrix_nuclei": len(cells),
        "barcode_source": json.loads((reference_dir / "sources.json").read_text()),
        "limits": "First 10,000 spots, not random. Exact match pilot, not a nucleus QC retention estimate."
    })
    print(json.dumps(modes, indent=2), flush=True)

def mapping_command(read_dir, output, threads):
    return [TOOLS / "chromap", "--preset", "atac", "-r", REF / "TAIR10.fa",
            "-x", REF / "TAIR10.chromap.index", "-1", read_dir / f"{RUN}_3.fastq",
            "-2", read_dir / f"{RUN}_4.fastq", "-b", read_dir / f"{RUN}_2.fastq",
            "--read-format", "r1:0:-1,r2:0:-1,bc:8:23:-",
            "--barcode-whitelist", WORK / "atac_whitelist.txt",
            "--barcode-translate", WORK / "barcode_translation.tsv",
            "--bc-error-threshold", "1", "--bc-probability-threshold", "0.9",
            "-q", "30", "-l", "2000", "-t", str(threads), "-o", output]

def pilot():
    output = WORK / "pilot.fragments.tsv"
    run(mapping_command(SRA.parent / "preview_10000", output, 4), "03_mapping_pilot")
    gex = set()
    with gzip.open(ROOT / "work/barcode_reference/gex_737K-arc-v1.txt.gz", "rt") as handle:
        gex = set(s.strip() + "-1" for s in handle)
    rows = output.read_text().splitlines()
    assert len(rows) > 1000
    contigs = Counter()
    for line in rows:
        chrom, start, end, barcode, multiplicity = line.split("\t")
        assert int(start) >= 0 and int(end) > int(start)
        assert barcode in gex and int(multiplicity) >= 1
        contigs[chrom] += 1
    # Independent oracle: explicitly reverse-complement barcode bases AND
    # reverse their quality scores, rather than relying on --read-format.
    manual_barcode = WORK / "pilot.manual_barcode.fastq"
    complement = str.maketrans("ACGTN", "TGCAN")
    with manual_barcode.open("w") as handle:
        for header, sequence, quality in fastq_records(SRA.parent / "preview_10000" / f"{RUN}_2.fastq"):
            bases = sequence[8:24].translate(complement)[::-1]
            quals = quality[8:24][::-1]
            handle.write(f"{header}\n{bases}\n+\n{quals}\n")
    oracle_output = WORK / "pilot.manual_barcode.fragments.tsv"
    command = mapping_command(SRA.parent / "preview_10000", oracle_output, 4)
    command[command.index("-b") + 1] = manual_barcode
    command[command.index("--read-format") + 1] = "r1:0:-1,r2:0:-1,bc:0:-1"
    run(command, "03_barcode_orientation_oracle")
    assert sorted(rows) == sorted(oracle_output.read_text().splitlines()), "Native barcode extraction disagrees with base/quality reverse-complement oracle"
    save_json(OUT / "alignment_pilot.json", {"records": 10000, "unique_fragments": len(rows),
               "translated_barcodes_valid": True, "contigs": contigs,
               "native_barcode_format_equals_explicit_base_quality_oracle": True,
               "note": "Small deterministic technical pilot, not biological QC."})
    print("Pilot fragments:", len(rows), "contigs:", dict(contigs), flush=True)

def validate_fastq():
    assert (WORK / "export.complete.json").exists()
    binary = WORK / "validate_fastq"
    run([os.environ.get("P3_CXX", "clang++"), "-O3", "-std=c++17", Path(__file__).resolve().parent / "core2_validate_fastq.cpp",
         "-o", binary], "04_compile_validator")
    run([binary, FASTQ, RUN, str(EXPECTED_SPOTS)], "04_validate_full_fastq")
    report = json.loads((LOG / "04_validate_full_fastq.log").read_text())
    save_json(OUT / "full_fastq_validation.json", report)

def map_full():
    assert (OUT / "full_fastq_validation.json").exists()
    assert (OUT / "alignment_pilot.json").exists()
    output = WORK / "control1.fragments.unsorted.tsv"
    done = WORK / "mapping.complete.json"
    if not done.exists():
        assert not output.exists(), "Partial mapping output found; inspect before retrying"
        run(mapping_command(FASTQ, output, 6), "05_map_full")
        assert output.stat().st_size > 0
        save_json(done, {"time": stamp(), "output_bytes": output.stat().st_size})
    else:
        print("Reuse completed full mapping", flush=True)

def finalize():
    assert (WORK / "mapping.complete.json").exists()
    done = WORK / "finalize.complete.json"
    if done.exists():
        print("Reuse validated, indexed fragments", flush=True)
        return
    source = WORK / "control1.fragments.unsorted.tsv"
    sorted_path = WORK / "control1.fragments.sorted.tsv"
    sort_env = dict(ENV, LC_ALL="C")
    command = ["/usr/bin/sort", "-S", "1G", "-T", str(WORK), "-k1,1", "-k2,2n",
               "-k3,3n", "-k4,4", "-o", str(sorted_path), str(source)]
    with (LOG / "commands.jsonl").open("a") as handle:
        handle.write(json.dumps({"time": stamp(), "name": "06_sort", "argv": command, "LC_ALL": "C"}) + "\n")
    print(stamp(), "06_sort started", flush=True)
    with (LOG / "06_sort.log").open("w") as log:
        subprocess.run(command, env=sort_env, cwd=WORK, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(stamp(), "06_sort completed", flush=True)
    lengths = json.loads((OUT / "reference_provenance.json").read_text())["contig_lengths"]
    with gzip.open(ROOT / "work/barcode_reference/gex_737K-arc-v1.txt.gz", "rt") as handle:
        whitelist = set(x.strip() + "-1" for x in handle)
    matrix_cells = set((OUT / "matrix_barcodes.tsv").read_text().splitlines())
    contigs, fragments_per_barcode = Counter(), Counter()
    raw_rows = unique_rows = support_total = merged_rows = 0
    destination = OUT / "control1.fragments.tsv.gz"
    partial = destination.with_suffix(destination.suffix + ".part")
    with partial.open("wb") as target, (LOG / "06_bgzip.log").open("w") as log:
        compressor = subprocess.Popen([str(TOOLS / "bgzip"), "-@", "4", "-c"],
                                      stdin=subprocess.PIPE, stdout=target, stderr=log)
        previous = None
        copies = 0
        def emit(key, multiplicity):
            nonlocal unique_rows
            chrom, start, end, barcode = key
            compressor.stdin.write(f"{chrom}\t{start}\t{end}\t{barcode}\t{multiplicity}\n".encode())
            unique_rows += 1
            contigs[chrom] += 1
            fragments_per_barcode[barcode] += 1
        try:
            with sorted_path.open() as source_handle:
                for line in source_handle:
                    fields = line.rstrip("\n").split("\t")
                    assert len(fields) == 5
                    chrom, start_s, end_s, barcode, copies_s = fields
                    start, end, multiplicity = int(start_s), int(end_s), int(copies_s)
                    assert chrom in lengths and 0 <= start < end <= lengths[chrom]
                    assert barcode in whitelist and multiplicity >= 1
                    key = (chrom, start, end, barcode)
                    assert previous is None or key >= previous, "Fragments not sorted"
                    raw_rows += 1
                    support_total += multiplicity
                    if key == previous:
                        copies += multiplicity
                        merged_rows += 1
                    else:
                        if previous is not None:
                            emit(previous, copies)
                        previous, copies = key, multiplicity
            if previous is not None:
                emit(previous, copies)
            compressor.stdin.close()
            assert compressor.wait() == 0
        except BaseException:
            compressor.kill()
            compressor.wait()
            raise
    assert unique_rows > 1000000 and all(contigs.get(str(x), 0) > 0 for x in range(1, 6))
    partial.replace(destination)
    run([TOOLS / "bgzip", "-t", destination], "06_bgzip_integrity")
    run([TOOLS / "tabix", "-p", "bed", destination], "06_tabix_index")
    run([TOOLS / "tabix", "-l", destination], "06_tabix_contigs")
    indexed = set((LOG / "06_tabix_contigs.log").read_text().splitlines())
    assert set(str(x) for x in range(1, 6)) <= indexed
    absent = sorted(matrix_cells - fragments_per_barcode.keys())
    report = {"time": stamp(), "input_mapping_rows": raw_rows, "unique_fragments": unique_rows,
              "read_pair_support": support_total, "additional_identical_rows_collapsed": merged_rows,
              "fragments_by_contig": dict(contigs), "barcodes_with_fragments": len(fragments_per_barcode),
              "existing_matrix_nuclei": len(matrix_cells), "matrix_nuclei_with_fragments": len(matrix_cells) - len(absent),
              "matrix_nuclei_without_fragments": absent, "fragments_sha256": sha256(destination),
              "index_sha256": sha256(Path(str(destination) + ".tbi")),
              "bed_coordinates": "0-based half-open, Tn5 shifted once by Chromap ATAC preset",
              "barcode_format": "paired GEX barcode plus -1; all valid-whitelist barcodes retained",
              "scope": "Upstream technical validation only; no nucleus QC thresholds, no cell calling."}
    save_json(OUT / "fragment_validation.json", report)
    save_json(done, {"time": stamp(), "output": str(destination)})
    print(json.dumps(report, indent=2), flush=True)

def count():
    assert (WORK / "finalize.complete.json").exists()
    run([os.environ.get("P3_RSCRIPT", "Rscript"), Path(__file__).resolve().parent / "core2_control1_count_peaks.R", ROOT], "07_count_peaks")
    report()

def report():
    assert (OUT / "count_matrix_manifest.csv").exists()
    versions = {}
    for name, command in {
        "chromap": [TOOLS / "chromap", "--version"],
        "fasterq_dump": [SRA_TOOLS / "fasterq-dump", "--version"],
        "samtools": [TOOLS / "samtools", "--version"],
        "bgzip": [TOOLS / "bgzip", "--version"],
        "tabix": [TOOLS / "tabix", "--version"],
    }.items():
        result = subprocess.run([str(x) for x in command], capture_output=True, text=True, env=ENV)
        assert result.returncode == 0
        versions[name] = (result.stdout + result.stderr).strip()
    artifacts = [OUT / "control1.fragments.tsv.gz", OUT / "control1.fragments.tsv.gz.tbi",
                 OUT / "control1_common_peak_counts.rds", OUT / "common_peaks.rds",
                 OUT / "common_peaks.tsv", OUT / "matrix_barcodes.tsv"]
    for path in artifacts:
        assert path.is_file() and path.stat().st_size > 0
    source_scripts = [Path(__file__).resolve().parent / name for name in
                      ("core2_control1_upstream.py", "core2_control1_count_peaks.R",
                       "core2_validate_fastq.cpp", "core2_control1_run.R")]
    manifest = {
        "status": "upstream_complete_before_nucleus_QC", "completed_utc": stamp(),
        "sample": "control_rep1", "sra_accession": RUN,
        "sra_path": str(SRA), "sra_bytes": SRA.stat().st_size, "sra_sha256": sha256(SRA),
        "sra_public_url": f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{RUN}/{RUN}",
        "software": versions,
        "method": {"aligner": "Chromap; an ATAC-only reconstruction, not an identical Cell Ranger ARC rerun",
                   "reference": "TAIR10 nuclear chromosomes plus Mt and Pt",
                   "preset": "atac", "MAPQ_min": 30, "max_insert_size": 2000,
                   "barcode_hamming_max": 1, "barcode_correction_probability_min": 0.9,
                   "barcode_format": "reverse complement of I2 bases 9-24 (1-based), then paired ATAC-to-GEX translation",
                   "duplicates": "cell-level deduplication; identical final BED coordinate+barcode rows additionally collapsed",
                   "Tn5_shift": "once, inside Chromap preset: +4 start / -5 end",
                   "count_unit": "one unique fragment per overlapping peak; duplicate_count not multiplied",
                   "common_peak_universe": "union/reduce of original four peak sets, 27677 intervals",
                   "count_matrix_barcodes": "all 3431 barcodes from the supplied Control 1 matrix, original order retained",
                   "nucleus_QC": "not performed; no TSS/FRiP/nucleosome thresholds, no doublet filtering, no clustering"},
        "artifacts": [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in artifacts],
        "scripts": [{"path": str(p), "sha256": sha256(p)} for p in source_scripts],
        "sources": {
            "chromap": "https://github.com/haowenz/chromap",
            "chromap_paper": "https://doi.org/10.1038/s41467-021-26865-w",
            "sra_export": "https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump",
            "signac_common_peaks": "https://stuartlab.org/signac/articles/merging",
            "barcode_pairing": "https://scg-lib-structs.readthedocs.io/en/latest/multi/10x_multiome.html"},
        "limits": ["Mapping and counting differ from the original Cell Ranger ARC pipeline; original cut-site counts are not numerically interchangeable with new fragment-overlap counts.",
                   "Paired barcode references are from the documented Teichmann laboratory mirror; hashes and provenance saved separately.",
                   "Existing matrix barcodes were supplied by the study; retaining them is not new cell calling and does not certify nucleus quality.",
                   "Reconstructed fragments preserve organellar alignments; the count matrix is restricted to the common nuclear peak universe.",
                   "No other sample was processed and no Core 1 results were modified."]
    }
    save_json(OUT / "upstream_completion.json", manifest)
    print("UPSTREAM COMPLETE; nucleus QC remains pending:", OUT, flush=True)

def main():
    setup()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["export", "reference", "prepare", "pilot", "validate_fastq", "map_full", "finalize", "count", "report", "all"])
    args = parser.parse_args()
    if args.stage == "all":
        for stage in (reference, prepare, export, validate_fastq, pilot, map_full, finalize, count):
            stage()
    else:
        globals()[args.stage]()

if __name__ == "__main__":
    main()
