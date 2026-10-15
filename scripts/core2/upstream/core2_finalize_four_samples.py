# Repository edition: paths/tools configured by P3_ROOT and P3_CORE2_WORK.
#!/usr/bin/env python3
"""Audit saved artifacts and record the completed four-sample pre-QC handoff."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from datetime import datetime, timezone

ROOT = Path(os.environ["P3_CORE2_WORK"]).expanduser().resolve()
OUT = ROOT / "outputs"
SAMPLES = [("control_rep1", "control1", "SRR24921616", 130436487),
           ("control_rep2", "control2", "SRR24921615", 129424454),
           ("stress_rep1", "stress1", "SRR24921614", 133450323),
           ("stress_rep2", "stress2", "SRR24921613", 129510415)]

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            h.update(block)
    return h.hexdigest()

metadata_path = ROOT / "work/remaining_run_metadata.tsv"
with metadata_path.open() as handle:
    public = {x["run_accession"]: x for x in csv.DictReader(handle, delimiter="\t")}
method_keys = ["aligner", "reference", "preset", "MAPQ_min", "max_insert_size",
               "barcode_hamming_max", "barcode_correction_probability_min", "barcode_format",
               "duplicates", "Tn5_shift", "count_unit", "common_peak_universe", "nucleus_QC"]
baseline = None
samples = []
for sample, tag, run, spots in SAMPLES:
    directory = OUT / (tag + "_upstream")
    manifest = json.loads((directory / "upstream_completion.json").read_text())
    assert manifest["sample"] == sample and manifest["sra_accession"] == run
    assert manifest["status"] == "upstream_complete_before_nucleus_QC"
    assert int(public[run]["read_count"]) == spots
    condition = "control" if sample.startswith("control") else "case"
    rep = sample[-1]
    assert f"{condition},rep{rep},scATAC-seq" in public[run]["experiment_title"]
    full = json.loads((directory / "full_fastq_validation.json").read_text())
    # Recheck the recorded full-stream validation, in addition to its exit code.
    assert full["synchronized"] and full["every_spot_id_checked"]
    assert full["sequence_quality_lengths_checked"]
    assert len(full["streams"]) == 4
    assert [x["records"] for x in full["streams"]] == [spots] * 4
    assert [x["sequence_length"] for x in full["streams"]] == [8, 24, 50, 50]
    method = {k: manifest["method"][k] for k in method_keys}
    if baseline is None:
        baseline = method
    assert method == baseline, (sample, "Method mismatch")
    for artifact in manifest["artifacts"]:
        path = Path(artifact["path"])
        assert path.stat().st_size == artifact["bytes"]
        assert digest(path) == artifact["sha256"], (sample, path)
    for script in manifest["scripts"]:
        assert digest(Path(script["path"])) == script["sha256"], (sample, "Script changed after run")
    assert digest(directory / "common_peaks.tsv") == digest(OUT / "control1_upstream/common_peaks.tsv")
    fragments = json.loads((directory / "fragment_validation.json").read_text())
    with (directory / "count_matrix_manifest.csv").open() as handle:
        matrix = next(csv.DictReader(handle))
    assert int(matrix["peaks"]) == 27677 and matrix["nucleus_QC_applied"] == "FALSE"
    samples.append({"sample": sample, "run": run, "spots": spots,
        "common_peaks": 27677, "nuclei_before_QC": int(matrix["matrix_nuclei_retained"]),
        "unique_fragments_all_whitelist_barcodes_including_organelles": fragments["unique_fragments"],
        "nonzero_matrix_entries": int(matrix["nonzero_entries"]),
        "independent_peak_rows_checked": int(matrix["independent_peak_rows_checked"]),
        "fragments": str(directory / (tag + ".fragments.tsv.gz")),
        "matrix": str(directory / (tag + "_common_peak_counts.rds")),
        "checksums_verified": True, "nucleus_QC_applied": False})

subprocess.run([os.environ.get("P3_RSCRIPT", "Rscript"), str(Path(__file__).resolve().parent / "core2_verify_four_matrices.R"), str(ROOT)], check=True)
handoff = {"status": "all_four_samples_upstream_complete_before_nucleus_QC",
           "completed_utc": datetime.now(timezone.utc).isoformat(),
           "method_identical_across_samples": True,
           "common_peak_coordinates_and_order_identical": True,
           "original_barcode_order_preserved_per_sample": True,
           "public_metadata_url": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP443792&result=read_run&fields=run_accession,experiment_title,read_count&format=tsv",
           "samples": samples,
           "next_step": "Nucleus QC only after reviewing metrics; no filtering or biological analysis performed."}
(OUT / "core2_four_sample_upstream_completion.json").write_text(json.dumps(handoff, indent=2) + "\n")
(OUT / "core2_run_identity_metadata.tsv").write_text(metadata_path.read_text())
with (OUT / "core2_four_sample_upstream_summary.csv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(samples[0]))
    writer.writeheader()
    writer.writerows(samples)

lines = ["# Core 2：四样本 upstream 交付", "",
    "已完成四个 ATAC 样本的 fragments 重建和共同区域重新计数，停在正式细胞核 QC 之前。", "",
    "| 样本 | SRA | 共同区域数 | 原矩阵核数（全部保留） | 全部有效 barcode 的去重片段数* |",
    "|---|---|---:|---:|---:|"]
for x in samples:
    lines.append(f'| {x["sample"]} | {x["run"]} | 27,677 | {x["nuclei_before_QC"]:,} | {x["unique_fragments_all_whitelist_barcodes_including_organelles"]:,} |')
lines += ["", "*片段文件包含核染色体、Mt/Pt 和原矩阵以外的有效 barcode；该列不是保留细胞核数。", "",
    "## 处理与核对", "",
    "- 四个样本使用同一 TAIR10 序列/索引、配对 barcode 表、Chromap ATAC 参数和 fragment-overlap 计数规则。",
    "- 另外三个 SRA 下载后通过 vdb-validate；全量四路 FASTQ 逐条验证。每个样本单独验证 barcode 提取方向和配对转换。",
    "- Fragments 全量检查坐标、barcode、重复和排序；BGZF 完整性及 tabix 索引检查通过。",
    "- 每份矩阵通过已知计数的小例子和至少 40 个真实区域 × 全部核列的独立重新计数检查。",
    "- 在新的 R 会话中重读四份矩阵：27,677 行的坐标和顺序完全一致；各样本原始 barcode 列及顺序完整保留。",
    "- 最终文件 SHA-256 与各样本完成清单一致；原 Control 1 结果及原项目/Core 1 输出未修改。", "",
    "## 文件位置", "",
    "各样本目录为 control1_upstream、control2_upstream、stress1_upstream、stress2_upstream。每个目录包含：", "",
    "- `<样本>.fragments.tsv.gz` 和 `.tbi`：后续 QC 的片段文件及索引。",
    "- `<样本>_common_peak_counts.rds`：共同 peak × 原矩阵 barcode 的计数矩阵。",
    "- `common_peaks.rds`、`common_peaks.tsv`、`matrix_barcodes.tsv`：坐标与核身份。",
    "- `upstream_completion.json`、检查清单和 `logs/`：来源、版本、参数、验证与文件指纹。", "",
    "汇总：core2_four_sample_upstream_summary.csv；重新读盘的矩阵核对：core2_four_sample_matrix_checks.csv。", "",
    "## 下一次从这里继续", "",
    "计算每个核的片段数量、TSS enrichment、核小体信号、FRiP 和细胞器比例，再根据四样本分布制定 QC。尚未执行这些指标的核筛选、doublet 去除、TF-IDF/LSI、UMAP 或差异分析。",
    "合并样本时使用样本名 + barcode 标识核，不能把跨样本同名 barcode 当成同一个核。", "",
    "## 方法边界", "",
    "这是 Chromap 的 ATAC-only 重建，不能称为原 Cell Ranger ARC 全流程的完全复现。新矩阵按去重 fragment 与区域重叠计数，原 ARC 切点计数及旧阈值不可直接照搬。一个片段可重叠多个区域；精确 FRiP 必须按片段是否命中任一 peak 统计，不能简单把矩阵列和当作唯一片段数。",
    "本次保留研究提供的 barcode 集合，未重新做 cell calling；保留不代表这些核已通过 QC。原始 SRA 和 FASTQ 中间文件目前均保留。", "",
    "Control 1 简明说明见 core2_control1_upstream_guide.pdf 和同名 Markdown。新增三样本执行脚本为 core2_sample_upstream.py、core2_count_common_peaks.R；下载和顺序启动脚本为 core2_download_remaining.py、core2_run_remaining.py。", ""]
(OUT / "core2_four_sample_upstream_README.md").write_text("\n".join(lines))
print("ALL FOUR SAMPLES VERIFIED; STOPPED BEFORE NUCLEUS QC", flush=True)
