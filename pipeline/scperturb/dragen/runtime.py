#!/usr/bin/env python3
"""Run a DRAGEN command in the rootless EL9 runtime."""

from __future__ import annotations

import os
import subprocess
import gzip
import csv
from pathlib import Path


def run_dragen(config: dict, command: list[str], *, stage: Path, intermediate: Path, capture: bool = False) -> subprocess.CompletedProcess:
    """Run the image-bundled DRAGEN binary directly; never launch a nested container."""
    intermediate.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    runtime = Path(config["tmp_root"]) / "dragen_runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    env.update({"HOME": str(runtime / "home"), "TMPDIR": str(runtime), "TMP": str(runtime),
                "TEMP": str(runtime), "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8"})
    Path(env["HOME"]).mkdir(parents=True, exist_ok=True)
    for key in ("http_proxy", "https_proxy", "all_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"):
        env.pop(key, None)
    env["NO_PROXY"] = "*"
    env["no_proxy"] = "*"
    return subprocess.run(command, env=env, cwd=config["project_root"], check=False, text=True,
                          stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.STDOUT if capture else None)


def make_fastq_subset(read1: Path, read2: Path, out1: Path, out2: Path, limit: int) -> int:
    out1.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(read1, "rt", encoding="ascii") as left, gzip.open(read2, "rt", encoding="ascii") as right, \
         gzip.open(out1, "wt", encoding="ascii") as left_out, gzip.open(out2, "wt", encoding="ascii") as right_out:
        while count < limit:
            a = [left.readline() for _ in range(4)]
            b = [right.readline() for _ in range(4)]
            if not a[0] and not b[0]:
                break
            if not a[0] or not b[0] or a[0].split()[0] != b[0].split()[0]:
                raise ValueError(f"paired FASTQ records are out of sync at record {count + 1}")
            if any(not value for value in a + b):
                raise ValueError(f"truncated FASTQ record at record {count + 1}")
            left_out.writelines(a)
            right_out.writelines(b)
            count += 1
    if count == 0:
        raise ValueError(f"no records found in FASTQ pair: {read1}, {read2}")
    return count


def make_probe_fastq_list(config: dict, sample: dict, scratch: Path, limit: int) -> Path:
    probe = scratch / "license_probe" / sample["sample_id"]
    probe.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, read1_key, read2_key in (
        ("GEX_1", "gex_r1", "gex_r2"), ("CRISPR_1", "crispr_r1", "crispr_r2")
    ):
        r1, r2 = probe / f"{label}_R1.fastq.gz", probe / f"{label}_R2.fastq.gz"
        observed = make_fastq_subset(Path(sample[read1_key]), Path(sample[read2_key]), r1, r2, limit)
        rows.append((label, sample["sample_id"], sample["rglb"], sample["lane"], str(r1), str(r2)))
        print(f"PILOT_FASTQ\t{label}\treads={observed}")
    fastq_list = probe / "fastq_list.csv"
    with fastq_list.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("RGID", "RGSM", "RGLB", "Lane", "Read1File", "Read2File"))
        writer.writerows(rows)
    return fastq_list
