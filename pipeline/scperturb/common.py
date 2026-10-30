#!/usr/bin/env python3
"""Shared helpers for the native scPerturb Snakemake stages."""

from __future__ import annotations

import base64
import csv
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


STAGE_SUBDIRS = ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata")


def load_json(path: str | Path) -> dict:
    resolved = Path(path).resolve()
    with resolved.open(encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"expected JSON object: {resolved}")
    data["_config_path"] = str(resolved)
    data["_config_sha256"] = hashlib.sha256(resolved.read_bytes()).hexdigest()
    return data


def load_analysis(config: dict) -> dict:
    data = load_json(config["analysis_config"])
    return data


def sample_rows(config: dict) -> list[dict[str, str]]:
    with Path(config["sample_sheet"]).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("sample sheet is empty")
    required = {"sample_id", "project_id", "lane", "rglb", "gex_r1", "gex_r2", "crispr_r1", "crispr_r2"}
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"sample sheet missing columns: {','.join(missing)}")
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("sample sheet contains duplicate sample_id")
    return rows


def sample_row(config: dict, sample_id: str) -> dict[str, str]:
    for row in sample_rows(config):
        if row["sample_id"] == sample_id:
            return row
    raise KeyError(f"unknown sample: {sample_id}")


def project_path(config: dict, path: str | Path) -> Path:
    """Resolve a configured project-relative path against the container project root."""
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path(config["project_root"]) / candidate


def delivery_name(config: dict, sample_id: str) -> str:
    """Return the configurable delivery basename; old runs default to delivery."""
    name = config.get("delivery_name", "delivery").format(sample_id=sample_id)
    if not name or name in {".", ".."} or Path(name).name != name:
        raise ValueError("delivery_name must be a single folder name")
    return name


def stage_root(config: dict, sample_id: str, stage_dir: str) -> Path:
    return Path(config["results_root"]) / sample_id / stage_dir


def ensure_stage(root: Path) -> None:
    for name in STAGE_SUBDIRS:
        (root / name).mkdir(parents=True, exist_ok=True)


def ensure_runtime_dirs(config: dict, sample_id: str) -> None:
    root = Path(config["tmp_root"]) / sample_id
    for name in ("runtime", "cache", "matplotlib", "numba", "dragen_intermediate"):
        (root / name).mkdir(parents=True, exist_ok=True)
    Path(config["log_root"]).mkdir(parents=True, exist_ok=True)
    Path(config["results_root"]).mkdir(parents=True, exist_ok=True)


def require_file(path: str | Path, label: str, *, nonempty: bool = True) -> Path:
    resolved = Path(path).resolve()
    if not resolved.is_file() or (nonempty and resolved.stat().st_size == 0):
        raise FileNotFoundError(f"missing {label}: {resolved}")
    return resolved


def file_digest(path: str | Path, algorithm: str = "sha256", chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: str | Path, header: Iterable[object], rows: Iterable[Iterable[object]]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(list(header))
        writer.writerows(rows)
    return target


def write_json(path: str | Path, payload: dict) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    building = target.with_suffix(target.suffix + ".building")
    building.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    building.replace(target)
    return target


def output_manifest(stage: Path, outputs: Iterable[str | Path]) -> Path:
    rows = []
    for raw in sorted({Path(value).resolve() for value in outputs}, key=str):
        if not raw.is_file():
            raise FileNotFoundError(f"manifest output missing: {raw}")
        size = raw.stat().st_size
        digest = file_digest(raw) if size <= 64 * 1024 * 1024 else "size_verified_only"
        rows.append((str(raw), size, digest))
    target = stage / "05_metadata" / "output_manifest.tsv"
    return write_tsv(target, ("path", "size_bytes", "sha256_or_size_check"), rows)


def completion_marker(stage: Path) -> Path:
    return stage / "05_metadata" / "stage_complete.json"


def write_completion(stage: Path, config: dict, sample_id: str, outputs: Iterable[str | Path], *, interface: str) -> Path:
    marker = completion_marker(stage)
    payload = {
        "schema_version": 1,
        "status": "complete",
        "sample_id": sample_id,
        "stage": stage.name,
        "interface": interface,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "config_sha256": config.get("_config_sha256", ""),
        "outputs": [str(Path(value).resolve()) for value in outputs],
    }
    return write_json(marker, payload)


def git_provenance(project_root: str | Path, stage: Path) -> Path:
    root = Path(project_root)
    head = os.environ.get("SC_PERTURB_GIT_HEAD")
    status_b64 = os.environ.get("SC_PERTURB_GIT_STATUS_BASE64")
    if head is not None and status_b64 is not None:
        dirty = bool(status_b64)
        rows = [
            ("head", 0, head),
            ("status", 0, "dirty" if dirty else "clean"),
            ("status_porcelain_base64", 0, status_b64),
        ]
    else:
        rows = []
        try:
            for args, label in ((["rev-parse", "HEAD"], "head"), (["status", "--porcelain=v1"], "status")):
                result = subprocess.run(["git", "-C", str(root), *args], check=False, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                rows.append((label, result.returncode, result.stdout.strip()))
        except FileNotFoundError:
            rows = [("head", 127, "Git CLI unavailable; run through pipeline/docker/run_snakemake.sh"),
                    ("status", 127, "Git CLI unavailable; host provenance was not supplied")]
    return write_tsv(stage / "05_metadata" / "git_provenance.tsv", ("field", "returncode", "value"), rows)


def set_thread_env(config: dict, sample_id: str) -> dict[str, str]:
    ensure_runtime_dirs(config, sample_id)
    root = Path(config["tmp_root"]) / sample_id
    env = os.environ.copy()
    env.update({
        "TMPDIR": str(root / "runtime"),
        "TMP": str(root / "runtime"),
        "TEMP": str(root / "runtime"),
        "XDG_CACHE_HOME": str(root / "cache"),
        "MPLCONFIGDIR": str(root / "matplotlib"),
        "NUMBA_CACHE_DIR": str(root / "numba"),
        "OMP_NUM_THREADS": str(config["resources"]["analysis_threads"]),
        "OPENBLAS_NUM_THREADS": str(config["resources"]["analysis_threads"]),
        "MKL_NUM_THREADS": str(config["resources"]["analysis_threads"]),
        "NUMBA_NUM_THREADS": str(config["resources"]["analysis_threads"]),
    })
    return env
