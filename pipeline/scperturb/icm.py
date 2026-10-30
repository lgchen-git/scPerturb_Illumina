"""Local handoff and returned-export contract for ICM Perturb-seq."""

from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path


HANDOFF_DIRNAME = "icm_perturbseq_handoff"
MANIFEST_NAME = "icm_export_manifest.tsv"
INPUT_SUFFIXES = (
    ".scRNA.filtered.matrix.mtx.gz",
    ".scRNA.filtered.barcodes.tsv.gz",
    ".scRNA.filtered.features.tsv.gz",
    ".scRNA.feature_barcode_reference.csv",
    ".scRNA.positive_cell_guide_assignments.csv",
)
RETURN_ASSETS = {
    "expression_overview": {
        "filename": "expression_overview.png",
        "title": "ICM expression and cell composition overview",
        "caption": "Static plot exported from the Illumina Connected Multiomics Perturb-seq workflow.",
        "report_file": "icm_expression_overview.png",
    },
    "guide_capture_overview": {
        "filename": "guide_capture_overview.png",
        "title": "ICM CRISPR guide capture overview",
        "caption": "Static guide-capture plot exported from the Illumina Connected Multiomics Perturb-seq workflow.",
        "report_file": "icm_guide_capture_overview.png",
    },
    "analysis_summary": {
        "filename": "analysis_summary.pdf",
        "title": "ICM Perturb-seq analysis summary",
        "caption": "Summary exported from the Illumina Connected Multiomics Perturb-seq workflow.",
    },
}

HANDOFF_README = """# Illumina Connected Multiomics Perturb-seq handoff

This folder contains the five DRAGEN outputs accepted by the ICM Perturb-seq
workflow. Upload these files manually in your ICM account and select
**Default: Illumina Single Cell Transcriptomics Perturb-seq**.

1. Upload every file under `input/` without renaming or editing it.
2. Complete the official ICM workflow and review its outputs.
3. Export the expression/cell-composition overview and guide-capture overview
   as PNG, and export the analysis summary as PDF.
4. Place them in `returned/` using the filenames in `icm_export_manifest.tsv`.
5. Run the report refresh only after all three files are present.

The package builder only copies local files. It does not connect to ICM, upload
data, or store credentials. Uploaded data is handled by the ICM account and
service.
"""

RETURN_README = """# ICM return files

Copy the completed ICM exports into this directory using the filenames listed
in `icm_export_manifest.tsv`. The report refresh validates these files and
embeds the two PNGs in the standalone offline HTML. The PDF is linked locally
from the report's reproducibility section.
"""


def handoff_root(stage2: Path) -> Path:
    return stage2 / "04_tool_output" / HANDOFF_DIRNAME


def validate_stage2(stage2: Path, sample_id: str) -> Path:
    marker = stage2 / "05_metadata/stage_complete.json"
    if not marker.is_file():
        raise FileNotFoundError(f"ICM handoff requires completed DRAGEN stage: {marker}")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    if payload.get("status") != "complete" or payload.get("sample_id") != sample_id:
        raise ValueError(f"invalid DRAGEN completion marker for {sample_id}: {marker}")
    dragen_dir = stage2 / "04_tool_output/dragen"
    missing = [dragen_dir / f"{sample_id}{suffix}" for suffix in INPUT_SUFFIXES
               if not (dragen_dir / f"{sample_id}{suffix}").is_file()
               or (dragen_dir / f"{sample_id}{suffix}").stat().st_size == 0]
    if missing:
        raise FileNotFoundError("ICM handoff missing required DRAGEN output(s): " + ", ".join(map(str, missing)))
    return dragen_dir


def prepare_handoff(stage2: Path, sample_id: str, *, dry_run: bool = False) -> Path:
    dragen_dir = validate_stage2(stage2, sample_id)
    target = handoff_root(stage2)
    if dry_run:
        return target
    input_dir = target / "input"
    returned = target / "returned"
    input_dir.mkdir(parents=True, exist_ok=True)
    returned.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for suffix in INPUT_SUFFIXES:
        source = dragen_dir / f"{sample_id}{suffix}"
        destination = input_dir / source.name
        temporary = destination.with_name(destination.name + ".partial")
        try:
            shutil.copyfile(source, temporary)
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
        manifest_rows.append((destination.name, source.relative_to(stage2).as_posix(), destination.stat().st_size))
    with (target / "input_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("file", "source_relative_to_dragen_stage", "size_bytes"))
        writer.writerows(manifest_rows)
    (target / "README.md").write_text(HANDOFF_README, encoding="utf-8")
    (returned / "README.md").write_text(RETURN_README, encoding="utf-8")
    manifest = returned / MANIFEST_NAME
    if not manifest.exists():
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("role", "filename", "title", "caption"))
            writer.writerows((role, asset["filename"], asset["title"], asset["caption"])
                             for role, asset in RETURN_ASSETS.items())
    return target


def validate_returned(handoff: Path) -> dict[str, dict[str, str]]:
    returned = handoff / "returned"
    manifest = returned / MANIFEST_NAME
    if not manifest.is_file():
        raise FileNotFoundError(f"ICM exports not returned; add {manifest}")
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"role", "filename", "title", "caption"}.issubset(reader.fieldnames or []):
            raise ValueError(f"invalid ICM export manifest columns: {manifest}")
        rows = list(reader)
    by_role = {row["role"]: row for row in rows}
    if len(by_role) != len(rows) or set(by_role) != set(RETURN_ASSETS):
        raise ValueError(f"ICM export manifest must contain exactly these roles: {', '.join(RETURN_ASSETS)}")
    for role, asset in RETURN_ASSETS.items():
        row = by_role[role]
        if row.get("filename") != asset["filename"]:
            raise ValueError(f"ICM manifest filename for {role} must be {asset['filename']}")
        if not row.get("title", "").strip() or not row.get("caption", "").strip():
            raise ValueError(f"ICM export manifest title/caption is empty for {role}")
        path = returned / asset["filename"]
        if not path.is_file() or path.stat().st_size < 8:
            raise FileNotFoundError(f"ICM export is missing or empty: {path}")
        with path.open("rb") as handle:
            signature = handle.read(8)
        expected = b"\x89PNG\r\n\x1a\n" if path.suffix == ".png" else b"%PDF-"
        if not signature.startswith(expected):
            raise ValueError(f"ICM export has an invalid {path.suffix} signature: {path}")
    return by_role


def install_returned_assets(handoff: Path, delivery: Path) -> dict[str, dict[str, str]]:
    """Copy validated ICM exports into the offline delivery and update its figure manifest."""
    rows = validate_returned(handoff)
    figures = delivery / "02_Figures/figures"
    icm_dir = delivery / "04_Methods_Reproducibility/icm"
    figures.mkdir(parents=True, exist_ok=True)
    icm_dir.mkdir(parents=True, exist_ok=True)
    returned = handoff / "returned"
    manifest_rows = []
    for role, asset in RETURN_ASSETS.items():
        row = rows[role]
        destination = icm_dir / asset["filename"] if role == "analysis_summary" else figures / asset["report_file"]
        shutil.copyfile(returned / asset["filename"], destination)
        if role != "analysis_summary":
            manifest_rows.append({
                "file": asset["report_file"], "title": row["title"], "section": "icm_perturbseq",
                "caption": row["caption"], "source_analysis_version": "Illumina Connected Multiomics Perturb-seq",
                "report_revision": "v3.2.1",
            })
    shutil.copy2(handoff / "returned" / MANIFEST_NAME, icm_dir / MANIFEST_NAME)
    (icm_dir / "README.md").write_text(
        "# ICM Perturb-seq exports\n\n"
        "The PNGs and summary PDF in this folder were returned by the operator from the official ICM workflow. "
        "The PNGs are embedded in the offline report; no ICM dashboard or remote asset is loaded.\n",
        encoding="utf-8",
    )
    figure_manifest = delivery / "03_Tables/figure_manifest.tsv"
    with figure_manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        current = [row for row in reader if row.get("file") not in {x["file"] for x in manifest_rows}]
    if not fields:
        raise ValueError(f"figure manifest has no header: {figure_manifest}")
    with figure_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(current)
        writer.writerows({key: value for key, value in row.items() if key in fields} for row in manifest_rows)
    return rows
