#!/usr/bin/env python3
"""Generate an official DRAGEN sample report and verify shared sample metrics."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path


class ResourceParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.external_resources: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        resource = values.get("src") if tag in {"img", "script", "iframe", "source"} else None
        if tag == "link" and "stylesheet" in str(values.get("rel", "")).lower():
            resource = values.get("href")
        if resource and re.match(r"https?://", str(resource), flags=re.IGNORECASE):
            self.external_resources.append(str(resource))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_row_data(report: Path) -> dict[str, list[str]]:
    source = report.read_text(encoding="utf-8", errors="replace")
    decoder = json.JSONDecoder()
    values: dict[str, list[str]] = {}
    for match in re.finditer(r"rowData\s*:\s*\[", source):
        bracket = source.find("[", match.start())
        try:
            rows, _ = decoder.raw_decode(source[bracket:])
        except (json.JSONDecodeError, ValueError):
            continue
        if not isinstance(rows, list):
            continue
        for row in rows:
            if isinstance(row, dict) and row.get("metric") is not None and row.get("value") is not None:
                values.setdefault(str(row["metric"]), []).append(str(row["value"]))
    if not values:
        raise ValueError(f"official report has no readable metric rows: {report}")
    return values


def unique_official_value(values: dict[str, list[str]], metric: str) -> str:
    candidates = values.get(metric, [])
    if not candidates:
        raise ValueError(f"official report is missing overlapping metric: {metric}")
    parsed = {parse_number(value) for value in candidates}
    if len(parsed) != 1:
        raise ValueError(f"official report has ambiguous values for {metric}: {candidates}")
    return candidates[0]


def parse_number(value: str) -> float:
    text = value.strip().replace(",", "").replace("%", "")
    match = re.fullmatch(r"([-+]?(?:\d+(?:\.\d*)?|\.\d+))(?:\s*[KkMm])?", text)
    if not match:
        raise ValueError(f"not a numeric report value: {value!r}")
    number = float(match.group(1))
    if text.rstrip().endswith(("M", "m")):
        number *= 1_000_000
    elif text.rstrip().endswith(("K", "k")):
        number *= 1_000
    return number


def overlap_rows(args: argparse.Namespace, official: dict[str, list[str]]) -> list[dict[str, str]]:
    if args.mode == "scrna":
        qc_rows = read_tsv(args.qc_table)
        own = next((row for row in qc_rows if row.get("library_id") == args.sample_id), None)
        if own is None and len(qc_rows) == 1:
            own = qc_rows[0]
        if own is None:
            raise ValueError(f"no unique QC row for sample {args.sample_id}: {args.qc_table}")
        fields = [
            ("Passing Cells", "passing_cells", "cells", 0),
            ("Mean Reads per Cell", "mean_reads_per_passing_cell", "reads/cell", 1),
            ("% Sequencing Saturation", "sequencing_saturation_percent", "percent", 1),
            ("Median Molecules per Cell", "median_molecules_per_passing_cell", "molecules/cell", 0),
            ("Median Genes per Cell", "median_genes_per_passing_cell", "genes/cell", 0),
        ]
        source_path = args.qc_table
    else:
        qc = json.loads(args.guide_qc_json.read_text(encoding="utf-8"))
        assignment = qc.get("assignment_status", {})
        passing = int(qc["filtered_cells"])
        singlet, multi, unassigned = (int(assignment[key]) for key in ("singlet", "multi", "none"))
        positive = int(qc["positive_assignment_cells"])
        if singlet + multi != positive or singlet + multi + unassigned != passing:
            raise ValueError("cell-guide QC counts do not form a consistent partition")
        own = {
            "filtered_cells": str(passing),
            "positive_assignment_cells": str(positive),
            "multi_guide_cells": str(multi),
        }
        fields = [
            ("Passing Cells", "filtered_cells", "cells", 0),
            ("Cells with >= 1 Guide", "positive_assignment_cells", "cells", 0),
            ("Cells with >= 2 Guides", "multi_guide_cells", "cells", 0),
        ]
        source_path = args.guide_qc_json

    rows: list[dict[str, str]] = []
    failures: list[str] = []
    for official_name, source_key, unit, precision in fields:
        report_value = unique_official_value(official, official_name)
        own_value = own[source_key]
        expected = parse_number(own_value)
        observed = parse_number(report_value)
        if precision == 0:
            matches = expected.is_integer() and observed.is_integer() and int(expected) == int(observed)
            tolerance = "exact integer"
        else:
            matches = round(expected, precision) == round(observed, precision)
            tolerance = f"rounded to {precision} decimal place(s)"
        status = "PASS" if matches else "FAIL"
        rows.append({
            "sample_id": args.sample_id,
            "metric": official_name,
            "official_value": report_value,
            "workflow_value": own_value,
            "unit": unit,
            "comparison": tolerance,
            "source": str(source_path),
            "status": status,
        })
        if not matches:
            failures.append(f"{official_name}: official={report_value}; workflow={own_value}")
    write_tsv_atomic(args.parity, rows)
    if failures:
        raise ValueError("DRAGEN report parity failed: " + "; ".join(failures))
    return rows


def write_tsv_atomic(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".building")
    columns = ["sample_id", "metric", "official_value", "workflow_value", "unit", "comparison", "source", "status"]
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def validate_offline_report(path: Path, sample_id: str) -> None:
    if not path.is_file() or path.stat().st_size < 100_000:
        raise ValueError(f"official report is missing or unexpectedly small: {path}")
    source = path.read_text(encoding="utf-8", errors="replace")
    if "<html" not in source.lower() or sample_id not in source:
        raise ValueError(f"official report does not look like a report for {sample_id}: {path}")
    parser = ResourceParser()
    parser.feed(source)
    if parser.external_resources:
        raise ValueError("official report has external resources: " + ", ".join(parser.external_resources[:5]))


def promote_report(source: Path, destination: Path) -> None:
    """Copy across filesystems, then atomically replace on the destination filesystem."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=".building", dir=destination.parent, delete=False,
    ) as handle:
        staged = Path(handle.name)
    try:
        shutil.copy2(source, staged)
        staged.replace(destination)
    finally:
        staged.unlink(missing_ok=True)


def build_report(args: argparse.Namespace) -> list[dict[str, str]]:
    binary = Path(args.binary)
    report_root = Path(args.reports_root)
    manifest = Path(args.manifest)
    source_dir = args.dragen_output_dir.resolve()
    for path in (binary, manifest, source_dir):
        if not path.exists():
            raise FileNotFoundError(path)
    if not report_root.is_dir():
        raise FileNotFoundError(report_root)
    version = subprocess.run([str(binary), "--version"], check=True, text=True, capture_output=True).stdout.strip()
    if "4.5.4" not in version:
        raise RuntimeError(f"expected DRAGEN Reports 4.5.4, found {version!r}")
    temp_root = args.temp_root.resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"dragen-reports-{args.sample_id}-", dir=temp_root) as work:
        task = Path(work)
        input_root = task / "input"
        input_root.mkdir()
        (input_root / args.sample_id).symlink_to(source_dir, target_is_directory=True)
        temporary_report = task / f"{args.sample_id}_DRAGEN_report.html"
        command = [
            str(binary), "-f", "-m", str(manifest), "-d", str(input_root),
            "-o", str(task / "workflow.html"), "-s", args.sample_id,
            "-S", str(temporary_report), "-r", str(report_root),
        ]
        subprocess.run(command, check=True, cwd=task, text=True)
        validate_offline_report(temporary_report, args.sample_id)
        parity = []
        if not args.generate_only:
            official = read_row_data(temporary_report)
            parity = overlap_rows(args, official)
        promote_report(temporary_report, args.report)
    alignment_status = "parity=deferred_to_workflow" if args.generate_only else f"shared_metrics={len(parity)}"
    print(f"PASS\t{args.sample_id}\tDRAGEN Reports {version}\treport={args.report}\t{alignment_status}")
    return parity


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("scrna", "scperturb"), required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--dragen-output-dir", type=Path, required=True)
    parser.add_argument("--binary", default="/usr/bin/dragen-reports")
    parser.add_argument("--reports-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--parity", type=Path)
    parser.add_argument("--generate-only", action="store_true", help="Generate official HTML; workflow performs scoped parity separately")
    parser.add_argument("--temp-root", type=Path, required=True)
    parser.add_argument("--qc-table", type=Path)
    parser.add_argument("--guide-qc-json", type=Path)
    args = parser.parse_args()
    if not args.generate_only and not args.parity:
        parser.error("--parity is required unless --generate-only is selected")
    if not args.generate_only and args.mode == "scrna" and not args.qc_table:
        parser.error("--qc-table is required for scRNA mode")
    if not args.generate_only and args.mode == "scperturb" and not args.guide_qc_json:
        parser.error("--guide-qc-json is required for scPerturb mode")
    try:
        build_report(args)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
