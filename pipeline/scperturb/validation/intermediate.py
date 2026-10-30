#!/usr/bin/env python3
"""Validate an offline bioinformatics HTML report and optional delivery manifest."""

from __future__ import annotations

import argparse
import base64
import binascii
import csv
import json
import re
from collections import Counter
from pathlib import Path


ATTR_RE = r'''\b{key}\s*=\s*["']([^"']*)["']'''


def add_error(result: dict, message: str) -> None:
    result["errors"].append(message)


def attribute(attrs: str, key: str) -> str:
    match = re.search(ATTR_RE.format(key=re.escape(key)), attrs, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def has_class(attrs: str, class_name: str) -> bool:
    return class_name in attribute(attrs, "class").split()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_data_uris(text: str, result: dict) -> None:
    uris = re.findall(r'''\b(?:src|href)=["'](data:[^"']+)["']''', text, flags=re.IGNORECASE)
    result["data_uri_count"] = len(uris)
    for index, uri in enumerate(uris, start=1):
        if ";base64," not in uri:
            add_error(result, f"data URI {index} is not base64 encoded")
            continue
        payload = uri.split(";base64,", 1)[1]
        try:
            decoded = base64.b64decode(payload, validate=True)
        except (ValueError, binascii.Error):
            add_error(result, f"data URI {index} has invalid base64")
            continue
        if not decoded:
            add_error(result, f"data URI {index} is empty")


def validate_offline_resources(text: str, result: dict) -> None:
    hrefs = re.findall(r'''\bhref=["']([^"']+)["']''', text, flags=re.IGNORECASE)
    srcs = re.findall(r'''\bsrc=["']([^"']+)["']''', text, flags=re.IGNORECASE)
    srcsets = re.findall(r'''\bsrcset=["']([^"']+)["']''', text, flags=re.IGNORECASE)
    css_urls = [value.strip().strip('"\'') for value in re.findall(r"url\(([^)]+)\)", text, flags=re.IGNORECASE)]

    external_hrefs = [value for value in hrefs if not value.startswith(("#", "data:"))]
    external_srcs = [value for value in srcs if not value.startswith("data:")]
    external_srcsets = [
        value for value in srcsets
        if not value.lstrip().startswith("data:") or re.search(r"(?:https?://|file:)", value, flags=re.IGNORECASE)
    ]
    external_css_urls = [value for value in css_urls if value and not value.startswith(("#", "data:"))]

    result["external_hrefs"] = len(external_hrefs)
    result["external_srcs"] = len(external_srcs)
    result["external_srcsets"] = len(external_srcsets)
    result["external_css_urls"] = len(external_css_urls)
    if external_hrefs:
        add_error(result, "report contains external, local-file, or relative href values")
    if external_srcs:
        add_error(result, "report contains non-data src values")
    if external_srcsets:
        add_error(result, "report contains non-data or external srcset values")
    if external_css_urls:
        add_error(result, "report contains external, local-file, or relative CSS url() values")


def validate_structure(text: str, result: dict) -> None:
    ids = re.findall(r'''\bid=["']([^"']+)["']''', text, flags=re.IGNORECASE)
    duplicate_ids = sorted(value for value, count in Counter(ids).items() if count > 1)
    result["duplicate_id_count"] = len(duplicate_ids)
    if duplicate_ids:
        add_error(result, f"report contains duplicate IDs: {', '.join(duplicate_ids[:5])}")

    hrefs = re.findall(r'''\bhref=["'](#[^"']+)["']''', text, flags=re.IGNORECASE)
    known_ids = set(ids)
    stale = [href for href in hrefs if href[1:] not in known_ids]
    result["stale_anchor_count"] = len(stale)
    if stale:
        add_error(result, f"report contains stale anchors: {', '.join(stale[:5])}")

    nav_body = ""
    for attrs, body in re.findall(r"<nav\b([^>]*)>(.*?)</nav>", text, flags=re.IGNORECASE | re.DOTALL):
        if has_class(attrs, "sidebar"):
            nav_body = body
            break
    if not nav_body:
        add_error(result, "report sidebar is missing")

    section_ids = []
    for attrs in re.findall(r"<section\b([^>]*)>", text, flags=re.IGNORECASE):
        if has_class(attrs, "report-section"):
            section_id = attribute(attrs, "id")
            if section_id:
                section_ids.append(section_id)
    result["section_count"] = len(section_ids)
    nav_ids = [href[1:] for href in re.findall(r'''\bhref=["'](#[^"']+)["']''', nav_body, flags=re.IGNORECASE)]
    result["sidebar_anchor_count"] = len(nav_ids)
    if section_ids and nav_body and nav_ids != section_ids:
        add_error(result, "sidebar anchors must match report sections exactly and in order")

    if not re.search(r'''class=["'][^"']*\bfooter\b''', text, flags=re.IGNORECASE):
        add_error(result, "report footer is missing")


def table_body_rows(body: str) -> int:
    bodies = re.findall(r"<tbody\b[^>]*>(.*?)</tbody>", body, flags=re.IGNORECASE | re.DOTALL)
    if bodies:
        return sum(len(re.findall(r"<tr\b", tbody, flags=re.IGNORECASE)) for tbody in bodies)
    total = len(re.findall(r"<tr\b", body, flags=re.IGNORECASE))
    return max(0, total - int(bool(re.search(r"<thead\b", body, flags=re.IGNORECASE))))


def validate_tables(text: str, max_rows: int, result: dict) -> None:
    tables = re.findall(r"<table\b([^>]*)>(.*?)</table>", text, flags=re.IGNORECASE | re.DOTALL)
    result["table_count"] = len(tables)
    paginated_count = 0
    for index, (attrs, body) in enumerate(tables, start=1):
        rows = table_body_rows(body)
        paginated = "data-page-size" in attrs.lower() or has_class(attrs, "paginated-table") or has_class(attrs, "pagination-table")
        if not paginated:
            if rows > max_rows:
                add_error(result, f"table {index} has {rows} rows without pagination")
            continue

        paginated_count += 1
        table_id = attribute(attrs, "id")
        if not table_id:
            add_error(result, f"paginated table {index} has no ID")
            continue
        pager_match = re.search(
            rf'''<div\b(?=[^>]*\bdata-pager-for=["']{re.escape(table_id)}["'])[^>]*>(.*?)</div>''',
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
        if not pager_match:
            add_error(result, f"paginated table {table_id} has no matching pager")
            continue
        pager = pager_match.group(1)
        for control in ("data-page-prev", "data-page-next", "data-page-status"):
            if control not in pager:
                add_error(result, f"paginated table {table_id} pager is missing {control}")
    result["paginated_table_count"] = paginated_count


def figure_cards(text: str) -> list[tuple[str, str]]:
    cards = []
    for tag in ("article", "div"):
        pattern = rf"<{tag}\b([^>]*)>(.*?)</{tag}>"
        for attrs, body in re.findall(pattern, text, flags=re.IGNORECASE | re.DOTALL):
            if has_class(attrs, "figure-card"):
                cards.append((attrs, body))
    return cards


def validate_expandable_figures(text: str, required: bool, result: dict) -> None:
    cards = figure_cards(text)
    figure_cards_with_content = []
    non_expandable = []
    for index, (attrs, body) in enumerate(cards, start=1):
        if not re.search(r"<(?:img|svg)\b", body, flags=re.IGNORECASE):
            continue
        figure_cards_with_content.append(index)
        img_attrs = " ".join(re.findall(r"<img\b([^>]*)>", body, flags=re.IGNORECASE))
        expandable = (
            has_class(attrs, "expandable")
            or attribute(attrs, "data-expandable").lower() in {"1", "true", "yes"}
            or has_class(img_attrs, "expandable")
            or bool(re.search(r"\bonclick\s*=", img_attrs, flags=re.IGNORECASE))
            or bool(re.search(r'''<button\b[^>]*class=["'][^"']*(?:image-button|figure-expand|expandable)[^"']*["']''', body, flags=re.IGNORECASE))
        )
        if not expandable:
            non_expandable.append(index)
    result["figure_card_count"] = len(figure_cards_with_content)
    result["non_expandable_figure_count"] = len(non_expandable)
    if required and not figure_cards_with_content:
        add_error(result, "report has no figure cards with expandable content")
    elif required and non_expandable:
        add_error(result, f"figure cards are not expandable: {', '.join(map(str, non_expandable[:5]))}")


def resolve_manifest_target(value: str, manifest: Path) -> Path | None:
    if not value:
        return None
    target = Path(value)
    if not target.is_absolute():
        target = manifest.parent / target
    return target.resolve()


def validate_manifest(report: Path, manifest: Path, result: dict) -> None:
    if not manifest.is_file() or manifest.stat().st_size == 0:
        add_error(result, f"manifest missing or empty: {manifest}")
        return
    rows = read_manifest(manifest)
    result["manifest_rows"] = len(rows)
    resolved_targets = []
    for row in rows:
        target = resolve_manifest_target(row.get("target", ""), manifest)
        if target is not None:
            resolved_targets.append(target)
        if row.get("status") not in {"copied", "up_to_date"}:
            continue
        if target is None or not target.is_file() or target.stat().st_size == 0:
            add_error(result, f"manifest target missing for {row.get('description', 'unnamed item')}: {row.get('target', '')}")

    duplicates = sorted(str(target) for target, count in Counter(resolved_targets).items() if count > 1)
    result["duplicate_manifest_target_count"] = len(duplicates)
    if duplicates:
        add_error(result, f"manifest contains duplicate targets: {', '.join(duplicates[:3])}")

    report_indexes = [
        index for index, row in enumerate(rows)
        if resolve_manifest_target(row.get("target", ""), manifest) == report
    ]
    if len(report_indexes) != 1:
        add_error(result, f"manifest must contain exactly one final-report target; found {len(report_indexes)}")
    elif report_indexes[0] != len(rows) - 1:
        add_error(result, "final-report manifest row must be last")
    result["final_report_last"] = len(report_indexes) == 1 and report_indexes[0] == len(rows) - 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--require-section", action="append", default=[])
    parser.add_argument("--forbid-text", action="append", default=[])
    parser.add_argument("--max-unpaginated-rows", type=int, default=20)
    parser.add_argument("--require-expandable-figures", action="store_true")
    parser.add_argument("--require-print-css", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()

    report = args.report.resolve()
    result: dict[str, object] = {"report": str(report), "errors": []}
    text = ""
    if not report.is_file() or report.stat().st_size == 0:
        add_error(result, f"report missing or empty: {report}")
    else:
        text = report.read_text(encoding="utf-8", errors="replace")
        validate_offline_resources(text, result)
        validate_structure(text, result)
        validate_data_uris(text, result)
        validate_tables(text, args.max_unpaginated_rows, result)
        validate_expandable_figures(text, args.require_expandable_figures, result)
        result["has_print_css"] = bool(re.search(r"@media\s+print\b", text, flags=re.IGNORECASE))
        if args.require_print_css and not result["has_print_css"]:
            add_error(result, "report print CSS is missing")
        for section in args.require_section:
            if section not in text:
                add_error(result, f"required section is missing: {section}")
        for term in args.forbid_text:
            if term and term in text:
                add_error(result, f"forbidden text found in report: {term}")

    if args.manifest:
        manifest = args.manifest.resolve()
        validate_manifest(report, manifest, result)
        if manifest.is_file():
            manifest_text = manifest.read_text(encoding="utf-8", errors="replace")
            for term in args.forbid_text:
                if term and term in manifest_text:
                    add_error(result, f"forbidden text found in manifest: {term}")

    result["ok"] = not result["errors"]
    if args.as_json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(f"report: {result['report']}")
        print(f"ok: {result['ok']}")
        for error in result["errors"]:
            print(f"ERROR: {error}")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
