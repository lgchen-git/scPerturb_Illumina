#!/usr/bin/env python3
"""Validate rendered conserved-flow module artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path


class ReportParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.ids: list[str] = []
        self.tags_by_id: dict[str, str] = {}
        self.classes: set[str] = set()
        self.external: list[tuple[str, str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for key, value in attrs:
            if key == "id" and value:
                self.ids.append(value)
                self.tags_by_id[value] = tag
            if key == "class" and value:
                self.classes.update(value.split())
            if key in {"src", "href", "poster"} and value and not value.startswith(("#", "data:")):
                self.external.append((tag, key, value))


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_html(path: Path) -> tuple[str, ReportParser]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty HTML: {path}")
    text = path.read_text(encoding="utf-8", errors="strict")
    parser = ReportParser()
    parser.feed(text)
    duplicates = sorted({value for value in parser.ids if parser.ids.count(value) > 1})
    if duplicates:
        raise RuntimeError(f"Duplicate HTML IDs in {path}: {duplicates}")
    if parser.external:
        raise RuntimeError(f"External resources in {path}: {parser.external}")
    return text, parser


def main() -> int:
    cli = arguments()
    try:
        root = cli.output_dir
        nodes = tsv_rows(root / "tables/read_flow_nodes.tsv")
        links = tsv_rows(root / "tables/read_flow_links.tsv")
        payload_path = root / "tables/read_flow_payload.json"
        payload = json.loads(payload_path.read_text(encoding="utf-8"))
        if not payload.get("views"):
            raise RuntimeError("Payload has no views")
        node_layout = payload.get("node_layout", "snap")
        if node_layout not in {"snap", "main_top_losses"}:
            raise RuntimeError(f"Unsupported node layout in payload: {node_layout}")
        orientation = payload.get("orientation", "horizontal")
        if orientation not in {"horizontal", "vertical"}:
            raise RuntimeError(f"Unsupported orientation in payload: {orientation}")
        validation_path = root / "reports/sankey_validation.json"
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        if not validation.get("ok") or not validation.get("presentation_checked"):
            raise RuntimeError("Rendered output validation record is incomplete")
        report_path = root / "reports/read_flow_sankey.html"
        fragment_path = root / "reports/read_flow_sankey_section.html"
        report, report_parser = parse_html(report_path)
        fragment, fragment_parser = parse_html(fragment_path)
        for token in (
            "@media print", "bio-flow-view", "bio-flow-detail", "bio-flow-expand",
            "Plotly.react", "data.node_layout==='main_top_losses'?'fixed':'snap'",
            "nodeSpec.y=nodes.map(n=>n.y)", "raiseSankeyLabels(el)",
            "group.querySelector('.sankey-node-set')", "nodes.parentNode.appendChild(nodes)",
            "coordinateRoot.appendChild(layer)", "bio-flow-vertical-label-layer", "bio-flow-label-layer", "plotly_afterplot",
            "plotView(el,view,forcedWidth)", "el.innerHTML=view.vertical_svg",
            "Plotly.purge(el)", "Plotly.newPlot",
            "installVerticalInteractions(el,view)", "bio-flow-svg-tooltip",
            "max-height:78vh", "overflow-y:auto", "overscroll-behavior:contain",
            "Math.min(780,Math.floor(innerHeight*.72))", "max-height:none!important;overflow:visible!important",
            ".bio-flow-label-layer text{paint-order:stroke fill;stroke:#fff;stroke-width:4px",
        ):
            if token not in report:
                raise RuntimeError(f"Standalone HTML is missing {token}")
        for obsolete_id in ("bio-flow-connected", "bio-flow-tabs"):
            if obsolete_id in report_parser.ids or obsolete_id in fragment_parser.ids:
                raise RuntimeError(f"Obsolete multi-panel Sankey control remains: {obsolete_id}")
        obsolete_classes = {"bio-flow-tab", "bio-flow-tabs"}
        if obsolete_classes & (report_parser.classes | fragment_parser.classes):
            raise RuntimeError("Obsolete Sankey tabs remain in rendered HTML")
        if validation.get("node_layout") != node_layout:
            raise RuntimeError("Payload and validation node layouts differ")
        if validation.get("orientation") != orientation:
            raise RuntimeError("Payload and validation orientations differ")
        if node_layout == "main_top_losses" and "main_top_local_losses_layout" not in validation.get("checks", []):
            raise RuntimeError("Packed main/loss layout lacks validation record")
        if node_layout == "main_top_losses":
            for view_id, view in payload["views"].items():
                lookup = {node["node_id"]: node for node in view["nodes"]}
                for node in view["nodes"]:
                    if int(node["value"]) > 0 and not (0.0 < float(node["x"]) < 1.0 and 0.0 < float(node["y"]) < 1.0):
                        raise RuntimeError(f"{view_id}/{node['node_id']}: packed Plotly center is outside the visible domain")
                for stage in sorted({int(node["stage"]) for node in view["nodes"]}):
                    visible = [node for node in view["nodes"] if int(node["stage"]) == stage and int(node["value"]) > 0]
                    retained = [node for node in visible if node["node_kind"] == "retained"]
                    losses = sorted(
                        (node for node in visible if node["node_kind"] == "loss"),
                        key=lambda node: int(node["layout_order"]),
                    )
                    if retained and losses and max(float(node["y"]) for node in retained) > min(float(node["y"]) for node in losses):
                        raise RuntimeError(f"{view_id}/stage-{stage}: retained nodes are not above losses")
                    if [int(node["value"]) for node in losses] != sorted((int(node["value"]) for node in losses), reverse=True):
                        raise RuntimeError(f"{view_id}/stage-{stage}: loss nodes are not value-descending")
                for link in view["links"]:
                    target = lookup[link["target"]]
                    if target["node_kind"] != "loss" or int(target["value"]) == 0:
                        continue
                    source = lookup[link["source"]]
                    main_targets = [
                        node for node in view["nodes"]
                        if int(node["stage"]) == int(target["stage"]) and node["node_kind"] != "loss"
                    ]
                    if not float(source["x"]) < float(target["x"]):
                        raise RuntimeError(f"{view_id}/{target['node_id']}: loss target is not after its source")
                    if main_targets and not float(target["x"]) < min(float(node["x"]) for node in main_targets):
                        raise RuntimeError(f"{view_id}/{target['node_id']}: loss target is not in the local inter-stage lane")
        if "Plotly.react" not in fragment or len(fragment.encode("utf-8")) > 1_000_000:
            raise RuntimeError("Section fragment is missing interaction code or unexpectedly bundles Plotly")
        for token in (
            "raiseSankeyLabels(el)", "nodes.parentNode.appendChild(nodes)",
            "coordinateRoot.appendChild(layer)", "plotly_afterplot", "bio-flow-label-layer",
        ):
            if token not in fragment:
                raise RuntimeError(f"Section fragment is missing label-layer control: {token}")
        if "svg.main-svg" in report or "svg.main-svg" in fragment:
            raise RuntimeError("Interactive label overlays must stay inside the Sankey coordinate group")
        required_ids = {
            "bio-flow-view", "bio-flow-expand", "bio-flow-card", "bio-flow-detail",
            "bio-flow-plot", "bio-flow-table",
        }
        if not required_ids.issubset(report_parser.ids) or not required_ids.issubset(fragment_parser.ids):
            raise RuntimeError("Required Sankey controls are missing")
        if report_parser.tags_by_id.get("bio-flow-view") != "select" or fragment_parser.tags_by_id.get("bio-flow-view") != "select":
            raise RuntimeError("Sankey view control must be a dropdown select")
        if "bio-flow-orientation" in report_parser.ids or "bio-flow-orientation" in fragment_parser.ids:
            raise RuntimeError("Orientation must be fixed at build time, not exposed as a report control")
        svg_paths = sorted((root / "figures/read_flow_sankey").glob("*.svg"))
        if len(svg_paths) != len(payload["views"]):
            raise RuntimeError("Static SVG count does not match payload views")
        for path in svg_paths:
            svg = ET.parse(path).getroot()
            if svg.attrib.get("data-node-layout") != node_layout:
                raise RuntimeError(f"Static SVG node layout does not match payload: {path}")
            if svg.attrib.get("data-orientation") != orientation:
                raise RuntimeError(f"Static SVG orientation does not match payload: {path}")
            namespace = {"svg": "http://www.w3.org/2000/svg"}
            bands = [
                item for item in svg.findall(".//svg:path", namespace)
                if "data-link-id" in item.attrib
            ]
            blocks = [
                item for item in svg.findall(".//svg:rect", namespace)
                if "data-node-id" in item.attrib
            ]
            if not bands or not blocks:
                raise RuntimeError(f"Static SVG lacks value-scaled geometry: {path}")
            style = "".join(item.text or "" for item in svg.findall(".//svg:style", namespace))
            for token in ("paint-order:stroke fill", "stroke:#fff", "stroke-width:3px"):
                if token not in style:
                    raise RuntimeError(f"Static SVG label halo is missing {token}: {path}")
            children = list(svg)
            band_indexes = [
                index for index, item in enumerate(children)
                if item.tag.endswith("path") and "data-link-id" in item.attrib
            ]
            label_indexes = [
                index for index, item in enumerate(children)
                if item.tag.endswith("text") and item.attrib.get("class") == "label"
            ]
            if not label_indexes or max(band_indexes) >= min(label_indexes):
                raise RuntimeError(f"Static SVG labels are not layered above flow bands: {path}")
            ratios = [float(item.attrib["data-width"]) / int(item.attrib["data-value"]) for item in bands]
            ratios.extend(float(item.attrib["data-height"]) / int(item.attrib["data-value"]) for item in blocks)
            if max(ratios) - min(ratios) > max(ratios) * 1e-5:
                raise RuntimeError(f"Static SVG uses inconsistent value-to-pixel scales: {path}")
        transformations = payload.get("transformations", [])
        if transformations:
            transform_rows = tsv_rows(root / "tables/read_flow_transformations.tsv")
            if len(transform_rows) != len(transformations):
                raise RuntimeError("Transformation table does not match payload")
            if validation.get("default_mode") == "connected" and "cross_unit_transformations" not in validation.get("checks", []):
                raise RuntimeError("Connected layout lacks cross-unit validation")
        manifest_path = root / "reports/sankey_manifest.tsv"
        manifest = tsv_rows(manifest_path)
        for row in manifest:
            target = Path(row.get("target", ""))
            if not target.is_absolute():
                target = manifest_path.parent / target
            if row.get("status") not in {"up_to_date", "copied"}:
                raise RuntimeError(f"Invalid manifest status: {row}")
            if not target.is_file() or target.stat().st_size == 0:
                raise RuntimeError(f"Missing manifest target: {target}")
        summary = {
            "ok": True,
            "output_dir": str(root.resolve()),
            "views": len(payload["views"]),
            "nodes": len(nodes),
            "links": len(links),
            "transformations": len(transformations),
            "svg_files": len(svg_paths),
            "manifest_rows": len(manifest),
            "offline_resources": True,
            "unique_html_ids": True,
            "print_css": True,
            "node_layout": node_layout,
            "orientation": orientation,
            "compact_sankey_layout": node_layout == "main_top_losses",
            "value_scaled_svg": True,
            "dropdown_view_control": True,
            "build_time_orientation": True,
            "single_viewport_horizontal_layout": True,
            "scrollable_vertical_layout": True,
            "interactive_vertical_hover": True,
            "label_layer_checked": True,
        }
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0
    except (OSError, UnicodeError, ValueError, RuntimeError, ET.ParseError, json.JSONDecodeError) as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
