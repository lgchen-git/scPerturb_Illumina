#!/usr/bin/env python3
"""Validate a standalone V3 scPerturb delivery against its own source tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlparse
from scperturb.dragen.reports import validate_offline_report

SECTION_IDS = ["overview", "guide_assignment", "screen_usability", "coverage", "effects",
               "robustness", "cell_state", "pathways", "evidence", "limitations", "reproducibility"]
EVIDENCE_FILES = [f"fig{i:02d}_{name}.png" for i, name in (
    (1, "screen_usability_funnel"), (3, "target_coverage"),
    (4, "multi_guide_robustness"), (5, "representative_perturbation_effects"),
    (6, "pathway_ranking_qc"), (7, "reactome_ranking_sensitivity"),
)] + ["fig13_representative_de_effects.png", "fig14_dragen_barcode_rank.png"]

class ReportParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.sections, self.navigation, self.hrefs, self.sources, self.image_alts, self.text = [], [], [], [], [], []
        self.in_nav = False
        self.skip = 0
        self.table_id, self.in_body = None, False
        self.table_rows = {}
        self.table_headers, self.table_data = {}, {}
        self.in_head, self.current_row, self.current_header = False, None, None
        self.in_cell, self.current_cell = False, []
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in {"script", "style"}: self.skip += 1
        if tag == "nav": self.in_nav = True
        if tag == "section" and "report-section" in a.get("class", "").split(): self.sections.append(a.get("id"))
        if tag == "a" and a.get("href"):
            self.hrefs.append(a["href"])
            if self.in_nav and a["href"].startswith("#"): self.navigation.append(a["href"][1:])
        if tag in {"img", "script"}:
            key = "src"
            if a.get(key): self.sources.append(a[key])
        if tag == "img": self.image_alts.append(a.get("alt", ""))
        if tag == "table":
            self.table_id = a.get("id") or f"table_{len(self.table_rows)}"
            self.table_rows[self.table_id] = 0
            self.table_data[self.table_id] = []
        if tag == "thead": self.in_head = True
        if tag == "tbody": self.in_body = True
        if tag == "tr" and self.table_id:
            if self.in_head: self.current_header = []
            elif self.in_body:
                self.table_rows[self.table_id] += 1
                self.current_row = []
        if tag in {"td", "th"} and self.table_id:
            self.in_cell, self.current_cell = True, []
    def handle_endtag(self, tag):
        if tag in {"script", "style"} and self.skip: self.skip -= 1
        if tag == "nav": self.in_nav = False
        if tag in {"td", "th"} and self.in_cell:
            cell = " ".join("".join(self.current_cell).split())
            if self.in_head and self.current_header is not None:
                self.current_header.append(cell)
            elif self.current_row is not None:
                self.current_row.append(cell)
            self.in_cell, self.current_cell = False, []
        if tag == "tr" and self.table_id:
            if self.in_head and self.current_header is not None:
                self.table_headers[self.table_id] = self.current_header
                self.current_header = None
            elif self.in_body and self.current_row is not None:
                self.table_data[self.table_id].append(self.current_row)
                self.current_row = None
        if tag == "thead": self.in_head = False
        if tag == "tbody": self.in_body = False
        if tag == "table": self.table_id, self.in_body = None, False
    def handle_data(self, data):
        if not self.skip:
            self.text.append(data)
            if self.in_cell: self.current_cell.append(data)

def manifest_path(root: Path) -> Path:
    return root / "04_Methods_Reproducibility/delivery_manifest.tsv"


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def row_count(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return max(0, sum(1 for _ in handle) - 1)

def validate(root: Path, sample_id: str) -> dict:
    root = root.resolve()
    report = root / "01_Report" / f"{sample_id}_scPerturb_CRISPR_report.html"
    if not report.is_file() or report.stat().st_size == 0: raise FileNotFoundError(report)
    official_report = root / "05_DRAGEN_Report" / f"{sample_id}_DRAGEN_report.html"
    parity_path = root / "04_Methods_Reproducibility/validation/dragen_report_parity.tsv"
    html = report.read_text(encoding="utf-8")
    parser = ReportParser(); parser.feed(html)
    visible = " ".join(parser.text)
    checks = []
    def check(name, ok, detail): checks.append({"name": name, "status": "PASS" if bool(ok) else "FAIL", "detail": str(detail)})
    required_dirs = ["01_Report", "02_Figures", "03_Tables", "04_Methods_Reproducibility"]
    check("numbered_report_layout", all((root / d).is_dir() for d in required_dirs), required_dirs)
    check("report_navigation", parser.sections == SECTION_IDS and parser.navigation == SECTION_IDS,
          f"sections={parser.sections}; navigation={parser.navigation}")
    broken, external = [], []
    for href in parser.hrefs:
        parsed = urlparse(href)
        if parsed.scheme in {"http", "https"} or parsed.netloc: external.append(href)
        elif not href.startswith(("#", "mailto:", "data:")) and parsed.path:
            linked = (report.parent / unquote(parsed.path.split("#", 1)[0])).resolve()
            # This validator creates the manifest after checking the HTML links.
            if linked == manifest_path(root) and not linked.exists():
                continue
            if not linked.is_file(): broken.append(href)
    external.extend(x for x in parser.sources if not x.startswith("data:"))
    check("offline_assets_and_links", not broken and not external, f"broken={broken}; external={external}")
    try:
        validate_offline_report(official_report, sample_id)
        official_ok, official_detail = True, f"{official_report.relative_to(root)} ({official_report.stat().st_size} bytes)"
    except (OSError, ValueError) as exc:
        official_ok, official_detail = False, str(exc)
    check("official_dragen_report_offline", official_ok, official_detail)
    parity_rows = read_tsv(parity_path) if parity_path.is_file() else []
    from scperturb.dragen.alignment import SPECS, verify_alignment
    expected_overlap = {(scope, label) for _, scope, label, *_ in SPECS}
    try:
        aligned_count = verify_alignment(root)
        alignment_ok, alignment_detail = True, f"{aligned_count} scoped official metrics and exact barcode-rank coordinates"
    except (OSError, ValueError, KeyError) as exc:
        alignment_ok, alignment_detail = False, str(exc)
    check("dragen_report_shared_metrics", len(parity_rows) == len(expected_overlap)
          and all(row.get("sample_id") == sample_id and row.get("status") == "PASS" for row in parity_rows),
          f"rows={len(parity_rows)}; required={len(expected_overlap)}")
    check("official_metrics_in_final_html", alignment_ok, alignment_detail)
    check("logo_present", (root / "02_Figures/epi-logo.png").is_file() and "alt='Epi logo'" in html, "local logo source and embedded report logo")
    check("core_figures_present", all((root / "02_Figures/figures" / name).is_file() for name in EVIDENCE_FILES)
          and all(src.startswith("data:image/") for src in parser.sources if src.lower().endswith((".png", ".jpg", ".svg"))),
          f"evidence_and_official_files={EVIDENCE_FILES}; embedded_sources={len(parser.sources)}")
    icm_files = ["icm_expression_overview.png", "icm_guide_capture_overview.png"]
    icm_pdf = root / "04_Methods_Reproducibility/icm/analysis_summary.pdf"
    icm_alt_count = sum("ICM" in alt for alt in parser.image_alts)
    check("icm_perturbseq_exports", all((root / "02_Figures/figures" / name).is_file() for name in icm_files)
          and icm_pdf.is_file() and icm_pdf.stat().st_size > 0
          and "Official ICM Perturb-seq" in visible and icm_alt_count >= 2,
          f"images={icm_files}; summary_pdf={icm_pdf.is_file()}; ICM_alt_labels={icm_alt_count}")
    replaced_figures = {"fig02_gex_crispr_qc.png", "fig08_cell_state_composition.png",
                        "fig09_global_cell_state_umap.png", "fig10_guide_assignment_umap.png",
                        "fig11_ntc_vs_perturbation_umap.png", "fig12_representative_target_umap.png"}
    figure_manifest_rows = read_tsv(root / "03_Tables/figure_manifest.tsv")
    listed_figures = {row.get("file", "") for row in figure_manifest_rows}
    check("local_qc_and_exploratory_figures_replaced",
          not any((root / "02_Figures/figures" / name).exists() for name in replaced_figures)
          and not (listed_figures & replaced_figures),
          f"removed={sorted(replaced_figures)}")
    check("compact_body_tables", all(n <= 30 for n in parser.table_rows.values()), parser.table_rows)
    project = json.loads((root / "03_Tables/project_summary.json").read_text(encoding="utf-8"))
    classes = read_tsv(root / "03_Tables/target_evidence_class_summary.tsv")
    class_counts = {r["evidence_class"]: int(r["target_count"]) for r in classes}
    assignment = sum(int(project[k]) for k in ("singlet_cells", "multi_guide_cells", "unassigned_cells"))
    check("assignment_partition", assignment == int(project["filtered_cells"]),
          f"assignment={assignment}; passing={project['filtered_cells']}")
    check("evidence_partition", sum(class_counts.values()) == int(project["observed_targets"]),
          f"classes={class_counts}; observed={project['observed_targets']}")
    check("core_evidence_counts_visible",
          f"SUPPORTED 为 {class_counts.get('SUPPORTED', 0)} 个" in visible
          and f"PARTIAL_SUPPORT 为 {class_counts.get('PARTIAL_SUPPORT', 0)} 个" in visible,
          f"SUPPORTED={class_counts.get('SUPPORTED', 0)}; PARTIAL_SUPPORT={class_counts.get('PARTIAL_SUPPORT', 0)}")
    leaked_tokens = [token for token in ("{class_counts.get('SUPPORTED', 0)}",
                                          "{class_counts.get('PARTIAL_SUPPORT', 0)}",
                                          "pd.DataFrame", "Traceback (most recent call last)") if token in visible]
    check("no_unrendered_template_or_python_text", not leaked_tokens, f"found={leaked_tokens}")
    escaped_entities = re.findall(r"&(?:amp;)?(?:lt|gt|quot|#\d+);", visible)
    check("no_visible_double_escaped_entities", not escaped_entities,
          f"found={escaped_entities[:10]}")
    partial_confidence = read_tsv(root / "03_Tables/target_confidence.tsv")
    partial_confidence = [r for r in partial_confidence if r["confidence_class"] == "PARTIAL_SUPPORT"]
    target_summaries = {r["target_gene"]: r for r in read_tsv(root / "03_Tables/target_level_summary.tsv")}
    target_coverage = {r["target_gene"]: r for r in read_tsv(root / "03_Tables/target_coverage_full.tsv")}
    detail_headers = parser.table_headers.get("representative-effects", [])
    detail_rows = [dict(zip(detail_headers, row)) for row in parser.table_data.get("representative-effects", [])]
    rendered_candidates = {r.get("Target", ""): r for r in detail_rows}
    candidate_fields = {
        "Singlet cells": lambda r: str(int(float(r["target_singlet_cells"]))),
        "Active guides": lambda r: str(int(float(r["active_guides"]))),
        "Active guide IDs": lambda r: target_coverage[r["target_gene"]]["active_sgRNA_ids"],
        "Median guide rho": lambda r: f"{float(r['median_pairwise_spearman']):.3f}",
        "Direction concordance": lambda r: f"{100 * float(r['median_directional_concordance']):.2f}%",
        "Primary response gene": lambda r: r["representative_gene"],
        "Signed rank-biserial": lambda r: f"{float(r['representative_rank_biserial']):.3f}",
        "Target expression": lambda r: f"{100 * float(r['representative_pct_target']):.2f}%",
        "NTC expression": lambda r: f"{100 * float(r['representative_pct_control']):.2f}%",
        "Within-target BH FDR": lambda r: f"{float(r['representative_fdr']):.3g}",
        "Evidence class": lambda r: r["evidence_class"],
    }
    candidate_values_match = len(partial_confidence) == 4 and all(
        target in rendered_candidates and target in target_summaries and target in target_coverage
        and all(rendered_candidates[target].get(label) == value(target_summaries[target])
                for label, value in candidate_fields.items())
        for target in (r["target_gene"] for r in partial_confidence)
    )
    check("all_partial_targets_and_metrics_visible", candidate_values_match,
          f"partial_targets={[r['target_gene'] for r in partial_confidence]}; rendered={list(rendered_candidates)}")
    coverage = read_tsv(root / "03_Tables/target_coverage_full.tsv")
    eligible = read_tsv(root / "03_Tables/target_level_summary.tsv")
    methods = json.loads((root / "04_Methods_Reproducibility/methods/reanalysis.v3.json").read_text(encoding="utf-8"))
    min_cells = int(methods.get("min_target_cells", 20))
    eligible_from_coverage = sum(int(float(r["singlet_cells"])) >= min_cells for r in coverage)
    check("target_eligibility_matches_coverage", len(coverage) == int(project["observed_targets"])
          and len(eligible) == int(project["eligible_targets"]) == eligible_from_coverage,
          f"coverage={len(coverage)}; eligible={len(eligible)}; recomputed={eligible_from_coverage}")
    ntc = project["ntc_qc"]
    ntc_guides = read_tsv(root / "03_Tables/ntc_guide_qc.tsv")
    check("ntc_guide_summary_matches_source", len(ntc_guides) == int(ntc["ntc_guide_features"])
          and sum(str(r["active_for_effects"]).lower() in {"true", "1"} for r in ntc_guides) == int(ntc["active_ntc_guides_ge5_singlets"])
          and sum(int(r["singlet_cells"]) for r in ntc_guides) == int(ntc["ntc_singlet_cells"]),
          f"guides={len(ntc_guides)}; summary={ntc}")
    prep = json.loads((root / "04_Methods_Reproducibility/validation/v3_preparation.json").read_text(encoding="utf-8"))
    check("h5ad_mex_alignment", prep.get("h5ad_order_and_state", {}).get("h5ad_barcode_order_matches_filtered_mex") is True
          and prep.get("h5ad_order_and_state", {}).get("h5ad_barcode_order_matches_assignment") is True
          and prep.get("h5ad_order_and_state", {}).get("h5ad_gene_feature_order_matches_filtered_mex") is True,
          prep.get("h5ad_order_and_state"))
    check("exploratory_umap_replaced_by_icm", "ICM expression and cell composition overview" in " ".join(parser.image_alts)
          and "ICM CRISPR guide capture overview" in " ".join(parser.image_alts)
          and "ICM expression and cell composition overview" in visible
          and "ICM CRISPR guide capture overview" in visible
          and "GEX 与 CRISPR 捕获质量" not in visible
          and not any(x in visible for x in ("Existing UMAP with Leiden clusters", "Guide assignment UMAP",
                                             "NTC vs perturbation UMAP", "Representative target UMAP")),
          "ICM expression/capture figures replace local QC and exploratory UMAP panels")
    check("scientific_limits", all(x in visible for x in ["biological replicate inference", "UNKNOWN", "exploratory", "cannot prove editing success", "cannot prove editing failure"]),
          "Single sample, unknown perturbation mode, exploratory pathway and cell-state interpretation limits")
    check("leiden_cluster_language", "Leiden clusters are unannotated grouping labels" in visible
          and not any(x in visible for x in ["Cell-state composition", "Cell-state sensitivity", "细胞状态组成", "细胞状态敏感性"]),
          "Unannotated Leiden labels are described as clusters, not annotated cell states")
    check("exploratory_pathways", "Exploratory" in visible
          and "Pathway significance does not independently establish" in visible
          and "Reactome enrichment is exploratory for every target" in visible
          and "not a validated mechanism" in visible,
          "Reactome enrichment remains exploratory, including for low-evidence targets")
    check("no_unsupported_target_claim", not any(x in visible.lower() for x in ["validated target", "therapeutic target", "successful knockout", "knockout efficiency"]),
          "No unsupported validation or editing-efficiency language")
    pct_tokens = re.findall(r"(?<![\w.])([0-9][0-9,]*(?:\.[0-9]+)?)%", visible)
    bad_pct = [x for x in pct_tokens if "." not in x or len(x.rsplit(".", 1)[1]) != 2 or float(x.replace(",", "")) > 100]
    check("percentage_format", not bad_pct, f"percent_count={len(pct_tokens)}; invalid={bad_pct}")
    primary_de = root / "03_Tables/differential_expression/target_level_DE_all.tsv.gz"
    reactome = root / "03_Tables/pathway/reactome_ranking_sensitivity.tsv.gz"
    check("full_result_tables_retained", primary_de.is_file() and reactome.is_file()
          and row_count(primary_de) > 0 and row_count(reactome) > 0,
          f"DE_rows={row_count(primary_de) if primary_de.is_file() else 0}; pathway_rows={row_count(reactome) if reactome.is_file() else 0}")
    manifest = root / "04_Methods_Reproducibility/delivery_manifest.tsv"
    # The manifest is generated by this validator, so validate its intended path set
    # against all current delivery files (including the validation files written below).
    expected_validation = {
        "04_Methods_Reproducibility/validation/report_v3_2_1_validation.json",
        "04_Methods_Reproducibility/docs/REPORT_V3_2_1_VALIDATION.md",
    }
    current = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p != manifest}
    current |= expected_validation
    if manifest.is_file():
        listed = {row["path"] for row in read_tsv(manifest) if row.get("bytes") != "NON_RECURSIVE"}
        check("delivery_manifest_exact", listed == current,
              f"listed={len(listed)}; expected={len(current)}; missing={sorted(current-listed)[:5]}; extra={sorted(listed-current)[:5]}")
    else:
        check("delivery_manifest_exact", True, f"generated from {len(current)} delivery paths")
    failed = [c["name"] for c in checks if c["status"] != "PASS"]
    result = {"result": "PASS" if not failed else "FAIL", "sample_id": sample_id,
              "report_path": report.relative_to(root).as_posix(), "report_bytes": report.stat().st_size,
              "chapter_count": len(parser.sections), "summary": {"passing_cells": int(project["filtered_cells"]),
              "observed_targets": int(project["observed_targets"]), "eligible_targets": int(project["eligible_targets"]),
              "evidence_class_counts": class_counts, "ntc_singlet_cells": int(ntc["ntc_singlet_cells"])},
              "checks": checks}
    validation_dir = root / "04_Methods_Reproducibility/validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    (validation_dir / "report_v3_2_1_validation.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    lines = ["# V3.2.1 delivery validation", "", f"Result: **{result['result']}**", "",
             f"Report: {result['report_path']} ({result['report_bytes']:,} bytes); chapters={result['chapter_count']}.", "",
             "| Check | Status | Detail |", "|---|---|---|"]
    for item in checks:
        detail = str(item["detail"]).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {item['name']} | **{item['status']}** | {detail} |")
    (root / "04_Methods_Reproducibility/docs/REPORT_V3_2_1_VALIDATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_manifest(root)
    return result

def write_manifest(root: Path) -> None:
    manifest = manifest_path(root)
    files = sorted(p for p in root.rglob("*") if p.is_file() and p != manifest)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["path", "bytes", "format", "report_revision", "provenance"])
        for path in files:
            rel = path.relative_to(root).as_posix()
            writer.writerow([rel, path.stat().st_size, path.suffix.lstrip(".") or "text", "v3.2.1", "generated-scperturb-delivery"] )
        writer.writerow(["04_Methods_Reproducibility/delivery_manifest.tsv", "NON_RECURSIVE", "tsv", "v3.2.1", "self-inventory-entry"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()
    result = validate(args.root, args.sample_id)
    print(json.dumps({"result": result["result"], "checks": len(result["checks"]), "report_bytes": result["report_bytes"]}, ensure_ascii=False, indent=2))
    return 0 if result["result"] == "PASS" else 2

if __name__ == "__main__":
    raise SystemExit(main())
