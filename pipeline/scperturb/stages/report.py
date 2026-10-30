#!/usr/bin/env python3
"""Prepare, render, validate, and publish the standalone V3.2.1 report package."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from scperturb.common import delivery_name, load_json, output_manifest, sample_row, write_completion
from scperturb.icm import HANDOFF_DIRNAME, install_returned_assets
from scperturb.validation.delivery import validate
from scperturb.validation.review import prepare_review, finish_review

REPLACED_QC_AND_EXPLORATORY_FIGURES = {
    "fig02_gex_crispr_qc.png",
    "fig08_cell_state_composition.png",
    "fig09_global_cell_state_umap.png",
    "fig10_guide_assignment_umap.png",
    "fig11_ntc_vs_perturbation_umap.png",
    "fig12_representative_target_umap.png",
}


def remove_replaced_figures(root: Path) -> None:
    """Remove local QC/exploratory plots replaced by ICM exports from report assets."""
    figures = root / "02_Figures/figures"
    tables = root / "03_Tables"
    if not figures.is_dir():
        figures = root / "figures"
        tables = root / "tables"
    for filename in REPLACED_QC_AND_EXPLORATORY_FIGURES:
        (figures / filename).unlink(missing_ok=True)
    manifest = tables / "figure_manifest.tsv"
    if manifest.is_file():
        with manifest.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = [row for row in reader if row.get("file") not in REPLACED_QC_AND_EXPLORATORY_FIGURES]
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def move_tree(source: Path, destination: Path) -> None:
    if not source.exists(): raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))


def move_contents(source: Path, destination: Path) -> None:
    if not source.is_dir(): raise FileNotFoundError(source)
    destination.mkdir(parents=True, exist_ok=True)
    for item in sorted(source.iterdir()):
        target = destination / item.name
        if target.exists(): raise FileExistsError(target)
        shutil.move(str(item), str(target))
    source.rmdir()


def generate_official_report(run: dict, sample_id: str, delivery_root: Path, source_root: Path) -> tuple[Path, Path]:
    official_report = delivery_root / "05_DRAGEN_Report" / f"{sample_id}_DRAGEN_report.html"
    parity = delivery_root / "04_Methods_Reproducibility/validation/dragen_report_parity.tsv"
    command = [
        sys.executable, "-m", "scperturb.dragen.reports", "--mode", "scperturb",
        "--sample-id", sample_id,
        "--dragen-output-dir", str(source_root / "02_dragen_scrna_crispr/04_tool_output/dragen"),
        "--reports-root", "/opt/dragen-reports",
        "--manifest", "/opt/dragen-reports/manifests/scrna.json",
        "--report", str(official_report), "--generate-only",
        "--temp-root", str(Path(run["tmp_root"]) / sample_id / "dragen_report_build"),
        "--guide-qc-json", str(source_root / "03_cell_guide_qc/01_tables/cell_guide_qc.json"),
    ]
    subprocess.run(command, check=True, cwd=run["project_root"])
    return official_report, parity


def read_parity(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def refresh_report(run: dict, sample_id: str, stage: Path) -> int:
    """Refresh QC and report from completed delivery; no analysis stage is called."""
    from scperturb.dragen.alignment import prepare_alignment
    from scperturb.reporting.html import build
    final = stage / "04_tool_output" / delivery_name(run, sample_id)
    if not (stage / "05_metadata/stage_complete.json").is_file():
        raise RuntimeError("Report-only refresh requires a completed report stage")
    work = Path(run["tmp_root"]) / sample_id / "report_build" / f"official_alignment_{time.time_ns()}"
    package = work / "delivery"
    shutil.copytree(final, package)
    source = Path(run["results_root"]) / sample_id
    prepare_alignment(package, source / "02_dragen_scrna_crispr/04_tool_output/dragen",
                      source / "03_cell_guide_qc/01_tables/cell_guide_qc.json")
    remove_replaced_figures(package)
    install_returned_assets(
        source / "02_dragen_scrna_crispr/04_tool_output" / HANDOFF_DIRNAME,
        package,
    )
    for fragment in [package / "04_Methods_Reproducibility/assets/read_flow_sankey_section.html",
                     package / "03_Tables/read_flow/read_flow_sankey_section.html"]:
        if fragment.is_file():
            text = fragment.read_text().replace("默认展示 Read-pair 守恒流；Bin50 守恒流通过下拉框切换。", "按 GEX reads、CRISPR reads、molecules 和 cells 分统计单位展示；不同单位不可直接比较带宽。")
            fragment.write_text(text.replace("reads 到 bin50 是多对一聚合", "reads 到 molecules / cells 是聚合转换"))
    docs = package / "04_Methods_Reproducibility/docs"
    methods_path = package / "04_Methods_Reproducibility/methods/analysis_methods.md"
    if methods_path.is_file():
        methods = methods_path.read_text(encoding="utf-8")
        methods = methods.replace("cell-state sensitivity reads", "Leiden-cluster sensitivity reads")
        methods = methods.replace("## Cell-state composition sensitivity",
                                  "## Leiden-cluster composition and stratified sensitivity")
        methods = methods.replace(
            "描述 target 与 NTC singlet cells 在输入 H5AD Leiden labels 中的组成差异。对至少有 20 个 target cells 落在共有 Leiden 状态、且至少 2 个状态各有 ≥5 个 NTC singlets 的 target，基于 target cells 在共有状态内的比例，对各状态内 top 100 primary-effect gene 的 signed rank-biserial effects 加权平均；比较 primary 与分层效应方向。该结果不报告 p/FDR，仅为单样本探索性敏感性分析。条件不足的 target 标为 NOT_TESTABLE。",
            "描述 target 与 NTC singlet cells 在输入 H5AD 的未注释 Leiden cluster labels 中的组成差异；这些 cluster 不直接代表已注释细胞类型或细胞状态。对至少有 20 个 target cells 落在共有 clusters、且至少 2 个 clusters 各有 ≥5 个 NTC singlets 的 target，基于 target cells 在共有 clusters 内的比例，对各 cluster 内 top 100 primary-effect gene 的 signed rank-biserial effects 加权平均；比较 primary 与 cluster-stratified 效应方向。该结果不报告 p/FDR，仅为单样本探索性敏感性分析。条件不足的 target 标为 NOT_TESTABLE。",
        )
        methods = methods.replace(
            "全部 Reactome 信号均为探索性结果。",
            "全部 Reactome 信号均为探索性结果；低证据或非 SUPPORTED target 的 pathway enrichment 不构成已验证机制。",
        )
        methods_path.write_text(methods, encoding="utf-8")
    claim_map_path = package / "03_Tables/claim_evidence_map.tsv"
    if claim_map_path.is_file():
        with claim_map_path.open(encoding="utf-8", newline="") as handle:
            claim_rows = list(csv.DictReader(handle, delimiter="\t"))
            claim_fields = list(claim_rows[0]) if claim_rows else []
        changed = False
        for row in claim_rows:
            if row.get("claim") == "Cell state":
                row["claim"] = "Leiden cluster composition"
                row["evidence_basis"] = "Unannotated input H5AD Leiden clusters and existing embedding"
                row["interpretation_boundary"] = "Descriptive/exploratory grouping labels; do not directly identify cell states or establish mechanism"
                changed = True
        if changed:
            with claim_map_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=claim_fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(claim_rows)
    (docs / "REPORT_DRAGEN_ALIGNMENT_CHANGELOG.md").write_text(
        "# DRAGEN official alignment refresh\n\n"
        "- Replaced primary QC display with scoped official GEX, CRISPR, mapping and guide metrics.\n"
        "- Independent native-metric and barcode-summary parity covers all displayed official metrics.\n"
        "- Added official barcode-rank coordinates, static and interactive views.\n"
        "- Guide percentages explicitly distinguish passing-cell and feature-positive denominators.\n"
        "- Fixed NTC integer display and removed unrelated Bin50/read-pair wording.\n"
        "- DRAGEN, DE, evidence, guide/LOOGO, Reactome and UMAP analyses were not rerun or changed.\n")
    build(package, Path(run["project_root"]) / "pipeline/assets/plotly-3.5.1.min.js")
    return publish_report(run, sample_id, stage, package, work, status="refreshed")


def delivery_outputs(final: Path, sample_id: str) -> list[Path]:
    return [final / relative for relative in (
        f"01_Report/{sample_id}_scPerturb_CRISPR_report.html",
        f"05_DRAGEN_Report/{sample_id}_DRAGEN_report.html",
        "03_Tables/dragen_official_metrics.tsv", "03_Tables/dragen_barcode_rank.json",
        "02_Figures/figures/fig14_dragen_barcode_rank.png",
        "02_Figures/figures/icm_expression_overview.png",
        "02_Figures/figures/icm_guide_capture_overview.png",
        "04_Methods_Reproducibility/icm/analysis_summary.pdf",
        "04_Methods_Reproducibility/icm/icm_export_manifest.tsv",
        "04_Methods_Reproducibility/validation/dragen_report_parity.tsv",
        "04_Methods_Reproducibility/validation/report_v3_2_1_validation.json",
        "04_Methods_Reproducibility/docs/REPORT_V3_2_1_VALIDATION.md",
        "04_Methods_Reproducibility/delivery_manifest.tsv",
        f"04_Methods_Reproducibility/review/{sample_id}_DELIVERY_REVIEW.md",
        f"04_Methods_Reproducibility/review/{sample_id}_DELIVERY_REVIEW.tsv",
        f"04_Methods_Reproducibility/review/{sample_id}_DELIVERY_REVIEW.json",
        "04_Methods_Reproducibility/review/SC_PERTURB_DELIVERY_REVIEW_CRITERIA.md",
    )]


def publish_report(run: dict, sample_id: str, stage: Path, package: Path,
                   work: Path, *, status: str) -> int:
    """Validate a staged package, then replace delivery with rollback on move failure."""
    manifest = package / "04_Methods_Reproducibility/delivery_manifest.tsv"
    manifest.unlink(missing_ok=True)
    review = prepare_review(run, sample_id, package)
    validation = validate(package, sample_id)
    review = finish_review(package, review, validation)
    if validation["result"] != "PASS":
        raise RuntimeError(f"Report validation failed; existing delivery retained: {validation['checks']}")
    final = stage / "04_tool_output" / delivery_name(run, sample_id)
    final.parent.mkdir(parents=True, exist_ok=True)
    archive = work / "previous_delivery"
    had_previous = final.exists()
    if had_previous:
        shutil.move(str(final), str(archive))
    try:
        shutil.move(str(package), str(final))
    except Exception:
        if had_previous:
            shutil.move(str(archive), str(final))
        raise
    outputs = delivery_outputs(final, sample_id)
    stage_manifest = output_manifest(stage, outputs)
    marker = write_completion(stage, run, sample_id, [*outputs, stage_manifest],
                              interface="05-scperturb-report-official-aligned-v1")
    print(json.dumps({"status": status, "checks": len(validation["checks"]),
                      "delivery_review": review["overall_result"],
                      "shared_metrics": len(read_parity(outputs[5])),
                      "report": str(outputs[0]), "report_bytes": outputs[0].stat().st_size,
                      "chapters": validation["chapter_count"], "marker": str(marker),
                      "previous_delivery": str(archive) if had_previous else None}))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage-dir", required=True, type=Path)
    parser.add_argument("--refresh-report", action="store_true", help="Refresh completed report only; reuse all analysis results")
    args = parser.parse_args()
    run = load_json(args.config)
    sample = sample_row(run, args.sample_id)
    stage = args.stage_dir.resolve()
    stage.mkdir(parents=True, exist_ok=True)
    metadata = stage / "05_metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    marker_path = metadata / "stage_complete.json"
    report_final = stage / "04_tool_output" / delivery_name(run, args.sample_id) / "01_Report" / f"{args.sample_id}_scPerturb_CRISPR_report.html"
    source_root = Path(run["results_root"]) / args.sample_id
    if args.refresh_report:
        return refresh_report(run, args.sample_id, stage)
    if marker_path.is_file():
        previous = json.loads(marker_path.read_text(encoding="utf-8"))
        if previous.get("status") == "complete" and report_final.is_file():
            return refresh_report(run, args.sample_id, stage)
        raise RuntimeError(f"incomplete report marker exists; inspect before rerun: {marker_path}")

    primary_meta = source_root / "04_primary_reanalysis/05_metadata/reanalysis_root.json"
    if not primary_meta.is_file(): raise FileNotFoundError(primary_meta)
    primary = json.loads(primary_meta.read_text(encoding="utf-8"))
    v2_root = Path(primary["v2_root"])
    if not (v2_root / "05_metadata/stage_complete.json").is_file():
        raise RuntimeError("primary V2 reanalysis is not complete")

    project_tmp = Path(run["tmp_root"]) / args.sample_id / "report_build"
    task_root = project_tmp / f"{primary['fingerprint']}_{time.time_ns()}"
    prepared = task_root / "prepared"
    package = task_root / "delivery"
    task_root.mkdir(parents=True, exist_ok=False)
    control = str(load_json(run["analysis_config"])["guide_reference"]["control_gene"])
    analysis = load_json(run["analysis_config"])["pathway"]
    prepare_args = [sys.executable, "-m", "scperturb.reporting.prepare",
                    "--sample-id", args.sample_id, "--control-gene", control,
                    "--v2-root", str(v2_root), "--source-root", str(source_root),
                    "--output-root", str(prepared), "--runtime-root",
                    str(Path(run["tmp_root"]) / args.sample_id / "v3_runtime" / primary["fingerprint"]),
                    "--reactome-gmt", str(analysis["gmt_file"])]
    subprocess.run(prepare_args, check=True, cwd=run["project_root"])
    figure_args = [sys.executable, "-m", "scperturb.reporting.figures",
                   "--root", str(prepared), "--source-root", str(source_root),
                   "--sample-id", args.sample_id, "--control-gene", control]
    subprocess.run(figure_args, check=True, cwd=run["project_root"])

    (package / "01_Report").mkdir(parents=True, exist_ok=True)
    (package / "05_DRAGEN_Report").mkdir(parents=True, exist_ok=True)
    (package / "02_Figures").mkdir(parents=True, exist_ok=True)
    methods_root = package / "04_Methods_Reproducibility"
    (methods_root / "docs").mkdir(parents=True, exist_ok=True)
    (methods_root / "assets").mkdir(parents=True, exist_ok=True)
    (methods_root / "validation").mkdir(parents=True, exist_ok=True)
    move_tree(prepared / "tables", package / "03_Tables")
    move_tree(prepared / "figures", package / "02_Figures/figures")
    move_tree(prepared / "methods", methods_root / "methods")
    move_contents(prepared / "05_metadata", methods_root / "validation")
    for name in ("REPORT_V3_CHANGELOG.md", "REANALYSIS_SUMMARY.md"):
        move_tree(prepared / name, methods_root / "docs" / name)
    flow_fragment = package / "03_Tables/read_flow/read_flow_sankey_section.html"
    if flow_fragment.is_file():
        shutil.copy2(flow_fragment, methods_root / "assets/read_flow_sankey_section.html")
    shutil.copy2(Path(run["project_root"]) / "pipeline/assets/epi-logo.png", package / "02_Figures/epi-logo.png")
    plotly = Path(run["project_root"]) / "pipeline/assets/plotly-3.5.1.min.js"
    if not plotly.is_file(): raise FileNotFoundError(plotly)

    from scperturb.reporting.html import build, render_screen_funnel
    official_report, parity = generate_official_report(run, args.sample_id, package, source_root)
    from scperturb.dragen.alignment import prepare_alignment
    prepare_alignment(package, source_root / "02_dragen_scrna_crispr/04_tool_output/dragen",
                      source_root / "03_cell_guide_qc/01_tables/cell_guide_qc.json")
    render_screen_funnel(package)
    remove_replaced_figures(package)
    install_returned_assets(
        source_root / "02_dragen_scrna_crispr/04_tool_output" / HANDOFF_DIRNAME,
        package,
    )
    report = build(package, plotly)
    if not report.is_file() or report.stat().st_size < 100_000:
        raise RuntimeError("V3.2.1 report build produced an unexpectedly small or missing HTML")

    return publish_report(run, args.sample_id, stage, package, task_root, status="complete")

if __name__ == "__main__":
    raise SystemExit(main())
