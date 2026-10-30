#!/usr/bin/env python3
"""Build the standalone, offline customer report for scPerturb CRISPR V3."""

from __future__ import annotations

import argparse
import base64
import html
import json
import math
import re
from pathlib import Path

import pandas as pd


SECTIONS = [
    ("overview", "1 数据集概况与质量（Dataset overview / QC）"),
    ("guide_assignment", "2 Guide 检测与分配（Guide assignment）"),
    ("screen_usability", "3 筛选可用性（Screen usability）"),
    ("coverage", "4 靶点覆盖（Target coverage）"),
    ("effects", "5 靶点表达效应（Target-level expression effects）"),
    ("robustness", "6 Guide 稳健性评估（Guide robustness）"),
    ("cell_state", "7 ICM Perturb-seq 与 cluster 敏感性分析（ICM Perturb-seq / cluster sensitivity）"),
    ("pathways", "8 Reactome 通路分析（Exploratory pathways）"),
    ("evidence", "9 靶点证据汇总（Target evidence summary）"),
    ("limitations", "10 结果限制（Limitations）"),
    ("reproducibility", "11 下载与复现文件（Reproducibility）"),
]
from scperturb.reporting.style import NPG_EVIDENCE

CLASS_COLORS = dict(NPG_EVIDENCE)
PERCENT_ALREADY_0_TO_100 = {"pct_root", "pct_parent", "fraction_of_view_root"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--plotly-js", required=True, type=Path)
    return parser.parse_args()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def render_screen_funnel(root: Path) -> None:
    """Redraw the existing funnel counts with neutral final-stage wording."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tables = root / "03_Tables"
    project = read_json(tables / "project_summary.json")
    stages = read_tsv(tables / "screen_target_funnel.tsv")
    funnel = dict(zip(stages["stage"], pd.to_numeric(stages["count"], errors="raise").astype(int)))
    passing = int(project["filtered_cells"])
    cells = [
        ("Singlet", int(project["singlet_cells"]), "#00A087"),
        ("Multi-guide", int(project["multi_guide_cells"]), "#4DBBD5"),
        ("Unassigned", int(project["unassigned_cells"]), "#B8BDC6"),
    ]
    target_labels = [
        "Observed targets", "≥20 singlet cells", "≥2 active guides",
        "Guide-balanced", "Meets minimal consistency gate",
    ]
    target_keys = ["Observed targets", "≥20 singlet cells", "≥2 active guides", "Guide-balanced", "Guide-consistent"]
    target_values = [int(funnel[key]) for key in target_keys]
    target_colors = ["#8491B4", "#4DBBD5", "#00A087", "#3C5488", "#E64B35"]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 12, "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.6), gridspec_kw={"width_ratios": [0.95, 1.35]})
    ax = axes[0]
    vals = [v / passing * 100 for _, v, _ in cells]
    y = list(range(len(cells)))
    bars = ax.barh(y, vals, color=[c for _, _, c in cells], height=0.62)
    ax.set_yticks(y, [label for label, _, _ in cells])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100], ["0.00%", "25.00%", "50.00%", "75.00%", "100.00%"])
    ax.set_xlabel("Share of passing cells")
    ax.set_title(f"Assignment composition (n={passing:,})", pad=10)
    for bar, (_, count, _) in zip(bars, cells):
        ax.text(min(bar.get_width() + 1.2, 81), bar.get_y() + bar.get_height() / 2,
                f"{count:,} ({count/passing*100:.2f}%)", va="center", ha="left", fontweight="bold", fontsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#e7ebef", linewidth=0.7)
    ax.set_axisbelow(True)

    ax = axes[1]
    y = list(range(len(target_labels)))
    bars = ax.barh(y, target_values, color=target_colors, height=0.62)
    ax.set_yticks(y, target_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, max(target_values) * 1.13)
    ax.set_xlabel("Observed targets")
    ax.set_title("Target analysis funnel", pad=10)
    for bar, count in zip(bars, target_values):
        ax.text(max(bar.get_width() + 8, 8), bar.get_y() + bar.get_height() / 2,
                f"{count:,}", va="center", ha="left", fontweight="bold", fontsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#e7ebef", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout(pad=1.15, w_pad=2.2)
    out = root / "02_Figures/figures/fig01_screen_usability_funnel.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def esc(value) -> str:
    if value is None:
        return "—"
    text = str(value)
    if text.lower() in {"nan", "none", "nat", "na"}:
        return "—"
    return html.escape(text, quote=True)


def numeric(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def fmt(value, key="") -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "—"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    k = key.lower()
    bool_fields = {"stable", "loogo_stable", "direction_concordant", "main_candidate", "active_for_effects"}
    if str(value).lower() in {"true", "1"} and (k in bool_fields or "eligible" in k or "candidate" in k):
        return "Yes"
    if str(value).lower() in {"false", "0"} and (k in bool_fields or "eligible" in k or "candidate" in k):
        return "No"
    x = numeric(value)
    if x is None:
        return esc(value)
    if any(token in k for token in ("p_value", "fdr", "q_value", "fwer")):
        return f"{x:.3g}"
    if re.search(r"(?:^|_)(?:fraction|share|concordance|pct|percent|rate)(?:_|$)", k):
        percent_value = x if k in PERCENT_ALREADY_0_TO_100 else 100 * x
        return f"{percent_value:.2f}%"
    if re.search(r"(?:^|_)(?:cell|cells|guide|guides|target|targets|cluster|clusters|gene|genes|pathway|pathways|feature|features|count|counts|row|rows|permutation|permutations|overlap|umi|umis|pair|pairs|symbol|symbols|size|reads|molecules|barcodes|barcode)(?:_|$)", k):
        return f"{int(round(x)):,}"
    if any(token in k for token in ("effect", "rho", "spearman", "nes", "es", "log2fc", "distance", "jensen", "median_abs")):
        return f"{x:.3f}"
    return f"{x:,.3f}" if abs(x) < 1000 else f"{x:,.0f}"


def metric_card(label: str, value: str, sub: str = "") -> str:
    return f"<article class='metric-card'><span>{esc(label)}</span><strong>{esc(value)}</strong><small>{esc(sub)}</small></article>"


def table_html(frame: pd.DataFrame, columns: list[tuple[str, str]], table_id: str,
               page_size: int = 20, paginated: bool = True, max_rows: int | None = None) -> str:
    if max_rows is not None:
        frame = frame.head(max_rows)
    cls = "data-table paginated-table" if paginated else "data-table"
    attr = f" data-page-size='{int(page_size)}'" if paginated else ""
    headers = "".join(f"<th scope='col'>{esc(label)}</th>" for _, label in columns)
    body = []
    for record in frame.to_dict("records"):
        cells = []
        for key, _ in columns:
            value = record.get(key, "")
            if key in {"confidence_class", "evidence_class"}:
                cls_name = str(value)
                color = CLASS_COLORS.get(cls_name, "#777777")
                cell = f"<span class='class-chip' style='--class-color:{color}'>{esc(cls_name)}</span>"
            else:
                format_key = key
                if table_id == "threshold-sensitivity" and key in {"targets", "changed_from_primary", "SUPPORTED", "PARTIAL_SUPPORT", "DISCORDANT", "GUIDE_DOMINATED", "LOW_GUIDE_SUPPORT", "LOW_CELL", "INCONCLUSIVE"}:
                    format_key = "target_count"
                if table_id == "ntc-metrics" and key == "value":
                    format_key = "cell_count"
                elif key == "value" and record.get("metric"):
                    format_key = re.sub(r"[^a-zA-Z0-9]+", "_", str(record["metric"]).lower()).strip("_")
                # fmt() escapes non-numeric source strings; decode that formatting
                # layer before applying the single HTML-context escape here.
                cell = esc(html.unescape(fmt(value, format_key)))
            cells.append(f"<td>{cell}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    if not body:
        body = [f"<tr><td colspan='{len(columns)}' class='empty'>当前筛选条件下无数据。</td></tr>"]
    table = (f"<div class='table-scroll'><table id='{esc(table_id)}' class='{cls}'{attr}>"
             f"<thead><tr>{headers}</tr></thead><tbody>{''.join(body)}</tbody></table></div>")
    if paginated:
        table += (f"<div class='table-controls'><label>筛选 <input type='search' data-filter-for='{esc(table_id)}' placeholder='输入关键词'></label>"
                  f"<div class='pager' data-pager-for='{esc(table_id)}'><button type='button' data-prev>上一页</button>"
                  f"<span data-status></span><button type='button' data-next>下一页</button></div></div>")
    return table


def figure_card(root: Path, filename: str, title: str, caption: str, alt: str, show_title: bool = True) -> str:
    image = base64.b64encode((root / "02_Figures" / "figures" / filename).read_bytes()).decode("ascii")
    heading = f"<h3>{esc(title)}</h3>" if show_title else ""
    return (f"<article class='figure-card'><button class='figure-open' type='button' data-open-figure aria-label='Enlarge {esc(title)}'>"
            f"<img src='data:image/png;base64,{image}' alt='{esc(alt)}'></button>{heading}<p>{esc(caption)}</p></article>")


def format_table_columns(frame: pd.DataFrame, preferred: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return [(key, label) for key, label in preferred if key in frame.columns]


def build(root: Path, plotly_path: Path) -> Path:
    from scperturb.dragen.alignment import alignment_html
    official_qc_html, official_guide_html = alignment_html(root)
    tables = root / "03_Tables"
    project = read_json(tables / "project_summary.json")
    summary = read_json(tables / "reanalysis_summary.json")
    config = read_json(root / "04_Methods_Reproducibility/methods/reanalysis.v3.json")
    target_summary = read_tsv(tables / "target_level_summary.tsv")
    target_full = read_tsv(tables / "target_coverage_full.tsv")
    funnel = read_tsv(tables / "screen_target_funnel.tsv")
    coverage_dist = read_tsv(tables / "target_coverage_distribution.tsv")
    class_summary = read_tsv(tables / "target_evidence_class_summary.tsv")
    representatives = read_tsv(tables / "representative_target_effects.tsv")
    ntc_guides = read_tsv(tables / "ntc_guide_qc.tsv")
    consistency = read_tsv(tables / "guide_consistency_summary.tsv")
    loo = read_tsv(tables / "leave_one_guide_out_summary.tsv")
    sensitivity = read_tsv(tables / "target_threshold_sensitivity.tsv")
    state_comp = read_tsv(tables / "cell_state_composition.tsv")
    state_dist = read_tsv(tables / "cell_state_distance_to_ntc.tsv")
    state_overlap = read_tsv(tables / "cell_state_overlap.tsv")
    state_effect = read_tsv(tables / "cell_state_effect_sensitivity.tsv")
    gsea_qc = read_tsv(tables / "gsea_ranking_sensitivity_qc.tsv") if (tables / "gsea_ranking_sensitivity_qc.tsv").exists() else pd.DataFrame()
    pathway_candidates = read_tsv(tables / "pathway/reactome_exploratory_candidates.tsv") if (tables / "pathway/reactome_exploratory_candidates.tsv").exists() else pd.DataFrame()
    pathway_compare = pd.DataFrame()  # Full results remain delivered as TSV.GZ.
    self_expression = read_tsv(tables / "target_self_expression_effect.tsv")
    claims = read_tsv(tables / "claim_evidence_map.tsv")
    figure_manifest = read_tsv(tables / "figure_manifest.tsv")
    confidence = read_tsv(tables / "target_confidence.tsv")
    read_flow_tables = tables / "read_flow"
    read_flow_nodes = read_tsv(read_flow_tables / "read_flow_nodes.tsv")
    read_flow_transforms = read_tsv(read_flow_tables / "read_flow_transformations.tsv")
    umap_target_counts = read_tsv(tables / "umap_representative_target_counts.tsv")
    de_plot_genes = read_tsv(tables / "representative_DE_plot_genes.tsv")
    metrics = project
    ntc = project["ntc_qc"]
    funnel_values = {r["stage"]: int(r["count"]) for r in funnel.to_dict("records")}
    passing = int(metrics["filtered_cells"])
    sensitivity_changes = dict(zip(sensitivity["scenario"],
                                   pd.to_numeric(sensitivity["changed_from_primary"], errors="coerce").fillna(0).astype(int)))
    sensitivity_text = (
        f"在本次单因素敏感性网格中，把 target 最低细胞数从 20 改为 10 或 30，分别有 {sensitivity_changes.get('min_target_cells=10', 0)} 或 {sensitivity_changes.get('min_target_cells=30', 0)} 个 target 改变分类；"
        f"活跃 guide 最低数改为 1 或 3 时分别有 {sensitivity_changes.get('min_active_guides=1', 0)} 或 {sensitivity_changes.get('min_active_guides=3', 0)} 个改变；"
        f"最大 guide 占比阈值改为 0.4 或 0.6 时分别有 {sensitivity_changes.get('max_guide_fraction=0.4', 0)} 或 {sensitivity_changes.get('max_guide_fraction=0.6', 0)} 个改变。"
        "所测试的 Spearman、方向一致性和 LOOGO 数值边界未改变分类。结果对细胞覆盖与 guide 支持门槛更敏感；主报告仍保留预设 primary thresholds。"
    )

    feature_cell = read_tsv(tables / "cell_guide_feature_umi.tsv.gz")
    assign = read_tsv(tables / "guide_assignment_all.tsv.gz")
    feature_cell["total_feature_umi"] = pd.to_numeric(feature_cell["total_feature_umi"], errors="coerce")
    feature_cell = feature_cell.merge(assign[["barcode", "guide_status"]], on="barcode", how="left", validate="one_to_one")
    feature_medians = feature_cell.groupby("guide_status")["total_feature_umi"].median().to_dict()

    class_counts = {r["evidence_class"]: int(float(r["target_count"]))
                    for r in class_summary.to_dict("records")}
    class_text = " · ".join(f"{r['evidence_class']} {fmt(r['target_count'], 'target_count')}" for r in class_summary.to_dict("records") if int(float(r["target_count"])) > 0)
    cards = "".join([
        metric_card("Passing cells", fmt(metrics["filtered_cells"], "cells"), "DRAGEN passing-cell 统计口径"),
        metric_card("Singlet cells", fmt(metrics["singlet_cells"], "cells"), "仅分配到 1 条 guide"),
        metric_card("Unassigned cells", f"{fmt(metrics['unassigned_cells'], 'cells')} · {fmt(metrics['unassigned_rate_of_filtered_cells'], 'fraction')}", "占 passing cells"),
        metric_card("NTC singlet cells", fmt(ntc["ntc_singlet_cells"], "cells"), "target-vs-NTC DE reference"),
        metric_card("Biological samples", "1", "单样本分析"),
        metric_card("Observed targets", fmt(metrics["observed_targets"], "targets"), "至少 1 个 singlet assignment"),
        metric_card("Eligible targets", fmt(metrics["eligible_targets"], "targets"), "≥20 singlet cells；不是 hit 判定"),
        metric_card("Final evidence", f"{class_counts.get('SUPPORTED', 0)} SUPPORTED · {class_counts.get('PARTIAL_SUPPORT', 0)} PARTIAL_SUPPORT", "其余为 EXPLORATORY 汇总级别"),
    ])
    funnel_counts = dict(zip(funnel["stage"].astype(str), pd.to_numeric(funnel["count"], errors="raise").astype(int)))
    supported_n = class_counts.get("SUPPORTED", 0)
    partial_n = class_counts.get("PARTIAL_SUPPORT", 0)
    conclusion = (f"<div class='conclusion'><strong>项目结论</strong><p>原生 DRAGEN CRISPR 流程已完成并获得有效 guide 分配。"
                 f"{passing:,} 个质控通过细胞中，单 guide 细胞为 {int(metrics['singlet_cells']):,} ({fmt(metrics['singlet_rate_of_filtered_cells'], 'fraction')})，"
                 f"未分配细胞为 {int(metrics['unassigned_cells']):,} ({fmt(metrics['unassigned_rate_of_filtered_cells'], 'fraction')})。"
                 f"{int(metrics['observed_targets']):,} 个观察靶点中，{int(metrics['eligible_targets']):,} 个达到 ≥20 个 singlet 细胞；"
                 f"{funnel_counts.get('Guide-balanced', 0):,} 个同时具有 ≥2 条活跃 guide 且最大 guide 占比 ≤0.5。"
                 f"按预设规则，{supported_n:,} 个靶点为 SUPPORTED，{partial_n:,} 个为 PARTIAL_SUPPORT。主要限制是靶点细胞覆盖、NTC 对照规模、单个生物学样本和未知扰动模式。"
                 f"因此本报告分别呈现细胞级表达差异与多 guide 靶点证据强度；达到分析条件或出现表达效应本身不等于靶点验证。</p></div>")

    fig = {row["file"]: row for row in figure_manifest.to_dict("records")}
    figure_titles = {
        "fig01_screen_usability_funnel.png": "细胞分配与筛选可用性（Cell assignment / target funnel）",
        "fig02_gex_crispr_qc.png": "GEX 与 CRISPR 捕获质量（GEX / CRISPR capture QC）",
        "fig03_target_coverage.png": "靶点覆盖分布与高覆盖靶点（Target coverage）",
        "fig04_multi_guide_robustness.png": "多 guide 一致性与稳健性（Multi-guide robustness）",
        "fig05_representative_perturbation_effects.png": "代表性扰动相关表达效应（Perturbation-associated effects）",
        "fig06_pathway_ranking_qc.png": "Reactome 排名质量控制（Ranking QC）",
        "fig07_reactome_ranking_sensitivity.png": "Reactome 方向敏感性（Direction sensitivity）",
        "fig08_cell_state_composition.png": "NTC 与高覆盖靶点的 Leiden cluster 组成（Leiden-cluster composition）",
        "fig13_representative_de_effects.png": "代表性 signed rank-biserial 效应（Signed effect）",
    }
    def fig_card(filename: str, alt: str | None = None, show_title: bool = True) -> str:
        row = fig[filename]
        title = figure_titles.get(filename, row["title"])
        caption = row["caption"]
        if filename == "fig01_screen_usability_funnel.png":
            caption = ("Passing-cell guide classes are mutually exclusive. Target stages use ≥20 singlet cells, ≥2 active guides, "
                       "largest-guide share ≤0.5, then the existing rho and direction-concordance criteria. "
                       "The last stage is labelled Meets minimal consistency gate; it includes PARTIAL_SUPPORT and is not equivalent to SUPPORTED evidence.")
        return figure_card(root, filename, title, caption, alt or title, show_title=show_title)

    qc = pd.DataFrame([
        ("GEX input reads", metrics["gex_reads"], "reads"),
        ("CRISPR feature reads", metrics["crispr_feature_reads"], "reads"),
        ("GEX reads in passing cells", float(metrics["fraction_reads_in_passing_cells"]), "fraction"),
        ("Feature reads in passing cells", float(metrics["feature_fraction_reads_in_passing_cells"]), "fraction"),
        ("Sequencing saturation", float(metrics["sequencing_saturation"]), "fraction; descriptive only"),
        ("GEX MAPQ ≥40", float(metrics["mapq40_rate_gex"]), "fraction"),
        ("GEX bases ≥Q30", float(metrics["q30_rate_gex"]), "fraction"),
        ("Median genes per passing cell", metrics["median_genes_per_cell"], "genes"),
        ("Median molecules per passing cell", metrics["median_molecules_per_cell"], "molecules"),
        ("Feature molecules in passing cells", metrics["feature_molecules_passing_cells"], "molecules"),
        ("Median feature UMIs per singlet cell", feature_medians.get("singlet", math.nan), "UMIs"),
    ], columns=["metric", "value", "unit_or_denominator"])
    qc["value"] = [
        fmt(value, 'fraction') if 'fraction' in str(unit) else fmt(value, 'count')
        for value, unit in zip(qc["value"], qc["unit_or_denominator"])
    ]
    ntc_metrics = pd.DataFrame([
        ("NTC guide features", ntc["ntc_guide_features"]),
        ("Unique NTC feature-positive cells", ntc["ntc_feature_positive_cells_unique"]),
        ("NTC singlet assignments", ntc["ntc_cells_with_assignment_singlet"]),
        ("NTC multi-guide assignments", ntc["ntc_cells_with_assignment_multi"]),
        ("Active NTC guides (≥5 singlets)", ntc["active_ntc_guides_ge5_singlets"]),
        ("NTC singlet cells used as DE reference", ntc["ntc_singlet_cells"]),
    ], columns=["metric", "value"])
    ntc_metrics["value"] = [f"{int(float(value)):,}" for value in ntc_metrics["value"]]
    ntc_funnel = (
        "<div class='ntc-flow' aria-label='NTC guide detection and assignment summary'>"
        f"<div class='flow-step'><span>NTC guide features</span><strong>{fmt(ntc['ntc_guide_features'], 'features')}</strong><small>feature definitions</small></div>"
        "<b aria-hidden='true'>→</b>"
        f"<div class='flow-step'><span>Unique feature-positive cells</span><strong>{fmt(ntc['ntc_feature_positive_cells_unique'], 'cells')}</strong><small>cell-level detection</small></div>"
        "<b aria-hidden='true'>→</b>"
        f"<div class='flow-step'><span>Assigned NTC cells</span><strong>{fmt(ntc['ntc_cells_with_assignment_singlet'], 'cells')} singlet · {fmt(ntc['ntc_cells_with_assignment_multi'], 'cells')} multi-guide</strong><small>assignment outcomes</small></div>"
        "<b aria-hidden='true'>→</b>"
        f"<div class='flow-step'><span>Active NTC guides</span><strong>{fmt(ntc['active_ntc_guides_ge5_singlets'], 'guides')}</strong><small>≥5 singlet cells / guide</small></div>"
        "</div><p class='small-note'>各节点使用不同统计单位或筛选口径，不构成严格的守恒漏斗；图示不推断阶段间数量差异的原因。</p>"
    )

    target_singlet_total = int(pd.to_numeric(target_full["singlet_cells"], errors="coerce").fillna(0).sum())
    observed_targets = int(metrics["observed_targets"])
    eligible_targets = int(metrics["eligible_targets"])
    dataflow_rows = [
        ("GEX library input", fmt(metrics["gex_reads"], "reads"), "reads", "—", "GEX library read total."),
        ("GEX reads associated with passing cells", "—", "reads", f"GEX input reads ({fmt(metrics['gex_reads'], 'reads')})",
         f"{fmt(metrics['fraction_reads_in_passing_cells'], 'fraction')} of GEX input reads."),
        ("CRISPR feature-library input", fmt(metrics["crispr_feature_reads"], "reads"), "reads", "—", "CRISPR feature-library read total."),
        ("CRISPR feature reads associated with passing cells", "—", "reads", f"CRISPR input reads ({fmt(metrics['crispr_feature_reads'], 'reads')})",
         f"{fmt(metrics['feature_fraction_reads_in_passing_cells'], 'fraction')} of CRISPR feature reads."),
        ("GEX molecules across raw MEX barcodes", fmt(metrics["gex_molecules_all_barcodes"], "molecules"), "molecules", "All barcodes represented in raw GEX MEX", "Molecule counts; not read counts."),
        ("GEX molecules in passing cells", fmt(metrics["gex_molecules_passing_cells"], "molecules"), "molecules",
         f"All raw GEX molecules ({fmt(metrics['gex_molecules_all_barcodes'], 'molecules')})",
         f"{fmt(float(metrics['gex_molecules_passing_cells'])/float(metrics['gex_molecules_all_barcodes']), 'fraction')} of all raw GEX molecules."),
        ("CRISPR feature molecules across raw MEX barcodes", fmt(metrics["feature_molecules_all_barcodes"], "molecules"), "molecules", "All barcodes represented in raw feature MEX", "Guide-feature molecule counts; not read counts."),
        ("CRISPR feature molecules in passing cells", fmt(metrics["feature_molecules_passing_cells"], "molecules"), "molecules",
         f"All raw feature molecules ({fmt(metrics['feature_molecules_all_barcodes'], 'molecules')})",
         f"{fmt(float(metrics['feature_molecules_passing_cells'])/float(metrics['feature_molecules_all_barcodes']), 'fraction')} of all raw feature molecules."),
        ("DRAGEN passing cells", fmt(metrics["filtered_cells"], "cells"), "cells", "Cells called passing by DRAGEN", "Cell-level denominator used for assignment fractions."),
        ("Singlet guide cells", fmt(metrics["singlet_cells"], "cells"), "cells", f"Passing cells ({fmt(metrics['filtered_cells'], 'cells')})",
         f"{fmt(metrics['singlet_rate_of_filtered_cells'], 'fraction')} of passing cells; one assigned guide."),
        ("Multi-guide cells", fmt(metrics["multi_guide_cells"], "cells"), "cells", f"Passing cells ({fmt(metrics['filtered_cells'], 'cells')})",
         f"{fmt(metrics['multi_guide_rate_of_filtered_cells'], 'fraction')} of passing cells."),
        ("Unassigned cells", fmt(metrics["unassigned_cells"], "cells"), "cells", f"Passing cells ({fmt(metrics['filtered_cells'], 'cells')})",
         f"{fmt(metrics['unassigned_rate_of_filtered_cells'], 'fraction')} of passing cells."),
        ("Perturbation target-level singlet cells", fmt(target_singlet_total, "cells"), "cells", f"All {fmt(observed_targets, 'target_count')} observed perturbation targets",
         "Sum of per-target singlet cells; NTC is excluded from the perturbation-target table."),
        ("Observed perturbation targets", fmt(observed_targets, "target_count"), "targets", "Targets with ≥1 singlet assignment", "Target count; not a cell count."),
        ("Eligible targets (≥20 singlet cells)", fmt(eligible_targets, "target_count"), "targets", f"Observed targets ({fmt(observed_targets, 'target_count')})",
         f"{fmt(eligible_targets/observed_targets, 'fraction')} of observed targets meet the minimum-cell criterion."),
        ("Active guides (≥5 singlet cells)", fmt(metrics["active_guides_ge5_singlets"], "guides"), "guides", "Guide features with ≥5 singlet cells", "Guide count; guide support is not a target count."),
        ("Median genes per passing cell", fmt(metrics["median_genes_per_cell"], "genes"), "genes/cell", f"Passing cells ({fmt(metrics['filtered_cells'], 'cells')})", "Per-cell median; gene features are not molecules."),
        ("Median GEX molecules per passing cell", fmt(metrics["median_molecules_per_cell"], "molecules"), "molecules/cell", f"Passing cells ({fmt(metrics['filtered_cells'], 'cells')})", "Per-cell median UMI/molecule count."),
    ]
    dataflow = pd.DataFrame(dataflow_rows, columns=["metric", "value", "unit", "denominator", "interpretation"])
    transform_rows = []
    for row in read_flow_transforms.to_dict("records"):
        transform_rows.append({
            "stage_relation": row["label"],
            "source_measure": f"{row['source_label']}: {fmt(row['source_value'], 'count')} {row['source_unit']}",
            "target_measure": f"{row['target_label']}: {fmt(row['target_value'], 'count')} {row['target_unit']}",
            "relation": row["relation"], "interpretation": row["reason"],
        })
    transformations = pd.DataFrame(transform_rows)
    root_by_view = read_flow_nodes[read_flow_nodes["stage"].astype(str).eq("0")].set_index("view_id")["label"].to_dict()
    root_value_by_view = read_flow_nodes[read_flow_nodes["stage"].astype(str).eq("0")].set_index("view_id")["value"].to_dict()
    read_flow_nodes["denominator"] = read_flow_nodes["view_id"].map(
        lambda view: f"{root_by_view.get(view, 'view root')} ({fmt(root_value_by_view.get(view, 0), 'count')} {read_flow_nodes.loc[read_flow_nodes['view_id'].eq(view) & read_flow_nodes['stage'].astype(str).eq('0'), 'unit'].iloc[0]})"
    )
    read_flow_nodes = read_flow_nodes.rename(columns={"label": "stage_label", "pct_root": "fraction_of_view_root"})
    read_flow_nodes["metric"] = read_flow_nodes["stage_label"]
    read_flow_nodes["value"] = read_flow_nodes["value"].map(lambda value: fmt(value, "count"))

    distribution = coverage_dist.copy()
    distribution["fraction_of_observed_targets"] = pd.to_numeric(distribution["fraction_of_observed_targets"], errors="coerce")
    target_top = target_summary.sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True]).head(20)
    target_full_sorted = target_full.sort_values(["singlet_cells", "target_gene"], ascending=[False, True])
    target_summary = target_summary.sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True])
    representative_source = representatives.copy()
    partial_targets = set(target_summary.loc[
        target_summary["evidence_class"].eq("PARTIAL_SUPPORT"), "target_gene"
    ].astype(str))
    partial_rows = target_summary[target_summary["target_gene"].astype(str).isin(partial_targets)].copy()
    representatives = pd.concat([
        partial_rows,
        representative_source[~representative_source["target_gene"].astype(str).isin(partial_targets)],
    ], ignore_index=True, sort=False)
    guide_ids = target_full.set_index("target_gene")["active_sgRNA_ids"].to_dict()
    representatives["active_sgRNA_ids"] = representatives["target_gene"].map(guide_ids).fillna("")
    reps = representatives
    de_summary_cols = format_table_columns(target_summary, [
        ("target_gene", "Target"), ("target_singlet_cells", "Singlet 细胞数"), ("active_guides", "活跃 guide 数"),
        ("largest_guide_share", "最大 guide 占比"), ("de_genes_fdr_lt_0_05", "DE 基因数（FDR<0.05）"),
        ("strongest_effect_abs_rank_biserial", "最大 |效应量|"), ("strongest_effect_gene", "最大效应对应基因"),
        ("median_pairwise_spearman", "Guide rho 中位数"), ("median_directional_concordance", "方向一致性"),
        ("loogo_stable", "LOOGO 稳定"), ("evidence_class", "证据分类")])
    representative_cols = format_table_columns(reps, [
        ("target_gene", "Target"), ("target_singlet_cells", "Singlet 细胞数"), ("active_guides", "活跃 guide 数"),
        ("largest_guide_share", "最大 guide 占比"), ("representative_gene", "代表性基因"),
        ("representative_rank_biserial", "Signed rank-biserial 效应量"), ("representative_pct_target", "Target 表达比例"),
        ("representative_pct_control", "NTC 表达比例"), ("representative_fdr", "Within-target BH FDR"),
        ("evidence_class", "证据分类"), ("selection_reason", "纳入展示的考虑因素")])

    primary_tie = pd.to_numeric(read_tsv(tables / "gsea_ranking_qc.tsv")["tie_fraction"], errors="coerce")
    secondary_tie = pd.to_numeric(gsea_qc.get("secondary_tie_fraction", pd.Series(dtype=float)), errors="coerce")
    ranking_cards = "".join([
        metric_card("Primary tie fraction", fmt(primary_tie.median(), 'fraction'), f"Range {fmt(primary_tie.min(), 'fraction')}–{fmt(primary_tie.max(), 'fraction')} across {len(primary_tie):,} targets"),
        metric_card("Primary zero-effect fraction", fmt(pd.to_numeric(read_tsv(tables / 'gsea_ranking_qc.tsv')['zero_effect_fraction'], errors='coerce').median(), 'fraction'), "Median across targets"),
        metric_card("Genes retained", fmt(pd.to_numeric(read_tsv(tables / 'gsea_ranking_qc.tsv')['expressed_genes_before_symbol_collapse'], errors='coerce').median(), "genes"), "Median expressed genes before symbol collapse"),
        metric_card("Secondary tie fraction", fmt(secondary_tie.median(), 'fraction') if len(secondary_tie.dropna()) else '—', "Normalized mean-expression difference"),
    ])

    # Evidence classes use one named palette throughout report text and figures.
    legend = "".join(f"<span class='legend-item'><i style='background:{CLASS_COLORS.get(r['evidence_class'], '#777777')}'></i>{esc(r['evidence_class'])}</span>"
                     for r in class_summary.to_dict("records"))
    state_testable = state_effect[state_effect.status.eq("TESTED_EXPLORATORY")]
    state_summary_text = (f"{len(state_testable):,} of {len(target_summary):,} eligible targets met the declared Leiden-cluster overlap criterion and had a top-gene stratified effect sensitivity; "
                          f"{int((state_effect.status == 'NOT_TESTABLE').sum()):,} were NOT_TESTABLE. Leiden clusters are unannotated grouping labels.")
    umap_targets_text = "、".join(
        f"{r['target_gene']}（{fmt(r['target_singlet_cells'], 'cells')} cells）"
        for r in umap_target_counts.to_dict("records")
    )
    plotly_js = plotly_path.read_text(encoding="utf-8")
    plotly_loader = "<script>(()=>{const p=atob('" + base64.b64encode(plotly_js.encode()).decode("ascii") + "');const s=document.createElement('script');s.text=p;document.head.appendChild(s)})();</script>"

    qc_rows = qc
    top_coverage = target_summary.sort_values(
        ["target_singlet_cells", "target_gene"], ascending=[False, True]
    ).head(15)
    coverage_cols = format_table_columns(top_coverage, [
        ("target_gene", "Target"), ("target_singlet_cells", "Singlet cells"),
        ("active_guides", "Active guides"), ("largest_guide_share", "Largest-guide share"),
        ("evidence_class", "Evidence class"),
    ])
    de_gene_rows = de_plot_genes.head(20)
    representative_cols = format_table_columns(representatives, [
        ("target_gene", "Target"), ("target_singlet_cells", "Singlet cells"),
        ("active_guides", "Active guides"), ("active_sgRNA_ids", "Active guide IDs"),
        ("largest_guide_share", "Largest-guide share"),
        ("median_pairwise_spearman", "Median guide rho"),
        ("median_directional_concordance", "Direction concordance"),
        ("de_genes_fdr_lt_0_05", "DE genes, FDR<0.05"), ("representative_gene", "Primary response gene"),
        ("representative_rank_biserial", "Signed rank-biserial"),
        ("representative_pct_target", "Target expression"), ("representative_pct_control", "NTC expression"),
        ("representative_fdr", "Within-target BH FDR"), ("evidence_class", "Evidence class"),
    ])
    target_lookup = target_summary.set_index("target_gene").to_dict("index")
    pathway_display = pathway_candidates.head(10).copy()
    if not pathway_display.empty:
        pathway_display["target_singlet_cells"] = pathway_display["target_gene"].map(
            {k: v.get("target_singlet_cells", "") for k, v in target_lookup.items()})
        pathway_display["active_guides"] = pathway_display["target_gene"].map(
            {k: v.get("active_guides", "") for k, v in target_lookup.items()})
        pathway_display["evidence_class"] = pathway_display["target_gene"].map(
            {k: v.get("evidence_class", "") for k, v in target_lookup.items()})

    state_lookup = state_effect.set_index("target_gene").to_dict("index")
    pathway_lookup = {}
    if not pathway_candidates.empty:
        p_sorted = pathway_candidates.copy()
        p_sorted["_p"] = pd.to_numeric(p_sorted["primary_FDR"], errors="coerce")
        p_sorted = p_sorted.sort_values(["_p", "pathway_id"], na_position="last")
        pathway_lookup = p_sorted.drop_duplicates("target_gene").set_index("target_gene").to_dict("index")
    focus_targets = list(dict.fromkeys(
        representatives["target_gene"].tolist() + pathway_display["target_gene"].tolist()))
    evidence_rows = []
    for target in focus_targets:
        row = target_lookup.get(target, {})
        st = state_lookup.get(target, {})
        pw = pathway_lookup.get(target, {})
        original_class = row.get("evidence_class", "")
        reporting_tier = original_class if original_class in {"SUPPORTED", "PARTIAL_SUPPORT"} else "EXPLORATORY"
        effect_value = numeric(row.get("representative_rank_biserial", ""))
        effect_gene = row.get("representative_gene", "")
        evidence_rows.append({
            "target_gene": target,
            "target_singlet_cells": row.get("target_singlet_cells", ""),
            "active_guides": row.get("active_guides", ""),
            "de_genes_fdr_lt_0_05": row.get("de_genes_fdr_lt_0_05", ""),
            "strongest_effect": f"{effect_gene}: {effect_value:.3f}" if effect_gene and effect_value is not None else "—",
            "median_pairwise_spearman": row.get("median_pairwise_spearman", ""),
            "median_directional_concordance": row.get("median_directional_concordance", ""),
            "loogo_stable": row.get("loogo_stable", ""),
            "state_status": st.get("status", "NOT_TESTABLE"),
            "state_rho": st.get("spearman_primary_vs_state_adjusted", ""),
            "pathway_primary_fdr": pw.get("primary_FDR", ""),
            "pathway_secondary_fdr": pw.get("secondary_FDR", ""),
            "evidence_class": original_class,
            "reporting_tier": reporting_tier,
        })
    evidence_table = pd.DataFrame(evidence_rows)
    html_tables = {
        "qc": table_html(qc_rows, [("metric", "Metric"), ("value", "Value"),
                                    ("unit_or_denominator", "Unit / denominator")],
                         "qc-table", paginated=False),
        "ntc_metrics": table_html(ntc_metrics, [("metric", "NTC metric"), ("value", "Value")],
                                  "ntc-metrics", paginated=False),
        "coverage_bins": table_html(distribution, [("bin", "Singlet cells / target"),
                                                     ("target_count", "Targets"),
                                                     ("fraction_of_observed_targets", "Observed targets")],
                                    "coverage-bins", paginated=False),
        "class_summary": table_html(class_summary, [
            ("evidence_class", "Evidence class"), ("target_count", "Targets"),
            ("fraction", "Observed targets"), ("rule", "Classification rule"),
        ], "class-summary", paginated=False),
        "top_coverage": table_html(top_coverage, coverage_cols, "top-coverage", paginated=False),
        "representative_de_genes": table_html(de_gene_rows, [
            ("target_gene", "Target"), ("target_singlet_cells", "Singlet cells"),
            ("active_guides", "Active guides"), ("gene_name", "Gene"),
            ("rank_biserial", "Signed rank-biserial"), ("fdr", "Within-target BH FDR"),
            ("pct_target", "Target expression"), ("pct_control", "NTC expression"),
            ("delta_pct", "Expression-fraction difference"),
        ], "representative-de-genes", paginated=False),
        "representatives": table_html(representatives, representative_cols,
                                      "representative-effects", paginated=False),
        "pathway": table_html(pathway_display, [
            ("target_gene", "Target"), ("target_singlet_cells", "Target cells"),
            ("active_guides", "Active guides"), ("evidence_class", "Evidence class"),
            ("pathway_id", "Reactome ID"),
            ("pathway_name", "Pathway"), ("primary_NES", "Primary NES"),
            ("primary_FDR", "Primary pathway FDR"), ("secondary_NES", "Sensitivity NES"),
            ("secondary_FDR", "Sensitivity pathway FDR"), ("direction_concordant", "Direction concordant"),
            ("coverage_fraction", "Coverage"),
            ("leading_edge_size", "Leading edge"),
        ], "reactome-candidates", paginated=False),
        "target_evidence": table_html(evidence_table, [
            ("target_gene", "Target"), ("target_singlet_cells", "Target cells"),
            ("active_guides", "Active guides"), ("de_genes_fdr_lt_0_05", "DE genes, FDR<0.05"),
            ("strongest_effect", "Representative effect"), ("median_pairwise_spearman", "Median guide rho"),
            ("median_directional_concordance", "Direction concordance"), ("loogo_stable", "LOOGO stable"),
            ("state_status", "Cluster-stratified status"), ("state_rho", "Primary / cluster-stratified rho"),
            ("pathway_primary_fdr", "Top Reactome primary FDR"),
            ("pathway_secondary_fdr", "Top Reactome sensitivity FDR"),
            ("evidence_class", "Original evidence class"), ("reporting_tier", "Reporting tier"),
        ], "target-evidence-summary", paginated=False),
    }
    state_testable = state_effect[state_effect.status.eq("TESTED_EXPLORATORY")]
    state_not_testable = int((state_effect.status == "NOT_TESTABLE").sum())
    state_rho = pd.to_numeric(state_testable["spearman_primary_vs_state_adjusted"], errors="coerce")
    state_change = pd.to_numeric(state_testable["median_abs_effect_change"], errors="coerce")
    state_summary_text = (
        f"{len(state_testable):,} targets were TESTED_EXPLORATORY and "
        f"{state_not_testable:,} were NOT_TESTABLE under the declared Leiden-cluster overlap criterion. "
        f"Among tested targets, median primary-versus-cluster-stratified Spearman correlation was "
        f"{fmt(state_rho.median(), 'spearman')}; median absolute effect change was "
        f"{fmt(state_change.median(), 'effect')}. Leiden clusters are unannotated grouping labels, "
        f"not annotated cell types or cell states."
    )
    consistency_targets = target_summary[
        pd.to_numeric(target_summary["active_guides"], errors="coerce").fillna(0) >= 2
    ]
    loo_stable = target_summary["loogo_stable"].astype(str).str.lower().isin(["true", "1", "yes"])
    rho_values = pd.to_numeric(consistency_targets["median_pairwise_spearman"], errors="coerce")
    concordance_values = pd.to_numeric(consistency_targets["median_directional_concordance"], errors="coerce")
    robustness_summary = (
        f"{len(consistency_targets):,} eligible targets had at least two active guides; "
        f"{int(funnel_values.get('Guide-balanced', 0)):,} met the existing guide-balance gate and "
        f"{int(funnel_values.get('Guide-consistent', 0)):,} passed guide-consistency criteria "
        f"(this does not equal SUPPORTED evidence). Median evaluated guide-guide Spearman rho "
        f"was {fmt(rho_values.median(), 'spearman')}; median directional concordance was "
        f"{fmt(concordance_values.median(), 'concordance')}. "
        f"{int(loo_stable.sum()):,} eligible targets were marked LOOGO-stable."
    )
    primary_qc = read_tsv(tables / "gsea_ranking_qc.tsv")
    primary_tie = pd.to_numeric(primary_qc["tie_fraction"], errors="coerce")
    primary_zero = pd.to_numeric(primary_qc["zero_effect_fraction"], errors="coerce")
    primary_retained = pd.to_numeric(primary_qc["expressed_genes_before_symbol_collapse"], errors="coerce")
    secondary_tie = pd.to_numeric(gsea_qc.get("secondary_tie_fraction", pd.Series(dtype=float)), errors="coerce")
    ranking_cards = "".join([
        metric_card("Primary ranking ties", fmt(primary_tie.median(), "fraction"),
                    f"Median across {len(primary_tie):,} targets"),
        metric_card("Zero-effect fraction", fmt(primary_zero.median(), "fraction"),
                    "Primary signed rank-biserial; median"),
        metric_card("Genes retained", fmt(primary_retained.median(), "genes"),
                    "Median after expression filtering"),
        metric_card("Secondary ranking ties", fmt(secondary_tie.median(), "fraction")
                    if len(secondary_tie.dropna()) else "—",
                    "Normalized mean-expression difference"),
    ])
    legend = "".join(
        f"<span class='legend-item'><i style='background:{CLASS_COLORS.get(r['evidence_class'], '#777777')}'></i>{esc(r['evidence_class'])}</span>"
        for r in class_summary.to_dict("records")
    )
    umap_targets_text = "、".join(
        f"{r['target_gene']}（{fmt(r['target_singlet_cells'], 'cells')} cells）"
        for r in umap_target_counts.to_dict("records")
    )
    plotly_js = plotly_path.read_text(encoding="utf-8")
    plotly_loader = "<script>(()=>{const p=atob('" + base64.b64encode(plotly_js.encode()).decode("ascii") + "');const s=document.createElement('script');s.text=p;document.head.appendChild(s)})();</script>"
    flow_fragment = (root / "04_Methods_Reproducibility/assets/read_flow_sankey_section.html").read_text(encoding="utf-8")
    flow_fragment = flow_fragment.replace("<section class='bio-flow-section' id='bio-read-flow'>",
                                          "<div class='bio-flow-section' id='bio-read-flow'>").replace("</section>", "</div>")
    flow_fragment = flow_fragment.replace("默认展示 Read-pair 守恒流；Bin50 守恒流通过下拉框切换。", "按 GEX reads、CRISPR reads、molecules 和 cells 分统计单位展示；不同单位不可直接比较带宽。")
    flow_fragment = flow_fragment.replace("reads 到 bin50 是多对一聚合", "reads 到 molecules / cells 是聚合转换")
    downloads = """
    <div class='two-col download-grid'>
      <article class='panel'><h3>Target results</h3><a href='../03_Tables/target_coverage_full.tsv'>Full target coverage</a> · <a href='../03_Tables/target_level_summary.tsv'>Eligible target summary</a></article>
      <article class='panel'><h3>Differential expression</h3><a href='../03_Tables/differential_expression/target_level_DE_all.tsv.gz'>Full gene-level DE (TSV.GZ)</a> · <a href='../03_Tables/differential_expression/target_DE_summary.tsv'>DE summary</a></article>
      <article class='panel'><h3>Guide robustness</h3><a href='../03_Tables/guide_level_effects.tsv.gz'>Guide effects</a> · <a href='../03_Tables/guide_consistency_summary.tsv'>Guide consistency</a> · <a href='../03_Tables/leave_one_guide_out_summary.tsv'>LOOGO</a></article>
      <article class='panel'><h3>Reactome</h3><a href='../03_Tables/pathway/reactome_GSEA_all.tsv.gz'>Primary</a> · <a href='../03_Tables/pathway/reactome_GSEA_mean_difference.tsv.gz'>Sensitivity</a> · <a href='../03_Tables/pathway/reactome_ranking_sensitivity.tsv.gz'>Combined</a></article>
      <article class='panel'><h3>Leiden cluster analysis</h3><a href='../03_Tables/cell_state_effect_sensitivity.tsv'>Cluster-stratified sensitivity</a> · <a href='../03_Tables/cell_state_composition.tsv'>Composition</a> · <a href='../03_Tables/umap_cell_metadata.tsv.gz'>UMAP metadata</a></article>
      <article class='panel'><h3>Methods / reproducibility</h3><a href='../04_Methods_Reproducibility/methods/analysis_methods.md'>Methods</a> · <a href='../04_Methods_Reproducibility/methods/reanalysis.v3.json'>Reanalysis config</a> · <a href='../04_Methods_Reproducibility/delivery_manifest.tsv'>Delivery manifest</a></article>
    </div>
    """
    class_counts = {r["evidence_class"]: int(float(r["target_count"]))
                    for r in class_summary.to_dict("records")}
    eligible_fraction = (float(metrics["eligible_targets"]) / float(metrics["observed_targets"])
                         if float(metrics["observed_targets"]) else math.nan)
    consistency_targets = target_summary[
        pd.to_numeric(target_summary["active_guides"], errors="coerce").fillna(0) >= 2
    ]
    consistency_rho = pd.to_numeric(consistency_targets["median_pairwise_spearman"], errors="coerce")
    consistency_direction = pd.to_numeric(consistency_targets["median_directional_concordance"], errors="coerce")
    concise_conclusion = (
        "<div class='conclusion'><strong>项目结论</strong><p>"
        "本次筛选的可解释靶点覆盖有限，guide 间一致性整体较弱。"
        f"在 {fmt(metrics['observed_targets'], 'targets')} 个 observed targets 中，"
        f"仅 {fmt(metrics['eligible_targets'], 'targets')} 个（{fmt(eligible_fraction, 'fraction')}）满足 ≥20 singlet cells 的分析门槛；"
        f"这 {len(consistency_targets):,} 个至少有 2 条活跃 guides 的 eligible targets 中，"
        f"仅 {funnel_counts.get('Guide-consistent', 0):,} 个通过既定一致性条件；guide-guide rho 中位数为 "
        f"{fmt(consistency_rho.median(), 'spearman')}，方向一致性中位数为 {fmt(consistency_direction.median(), 'concordance')}。"
        f"按预设 evidence rules，SUPPORTED 为 {class_counts.get('SUPPORTED', 0)} 个，"
        f"PARTIAL_SUPPORT 为 {class_counts.get('PARTIAL_SUPPORT', 0)} 个。"
        "Eligible 仅表示达到分析条件，不等于 CRISPR hit；PARTIAL_SUPPORT 也不是已验证靶点。"
        "单样本 DE 与 Reactome 结果仅用于描述和探索，不能确立因果机制。</p></div>"
    )
    outcome_flow = (
        "<div class='outcome-flow' aria-label='Target eligibility and evidence class summary'>"
        "<div class='outcome-path'>"
        f"<div class='outcome-step'><span>Observed targets</span><strong>{fmt(metrics['observed_targets'], 'targets')}</strong></div>"
        "<b aria-hidden='true'>→</b>"
        f"<div class='outcome-step'><span>Eligible (≥20 singlet cells)</span><strong>{fmt(metrics['eligible_targets'], 'targets')}</strong></div>"
        "</div><div class='outcome-outcomes'><span class='outcome-label'>Evidence classification across observed targets</span><div class='outcome-classes'>"
        f"<div><span>SUPPORTED</span><strong>{class_counts.get('SUPPORTED', 0)}</strong></div>"
        f"<div><span>PARTIAL_SUPPORT</span><strong>{class_counts.get('PARTIAL_SUPPORT', 0)}</strong></div>"
        f"<div><span>EXPLORATORY*</span><strong>{sum(v for k, v in class_counts.items() if k not in {'SUPPORTED', 'PARTIAL_SUPPORT'})}</strong></div>"
        "</div></div></div><p class='small-note'>Eligible 是分析门槛；证据分类覆盖所有 observed targets。EXPLORATORY* 是将既有非 SUPPORTED / PARTIAL_SUPPORT 细分类合并用于报告展示，不改变原始分类或规则。</p>"
    )
    phkb_row = target_lookup.get("PHKB", {})
    phkb_path = pathway_lookup.get("PHKB", {})
    phkb_context = ""
    if phkb_row and phkb_path:
        phkb_context = (
            f"PHKB 的 target coverage 为 {fmt(phkb_row.get('target_singlet_cells'), 'cells')} singlet cells、"
            f"{fmt(phkb_row.get('active_guides'), 'guides')} active guides，原始 evidence class 为 "
            f"{esc(phkb_row.get('evidence_class'))}；其 pathway FDR 应与该 target-level 支持度一并阅读。"
        )
    class_tier_cards = "".join([
        metric_card("SUPPORTED", fmt(class_counts.get("SUPPORTED", 0), "targets"), "当前分类结果"),
        metric_card("PARTIAL_SUPPORT", fmt(class_counts.get("PARTIAL_SUPPORT", 0), "targets"), "当前分类结果"),
        metric_card("EXPLORATORY", fmt(sum(v for k, v in class_counts.items() if k not in {"SUPPORTED", "PARTIAL_SUPPORT"}), "targets"), "报告展示汇总，不改变原始细分类"),
    ])
    sections_data = {
        "overview": (
            f"{concise_conclusion}<div class='metric-grid'>{cards}</div>{outcome_flow}"
            "<div class='exploratory-banner'><strong>Official ICM Perturb-seq</strong><span>表达与细胞组成图为 Illumina Connected Multiomics 默认 Perturb-seq workflow 导出的静态图；报告离线显示，不加载在线 ICM 内容。</span></div>"
            f"{fig_card('icm_expression_overview.png')}"
            f"{official_qc_html}"
            f"{figure_card(root, 'fig14_dragen_barcode_rank.png', 'DRAGEN barcode rank', 'Direct reuse of official Cells / Background coordinates; cell calls and molecule threshold are unchanged.', 'Official barcode rank')}"
            "<details class='appendix print-skip'><summary>交互查看 barcode rank</summary><div id='official-barcode-rank-plot' style='height:420px'></div></details>"
            f"<script type='application/json' id='official-barcode-rank-data'>{(tables / 'dragen_barcode_rank.json').read_text()}</script>"
            "<script>document.addEventListener('DOMContentLoaded',()=>{const traces=JSON.parse(document.getElementById('official-barcode-rank-data').textContent);const el=document.getElementById('official-barcode-rank-plot');el.parentElement.addEventListener('toggle',()=>{if(el.parentElement.open){Plotly.react(el,traces.map((t,i)=>({x:t.data.map(p=>p.x),y:t.data.map(p=>p.y),name:t.label,type:'scatter',mode:'lines',line:{color:['#00A087','#B8BDC6'][i]}})),{xaxis:{type:'log',title:{text:'Barcode rank'}},yaxis:{type:'log',title:{text:'GEX molecules per barcode'}},margin:{t:20},paper_bgcolor:'#fff'},{responsive:true});}});});</script>"
            "<p>统计单位：reads、molecules、barcodes 和 cells 不可跨单位作守恒比较。Read-flow 百分比按原始整数计数计算，可能与官方已舍入比例有精度差异；不据此反推计数。</p>"
            f"<details class='appendix print-skip sankey-details' id='read-flow-sankey'><summary>查看分统计单位的 read-flow Sankey</summary>{flow_fragment}</details>"
        ),
        "guide_assignment": (
            "<p class='lead'>Guide assignment 将 passing cells 分为 singlet、multi-guide 与 unassigned；各类数量之和等于 passing cells。</p>"
            f"<div class='metric-grid compact'>{metric_card('Singlet cells', fmt(metrics['singlet_cells'], 'cells'), fmt(metrics['singlet_rate_of_filtered_cells'], 'fraction') + ' of passing cells')}"
            f"{metric_card('Multi-guide cells', fmt(metrics['multi_guide_cells'], 'cells'), fmt(metrics['multi_guide_rate_of_filtered_cells'], 'fraction') + ' of passing cells')}"
            f"{metric_card('Unassigned cells', fmt(metrics['unassigned_cells'], 'cells'), fmt(metrics['unassigned_rate_of_filtered_cells'], 'fraction') + ' of passing cells')}"
            f"{metric_card('Active guides', fmt(metrics['active_guides_ge5_singlets'], 'guides'), '≥5 singlet cells / guide')}</div>"
            "<h3>ICM CRISPR guide capture</h3>"
            f"{fig_card('icm_guide_capture_overview.png')}"
            f"{official_guide_html}"
            "<h3>NTC guide 检测与 assignment 概况（NTC assignment）</h3>"
            f"{ntc_funnel}{html_tables['ntc_metrics']}"
            f"<div class='notice compact-notice'><strong>NTC reference size is limited</strong>"
            f"<p>Target-vs-NTC cell-level DE、effect 与 pathway analysis 实际使用 {fmt(ntc['ntc_singlet_cells'], 'cells')} 个满足当前 assignment 与 singlet 条件的 NTC cells，而不是全部 feature-positive cells；结果需谨慎解释。</p></div>"
            "<p class='file-note'>每条 NTC guide 的完整支持度见 <a href='../03_Tables/ntc_guide_qc.tsv'>ntc_guide_qc.tsv</a>。</p>"
        ),
        "screen_usability": (
            f"<p class='lead'>细胞 assignment 与 target eligibility 分属不同分析层级；图中最后阶段改称 “Meets minimal consistency gate”。</p>"
            f"{fig_card('fig01_screen_usability_funnel.png')}"
            "<p class='small-note'>Meets minimal consistency gate 表示满足既有 guide consistency 条件；它不等于 SUPPORTED target，也不构成独立验证。</p>"
        ),
        "coverage": (
            f"<p class='lead'>{fmt(metrics['observed_targets'], 'targets')} 个 observed targets 中，"
            f"{fmt(metrics['eligible_targets'], 'targets')} 个达到 ≥20 singlet cells 的分析门槛；eligible 不等于 hit。</p>"
            f"{fig_card('fig03_target_coverage.png')}"
            f"<article class='panel'><h3>按 singlet cell 数分层（Target coverage distribution）</h3>{html_tables['coverage_bins']}</article>"
            f"<h3>覆盖最高的代表性 targets（Top coverage）</h3>{html_tables['top_coverage']}"
            "<p class='file-note'>完整 target 结果见交付表：<a href='../03_Tables/target_coverage_full.tsv'>完整 observed-target coverage</a> 与 <a href='../03_Tables/target_level_summary.tsv'>eligible target summary</a>。</p>"
        ),
        "effects": (
            f"<p class='lead'>Target-vs-NTC cell-level DE 以 {fmt(ntc['ntc_singlet_cells'], 'cells')} 个 NTC singlets 为 reference。"
            "Within-target BH FDR 指 target 内基因层面的 BH 校正 FDR；不能替代 biological-replicate inference。图表展示 perturbation-associated expression effects，不代表 mechanism validation。</p>"
            f"{fig_card('fig05_representative_perturbation_effects.png')}"
            "<h3>代表性 signed rank-biserial effects</h3><p>FDR 显著基因数为 0 的 target 如实保留；不据此宣称显著扰动。</p>"
            f"{fig_card('fig13_representative_de_effects.png')}{html_tables['representative_de_genes']}"
            f"<h3>代表性 target-level summary</h3>{html_tables['representatives']}"
            "<p class='file-note'>完整结果见 <a href='../03_Tables/differential_expression/target_level_DE_all.tsv.gz'>全量 gene-level DE</a> 与 <a href='../03_Tables/differential_expression/target_DE_summary.tsv'>target DE summary</a>。</p>"
        ),
        "robustness": (
            "<p class='lead'>Median guide-guide correlation、direction concordance、guide balance 与 LOOGO 是并列的稳健性证据；不计算单一综合分数。</p>"
            f"{fig_card('fig04_multi_guide_robustness.png')}<p>{esc(robustness_summary)}</p>"
            f"<p class='small-note'>{esc(sensitivity_text)} 主报告保留预设阈值与当前 evidence classes。</p>"
            "<p class='file-note'>完整 guide consistency 与 leave-one-guide-out 明细保留在交付 TSV。</p>"
        ),
        "cell_state": (
            "<div class='exploratory-banner'><strong>Exploratory</strong><span>Leiden clusters are unannotated grouping labels and do not establish cell states.</span></div>"
            "<h3>Local Leiden-cluster sensitivity</h3>"
            f"<p class='lead'>{esc(state_summary_text)} 在可评估 targets 中，cluster stratification 后总体 effect estimate 与 primary analysis 较一致，提示对当前 cluster composition 的敏感性相对有限；但这不证明 effect 一定可靠。</p>"
            "<p>Leiden clusters 是输入 H5AD 中未注释的分组标签，不能直接称为细胞状态或细胞类型。分层结果是单样本探索性敏感性分析，不证明细胞状态或扰动机制。</p>"
            "<p class='file-note'>完整 cluster composition 与 stratified sensitivity 仍见对应 TSV。</p>"
        ),
        "pathways": (
            "<div class='exploratory-banner'><strong>Exploratory</strong><span>Pathway significance does not independently establish a supported CRISPR hit or validated mechanism.</span></div>"
            f"<p class='lead'>Primary signed rank-biserial ranking is retained. Median tie fraction is {fmt(primary_tie.median(), 'fraction')} across {len(primary_tie):,} targets, limiting ranking resolution.</p>"
            f"<div class='metric-grid compact'>{ranking_cards}</div>{fig_card('fig06_pathway_ranking_qc.png')}{fig_card('fig07_reactome_ranking_sensitivity.png')}"
            "<p>Secondary ranking uses normalized mean-expression difference to assess direction sensitivity; it does not replace the primary ranking. The table shows only the existing representative pathways.</p>"
            f"<p>{esc(phkb_context)}</p>{html_tables['pathway']}"
            "<p class='notice compact-notice'>Reactome enrichment is exploratory for every target. For PARTIAL_SUPPORT and other non-SUPPORTED targets, pathway enrichment is a hypothesis-generating expression association, not a validated mechanism; pathway FDR cannot promote a target to SUPPORTED.</p>"
            "<p class='file-note'>完整 target-pathway Reactome 结果见交付文件：<a href='../03_Tables/pathway/reactome_GSEA_all.tsv.gz'>primary</a>、<a href='../03_Tables/pathway/reactome_GSEA_mean_difference.tsv.gz'>sensitivity</a> 与 <a href='../03_Tables/pathway/reactome_ranking_sensitivity.tsv.gz'>combined results</a>。</p>"
        ),
        "evidence": (
            f"<p class='lead'>最终判断遵循预先定义的 evidence classification。SUPPORTED={class_counts.get('SUPPORTED', 0)}；PARTIAL_SUPPORT={class_counts.get('PARTIAL_SUPPORT', 0)}。下表包含全部 PARTIAL_SUPPORT targets，并展示覆盖、代表性表达响应、guide 一致性、LOOGO、cluster-stratified sensitivity 与 pathway 线索。</p>"
            f"<div class='metric-grid compact'>{class_tier_cards}</div>"
            "<p class='small-note'>报告级 EXPLORATORY 汇总把原始细分类中除 SUPPORTED 与 PARTIAL_SUPPORT 外的类别合并展示；原始类别及规则保留如下，分类本身未重算或修改。</p>"
            f"<div class='legend'>{legend}</div>{html_tables['class_summary']}"
            f"<h3>代表性 target 证据总览（Representative target evidence）</h3>{html_tables['target_evidence']}"
            "<p class='file-note'>完整 observed-target 分类与 eligible-target 汇总仍保留于交付表；正文 target-level summary 已纳入全部 PARTIAL_SUPPORT targets。</p>"
        ),
        "limitations": (
            "<ul><li>本项目只有 1 个 biological sample；cell-level p 值/FDR 不能替代 biological replicate inference。</li>"
            f"<li>NTC reference 仅 {fmt(ntc['ntc_singlet_cells'], 'cells')} 个 singlet cells，target-vs-NTC DE、effect 与 pathway 结果需谨慎解释。</li>"
            "<li>Perturbation mode 为 UNKNOWN；target RNA decrease cannot prove editing success, and unchanged RNA cannot prove editing failure.</li>"
            "<li>Eligible targets 只表示达到当前分析覆盖阈值，不等于 CRISPR hits；guide imbalance 与 target coverage 会影响 robustness。</li>"
            "<li>Reactome pathway signals are exploratory；对 PARTIAL_SUPPORT 等低证据 targets 的富集结果不能解释为已验证机制。UMAP 使用 input H5AD 既有 embedding，仅作描述，UMAP proximity 不能证明 perturbation mechanism。</li>"
            "<li>当前结果不足以单独建立因果机制或独立验证某个靶点。</li></ul>"
        ),
        "reproducibility": (
            "<p class='lead'>正文保留客户解读所需的核心结果；完整表格、分析方法与复现配置按类别提供。</p>"
            f"{downloads}"
            "<p class='file-note'>ICM 完整分析摘要：<a href='../04_Methods_Reproducibility/icm/analysis_summary.pdf'>打开本地 PDF</a>。该文件由操作员从 ICM 导出，与正文内嵌的两张 ICM PNG 一并交付。</p>"
        ),
    }
    nav = "".join(f"<a href='#{key}'>{esc(label)}</a>" for key, label in SECTIONS)
    section_html = "".join(f"<section id='{key}' class='report-section'><h2>{esc(label)}</h2>{sections_data[key]}</section>" for key, label in SECTIONS)
    css = r"""
:root{--ink:#17212b;--muted:#5d6872;--line:#d8dee4;--soft:#f4f7f9;--accent:#3c5488;--teal:#00a087;--warning:#8a5a00}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#fff;color:var(--ink);font:15px/1.56 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}
.layout{display:grid;grid-template-columns:250px minmax(0,1fr);align-items:start}.sidebar{position:sticky;top:0;height:100vh;overflow:auto;padding:22px 14px;border-right:1px solid var(--line);background:#fbfcfd}.brand{font-size:14px;font-weight:750;letter-spacing:.04em;margin:0 8px 18px;color:#2e4055}.sidebar a{display:block;padding:8px 10px;margin:2px 0;border-radius:6px;text-decoration:none;color:#344252;font-size:13px}.sidebar a:hover,.sidebar a:focus{background:#eaf0f5;outline:2px solid transparent}
.main{max-width:1420px;width:100%;padding:30px clamp(18px,3.4vw,56px) 60px;margin:0 auto}.hero{border-bottom:1px solid var(--line);padding:6px 0 22px;margin-bottom:24px}.report-logo{display:block;width:min(260px,60vw);height:auto;max-height:58px;object-fit:contain;object-position:left;margin:0 0 12px}.eyebrow{text-transform:uppercase;letter-spacing:.09em;color:var(--muted);font-size:11px;font-weight:700}.hero h1{font-size:clamp(27px,3vw,40px);line-height:1.17;margin:8px 0}.hero p{margin:4px 0;color:var(--muted);font-size:13px}.report-section{padding:12px 0 30px;margin-bottom:20px;scroll-margin-top:18px;border-bottom:1px solid var(--line)}.report-section h2{font-size:24px;line-height:1.25;margin:8px 0 18px;color:#202c38}.report-section h3{font-size:17px;line-height:1.35;margin:18px 0 9px}.lead{font-size:16px;max-width:1000px;color:#394956}.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(140px,1fr));gap:10px;margin:16px 0}.metric-grid.compact{grid-template-columns:repeat(4,minmax(125px,1fr))}.metric-card,.panel{border:1px solid var(--line);border-radius:8px;background:#fff;padding:13px 15px}.metric-card span{display:block;color:var(--muted);font-size:12px}.metric-card strong{display:block;font-size:22px;line-height:1.24;margin:5px 0;color:#26394b}.metric-card small{display:block;color:var(--muted);font-size:11px}.conclusion{border-left:5px solid #00a087;background:#eff8f5;padding:15px 18px;margin:14px 0 18px}.conclusion strong{font-size:18px}.conclusion p{margin:6px 0 0}.finding{border-left:4px solid #3c5488;background:#f2f5fa;padding:12px 15px;margin:14px 0}.finding p{margin:5px 0}.legend{display:flex;gap:8px 15px;flex-wrap:wrap;margin:10px 0}.legend-item{font-size:12px;display:inline-flex;gap:6px;align-items:center}.legend-item i{width:11px;height:11px;border-radius:50%;display:inline-block}.class-chip{display:inline-block;border-left:4px solid var(--class-color);padding:1px 6px;background:#f7f8fa;border-radius:3px;font-size:11px;white-space:nowrap}
.two-col{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:12px 0}.notice{border-left:4px solid #d6a33a;background:#fff9ea;padding:12px 15px;margin:12px 0}.notice strong{display:block}.notice p{margin:4px 0 0}.compact-notice{font-size:13px}.exploratory-banner{display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap;padding:12px 15px;margin:5px 0 18px;border:1px solid #c9d7e5;background:#f5f8fb;border-radius:7px;color:#344b61}.exploratory-banner strong{letter-spacing:.02em}.file-note{font-size:13px;color:var(--muted)}.file-list{padding-left:20px}.file-list li{margin:5px 0}.small-note{font-size:12px;color:var(--muted)}.outcome-flow{margin:14px 0 20px;padding:13px;border:1px solid var(--line);border-radius:8px;background:#fbfcfd}.outcome-path{display:flex;align-items:center;justify-content:center;gap:14px}.outcome-step{min-width:180px;padding:10px 14px;border:1px solid #cbd5df;border-radius:6px;background:#fff;text-align:center}.outcome-step span,.outcome-classes span{display:block;font-size:12px;color:var(--muted)}.outcome-step strong,.outcome-classes strong{display:block;font-size:23px;color:#26394b}.outcome-outcomes{margin-top:12px;padding-top:10px;border-top:1px solid var(--line)}.outcome-label{display:block;text-align:center;font-size:12px;color:var(--muted);margin-bottom:8px}.outcome-classes{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.outcome-classes>div{text-align:center;padding:8px;border-radius:6px;background:#fff;border:1px solid var(--line)}.ntc-flow{display:flex;align-items:stretch;justify-content:space-between;gap:8px;margin:14px 0}.ntc-flow>b{align-self:center;color:#7a8792;font-size:18px}.flow-step{flex:1;min-width:0;padding:12px;border:1px solid var(--line);border-radius:7px;background:#fbfcfd;text-align:center}.flow-step span,.flow-step small{display:block;font-size:11px;color:var(--muted)}.flow-step strong{display:block;font-size:18px;line-height:1.3;margin:5px 0;overflow-wrap:anywhere}
.figure-card{margin:14px 0 18px;border:1px solid var(--line);border-radius:8px;padding:12px 14px;background:#fff}.figure-open{display:block;width:100%;border:0;background:#fff;padding:0;cursor:zoom-in}.figure-open img{display:block;width:100%;height:auto;max-height:720px;object-fit:contain}.figure-card h3{margin:10px 0 3px;font-size:15px}.figure-card p{color:var(--muted);font-size:12px;margin:3px 0}.figure-modal{display:none;position:fixed;inset:0;z-index:3000;background:rgba(11,20,30,.92);align-items:center;justify-content:center;padding:28px}.figure-modal.open{display:flex}.figure-modal img{max-width:96vw;max-height:90vh;object-fit:contain;background:#fff}.figure-modal button{position:absolute;right:22px;top:15px;background:#fff;border:0;border-radius:50%;width:38px;height:38px;font-size:24px;cursor:pointer}
.table-scroll{width:100%;overflow:auto;border:1px solid var(--line);border-radius:7px;margin:8px 0}.data-table{border-collapse:collapse;width:100%;font-size:12px;min-width:680px}.data-table th{background:#f1f4f6;text-align:left;color:#293a49;font-weight:700;position:sticky;top:0}.data-table th,.data-table td{padding:7px 9px;border-bottom:1px solid #e4e8ec;vertical-align:top}.data-table tbody tr:nth-child(even){background:#fafbfc}.data-table td{max-width:360px;overflow-wrap:anywhere}.empty{text-align:center;color:var(--muted);padding:18px!important}.table-controls{display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;margin:7px 0 14px;font-size:12px;color:var(--muted)}.table-controls input{font:inherit;padding:6px 8px;margin-left:6px;border:1px solid #bdc7d0;border-radius:5px;min-width:220px}.pager{display:flex;align-items:center;gap:8px}.pager button{font:inherit;background:#fff;border:1px solid #b9c4ce;border-radius:5px;padding:5px 9px;cursor:pointer;color:#26394b}.pager button:disabled{opacity:.45;cursor:default}
.footer{margin-top:35px;border-top:1px solid var(--line);padding-top:14px;color:var(--muted);font-size:12px}.bio-flow-module-head h3{margin-top:0}
@media(max-width:980px){.layout{grid-template-columns:1fr}.sidebar{position:relative;height:auto;display:flex;gap:4px;overflow:auto;border-right:0;border-bottom:1px solid var(--line);padding:8px}.brand{display:none}.sidebar a{white-space:nowrap}.main{padding-top:18px}.metric-grid{grid-template-columns:repeat(2,minmax(130px,1fr))}}
@media(max-width:640px){.two-col{grid-template-columns:1fr}.metric-grid,.metric-grid.compact{grid-template-columns:repeat(2,minmax(120px,1fr))}.metric-card{padding:10px}.metric-card strong{font-size:19px}.report-section h2{font-size:21px}.table-controls input{min-width:160px}.ntc-flow{display:grid;grid-template-columns:1fr}.ntc-flow>b{transform:rotate(90deg);justify-self:center}.outcome-path{flex-wrap:wrap}.outcome-step{min-width:135px}.outcome-classes{grid-template-columns:1fr}}
@media print{body{font-size:10pt}.layout{display:block}.sidebar,.table-controls,.figure-modal,.bio-flow-expand{display:none!important}.main{max-width:none;padding:0}.hero{padding-bottom:10px}.report-section{break-inside:auto;border-bottom:1px solid #aaa;padding:5px 0;margin:0}.report-section h2{font-size:17pt;break-after:avoid}.figure-card,.panel,.notice,.conclusion{break-inside:avoid;box-shadow:none}.figure-open{cursor:default}.figure-open img{max-height:none}.appendix{display:block;break-inside:auto}.appendix>summary{display:none}.appendix>*:not(summary){display:block!important}.appendix.print-skip{display:none!important}.table-scroll{overflow:visible;border:0}.data-table{display:table!important;min-width:0;font-size:7.5pt}.data-table thead{display:table-header-group}.data-table tbody{display:table-row-group}.data-table tr{display:table-row!important;break-inside:avoid}.data-table th,.data-table td{padding:3px 4px}.footer{margin-top:10px}}
"""
    js = r"""
(()=>{
  document.querySelectorAll('.paginated-table').forEach(table=>{
    const size=Number(table.dataset.pageSize||20), rows=Array.from(table.tBodies[0].rows), input=document.querySelector(`[data-filter-for="${table.id}"]`), pager=document.querySelector(`[data-pager-for="${table.id}"]`);
    if(!pager)return; const status=pager.querySelector('[data-status]'), prev=pager.querySelector('[data-prev]'), next=pager.querySelector('[data-next]'); let page=0, filtered=rows;
    const render=()=>{const pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.min(page,pages-1);rows.forEach(r=>r.style.display='none');filtered.slice(page*size,(page+1)*size).forEach(r=>r.style.display='table-row');status.textContent=`${filtered.length? page*size+1:0}–${Math.min((page+1)*size,filtered.length)} 条，共 ${filtered.length} 条`;prev.disabled=page===0;next.disabled=page>=pages-1;};
    const apply=()=>{const q=(input?.value||'').trim().toLowerCase();filtered=rows.filter(r=>r.textContent.toLowerCase().includes(q));page=0;render();};
    input?.addEventListener('input',apply);prev.addEventListener('click',()=>{page--;render();});next.addEventListener('click',()=>{page++;render();});render();
  });
  const modal=document.getElementById('figure-modal'), modalImg=modal?.querySelector('img');
  document.querySelectorAll('[data-open-figure]').forEach(btn=>btn.addEventListener('click',()=>{modalImg.src=btn.querySelector('img').src;modal.classList.add('open');}));
  const close=()=>{modal?.classList.remove('open');if(modalImg)modalImg.removeAttribute('src');};
  modal?.querySelector('button')?.addEventListener('click',close);modal?.addEventListener('click',e=>{if(e.target===modal)close();});document.addEventListener('keydown',e=>{if(e.key==='Escape')close();});const flowDetails=document.getElementById('read-flow-sankey');flowDetails?.addEventListener('toggle',()=>{if(flowDetails.open)requestAnimationFrame(()=>window.bioReadFlowRefresh?.());});
})();
"""
    logo_data = base64.b64encode((root / "02_Figures" / "epi-logo.png").read_bytes()).decode("ascii")
    title = f"scPerturb CRISPR 筛选分析报告 | {project.get('sample_id', '')}"
    page = ("<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{esc(title)}</title><style>{css}</style>{plotly_loader}</head><body>"
            f"<div class='layout'><nav class='sidebar' aria-label='Report navigation'><div class='brand'>scPERTURB · CRISPR</div>{nav}</nav>"
            f"<main class='main'><header class='hero'><img class='report-logo' src='data:image/png;base64,{logo_data}' alt='Epi logo'><div class='eyebrow'>单细胞 CRISPR 扰动筛选分析</div><h1>{esc(title)}</h1>"
            f"<p>项目 {esc(project.get('project_id'))} · 原生 DRAGEN {esc(project.get('dragen_version'))} 已完成 · 1 个生物学样本</p>"
            f"<p>观察靶点 {fmt(project.get('observed_targets'), 'target_count')} · 达到主要分析门槛 {fmt(project.get('eligible_targets'), 'target_count')} · 扰动模式 UNKNOWN</p></header>"
            f"{section_html}<footer class='footer'>scPerturb CRISPR 分析 · 独立离线报告 · 支持打印为 PDF</footer></main></div>"
            "<div class='figure-modal' id='figure-modal' role='dialog' aria-modal='true' aria-label='放大图'><button type='button' aria-label='关闭'>×</button><img alt='放大图'></div>"
            f"<script>{js}</script></body></html>")
    sample_prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", str(project.get("sample_id") or "sample")).strip("._-") or "sample"
    output = root / "01_Report" / f"{sample_prefix}_scPerturb_CRISPR_report.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    return output


def main() -> int:
    args = parse_args()
    report = build(args.root.resolve(), args.plotly_js.resolve())
    print(json.dumps({"report": str(report), "bytes": report.stat().st_size, "sections": len(SECTIONS)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
