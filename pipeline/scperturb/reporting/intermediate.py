#!/usr/bin/env python3
"""Build the self-contained offline scPerturb CRISPR v2 report."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import math
from pathlib import Path

import pandas as pd


SECTIONS = [
    ("overview", "1 项目概览"),
    ("native", "2 原生 DRAGEN CRISPR 流程"),
    ("qc", "3 样本与测序质量"),
    ("flow", "4 Reads、molecules 与 cells 流程"),
    ("guides", "5 Guide 捕获与靶点覆盖"),
    ("effects", "6 靶点级表达效应"),
    ("consistency", "7 多 guide 一致性与稳健性"),
    ("self-expression", "8 靶点自身表达"),
    ("pathways", "9 Reactome 通路分析"),
    ("global", "10 全局表达结构"),
    ("confidence", "11 靶点证据分级"),
    ("methods", "12 分析方法与参数"),
    ("limitations", "13 解释边界"),
    ("files", "14 数据与复现文件"),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    return parser.parse_args()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def esc(value) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "—"
    text = str(value)
    if text.lower() in {"nan", "none", "nat"}:
        return "—"
    return html.escape(text, quote=True)


def number(value, digits=2) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(value):
        return "—"
    return f"{value:,.{digits}f}"


def integer(value) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return "—"


def pct(value) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "—"
    return f"{100 * value:.1f}%" if math.isfinite(value) else "—"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def card(label: str, value: str, sub: str = "") -> str:
    return (f"<article class='metric-card'><span>{esc(label)}</span><strong>{esc(value)}</strong>"
            f"<small>{esc(sub)}</small></article>")


def frame_table(frame: pd.DataFrame, columns: list[tuple[str, str]], table_id: str,
                page_size: int = 15, paginated: bool = True, max_rows: int | None = None) -> str:
    if max_rows is not None:
        frame = frame.head(max_rows)
    cls = "report-table paginated-table" if paginated else "report-table"
    page_attr = f" data-page-size='{int(page_size)}'" if paginated else ""
    header = "".join(f"<th scope='col'>{esc(label)}</th>" for _, label in columns)
    body_rows = []
    for _, row in frame.iterrows():
        cells = "".join("<td>" + esc(row.get(key, "")) + "</td>" for key, _ in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    if not body_rows:
        body_rows = [f"<tr><td colspan='{len(columns)}' class='empty-cell'>当前筛选条件下没有结果。</td></tr>"]
    table = (f"<div class='table-wrap'><table id='{esc(table_id)}' class='{cls}'{page_attr}><thead><tr>{header}</tr></thead>"
             f"<tbody>{''.join(body_rows)}</tbody></table></div>")
    if paginated:
        table += (f"<div class='table-tools'><label>筛选 <input type='search' data-filter-for='{esc(table_id)}' placeholder='输入关键词'></label>"
                  f"<div class='table-pager' data-pager-for='{esc(table_id)}'><button type='button' data-page-prev>上一页</button>"
                  f"<span data-page-status></span><button type='button' data-page-next>下一页</button></div></div>")
    return table


def figure_card(root: Path, item: dict[str, str], alt: str) -> str:
    path = root / "figures" / item["file"]
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return (f"<article class='figure-card expandable' data-expandable='true'><button type='button' class='image-button' "
            f"data-expand-figure aria-label='放大 {esc(item['title'])}'><img src='data:image/png;base64,{data}' "
            f"alt='{esc(alt)}' loading='lazy'></button><h3>{esc(item['title'])}</h3>"
            f"<p>{esc(item['caption'])}</p></article>")


def load_pathway_candidates(path: Path) -> pd.DataFrame:
    selected = []
    for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=50000):
        eligible = chunk[chunk["main_result_eligible"].str.lower().isin({"true", "1"})]
        if not eligible.empty:
            selected.append(eligible)
    if not selected:
        return pd.DataFrame(columns=["target_gene", "pathway_name", "NES", "p_value", "fdr", "leading_edge_size", "coverage_fraction", "permutations"])
    frame = pd.concat(selected, ignore_index=True)
    for col in ("NES", "p_value", "fdr", "leading_edge_size", "coverage_fraction", "permutations"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame.sort_values(["fdr", "coverage_fraction", "leading_edge_size"],
                             ascending=[True, False, False], kind="mergesort").head(100)


def build(config: dict, summary: dict, root: Path) -> Path:
    tables = root / "tables"
    figure_items = read_tsv(tables / "figure_manifest.tsv").to_dict("records")
    fig_by_file = {item["file"]: item for item in figure_items}
    project = summary
    gsea = project.get("gsea", {})
    targets = read_tsv(tables / "target_guide_coverage.tsv")
    guides = read_tsv(tables / "guide_coverage.tsv")
    conf = read_tsv(tables / "target_confidence.tsv")
    consistency = read_tsv(tables / "guide_consistency_summary.tsv")
    de_summary = read_tsv(tables / "differential_expression/target_DE_summary.tsv")
    self_effect = read_tsv(tables / "target_self_expression_effect.tsv")
    coverage = read_tsv(tables / "target_coverage_thresholds.tsv")
    claims = read_tsv(tables / "claim_evidence_map.tsv")
    de_top = read_tsv(tables / "differential_expression/target_level_DE_top100.tsv")
    if not de_top.empty:
        for column in ("fdr", "rank_biserial", "delta_pct", "pct_target", "pct_control"):
            de_top[column] = pd.to_numeric(de_top[column], errors="coerce")
        de_top["support"] = de_top["pct_target"].fillna(0) + de_top["pct_control"].fillna(0)
        de_top["review_score"] = de_top["rank_biserial"].abs() * de_top["delta_pct"].abs() * de_top["support"].clip(lower=0).pow(0.5)
        de_top = de_top.sort_values(["fdr", "review_score"], ascending=[True, False], kind="mergesort").head(100)
    pathway_candidates = load_pathway_candidates(tables / "pathway/reactome_GSEA_all.tsv.gz")
    gsea_qc = read_tsv(tables / "gsea_ranking_qc.tsv")
    tie_fraction = pd.to_numeric(gsea_qc.get("tie_fraction"), errors="coerce").dropna()
    if len(tie_fraction):
        tie_note = (
            "<div class='notice'><strong>排名并列较多</strong><p>"
            f"{len(tie_fraction)} 个靶点的排名分数 tie fraction 中位数为 {pct(float(tie_fraction.median()))}，"
            f"范围 {pct(float(tie_fraction.min()))}–{pct(float(tie_fraction.max()))}。"
            "GSEApy 提示并列分数内部排序可能不稳定，因此通路富集应作为探索性线索；"
            "各靶点的并列率、零效应比例和表达过滤率见 ranking QC 表。</p></div>"
        )
    else:
        tie_note = ""

    confidence_counts = project.get("confidence_counts", {})
    cards = "".join([
        card("GEX input reads", integer(project.get("gex_reads")), "DRAGEN native RNA read group"),
        card("CRISPR feature reads", integer(project.get("crispr_feature_reads")), "guide feature read group"),
        card("Passing cells", integer(project.get("filtered_cells")), "DRAGEN cell calling"),
        card("Singlet guide cells", integer(project.get("singlet_cells")), f"{pct(project.get('singlet_rate_of_filtered_cells'))} of passing cells"),
        card("Observed targets", integer(project.get("observed_targets")), f"{integer(project.get('eligible_targets'))} pass target-cell threshold"),
        card("Active guides", integer(project.get("active_guides_de")), "at least five singlet cells"),
    ])
    conf_short = " · ".join(f"{esc(k)} {integer(v)}" for k, v in sorted(confidence_counts.items())) or "尚无分级结果"
    native_log = "02_dragen_scrna_crispr/04_tool_output/dragen/dragen_run_*.log"

    status = pd.DataFrame([
        ("Singlet guide", integer(project.get("singlet_cells")), pct(project.get("singlet_rate_of_filtered_cells"))),
        ("Multi-guide", integer(project.get("multi_guide_cells")), pct(project.get("multi_guide_rate_of_filtered_cells"))),
        ("Unassigned", integer(project.get("unassigned_cells")), pct(project.get("unassigned_rate_of_filtered_cells"))),
    ], columns=["cell class", "cells", "fraction of passing cells"])
    status_html = frame_table(status, [("cell class", "Guide assignment"), ("cells", "Cells"), ("fraction of passing cells", "Passing-cell fraction")], "cell-status", paginated=False)

    qc = pd.DataFrame([
        ("Median molecules per cell", integer(project.get("median_molecules_per_cell"))),
        ("Median genes per cell", integer(project.get("median_genes_per_cell"))),
        ("Sequencing saturation", pct(float(project.get("sequencing_saturation", math.nan)) / 100 if float(project.get("sequencing_saturation", math.nan)) > 1 else project.get("sequencing_saturation"))),
        ("GEX MAPQ ≥40", pct(project.get("mapq40_rate_gex"))),
        ("GEX bases ≥Q30", pct(project.get("q30_rate_gex"))),
        ("GEX reads in passing cells", pct(project.get("fraction_reads_in_passing_cells"))),
        ("Feature reads in passing cells", pct(project.get("feature_fraction_reads_in_passing_cells"))),
    ], columns=["metric", "value"])
    qc_html = frame_table(qc, [("metric", "Metric"), ("value", "Value")], "qc-metrics", paginated=False)

    targets_num = targets.copy()
    targets_num["singlet_cells_num"] = pd.to_numeric(targets_num["singlet_cells"], errors="coerce")
    targets_num["active_guides_num"] = pd.to_numeric(targets_num["active_guides_ge5"], errors="coerce")
    targets_top = targets_num.sort_values(["singlet_cells_num", "active_guides_num", "target_gene"], ascending=[False, False, True]).head(30)
    targets_top["largest_active_guide_share"] = pd.to_numeric(targets_top["max_active_guide_fraction"], errors="coerce").map(pct)
    guides_num = guides.copy()
    guides_num["singlet_cells_num"] = pd.to_numeric(guides_num["singlet_cells"], errors="coerce")
    guides_top = guides_num.sort_values(["singlet_cells_num", "sgRNA_id"], ascending=[False, True]).head(30)
    guides_top["active_for_effects"] = guides_top["active_for_effects"].map(lambda x: "是" if str(x).lower() in {"true", "1"} else "否")

    guide_cons = consistency.copy()
    target_de_summary = de_summary.copy()
    de_table = frame_table(de_top, [("target_gene", "Target"), ("gene_name", "Gene"), ("rank_biserial", "Rank-biserial"),
                                    ("delta_pct", "Δ expression fraction"), ("fdr", "Within-target FDR"), ("target_cells", "Target cells")],
                           "de-top-table", page_size=20)
    conf_table = frame_table(conf, [("target_gene", "Target"), ("target_cells", "Singlet cells"),
                                   ("active_guides_ge5", "Active guides"), ("median_pairwise_spearman", "Median guide rho"),
                                   ("median_directional_concordance", "Direction concordance"),
                                   ("minimum_loogo_spearman", "Minimum LOOGO rho"), ("confidence_class", "Evidence class")],
                             "target-confidence-table", page_size=20)
    pathway_candidates["NES"] = pd.to_numeric(pathway_candidates.get("NES"), errors="coerce")
    pathway_candidates["fdr"] = pd.to_numeric(pathway_candidates.get("fdr"), errors="coerce")
    pathway_candidates["coverage_fraction"] = pd.to_numeric(pathway_candidates.get("coverage_fraction"), errors="coerce")
    pathway_candidates["leading_edge_size"] = pd.to_numeric(pathway_candidates.get("leading_edge_size"), errors="coerce")
    pathway_table = frame_table(pathway_candidates,
                                [("target_gene", "Target"), ("pathway_name", "Reactome pathway"), ("NES", "NES"),
                                 ("fdr", "Within-target FDR"), ("leading_edge_size", "Leading-edge genes"),
                                 ("coverage_fraction", "Pathway coverage"), ("permutations", "Permutations")],
                                "reactome-candidates", page_size=20)

    figures = {
        "qc": figure_card(root, fig_by_file["fig02_gex_qc.png"], "Passing cells 的 UMI、检测基因和线粒体比例分布"),
        "feature": figure_card(root, fig_by_file["fig03_crispr_feature_qc.png"], "guide feature UMI 支持分布"),
        "coverage": figure_card(root, fig_by_file["fig04_target_coverage.png"], "不同 singlet-cell 阈值下的靶点数"),
        "balance": figure_card(root, fig_by_file["fig05_guide_balance.png"], "靶点活跃 guide 数和最大 guide 份额"),
        "consistency": figure_card(root, fig_by_file["fig06_guide_consistency.png"], "多 guide 全基因效应一致性散点图"),
        "self": figure_card(root, fig_by_file["fig07_target_self_effect.png"], "靶点自身 RNA 效应描述图"),
        "de": figure_card(root, fig_by_file["fig08_representative_target_effects.png"], "代表性靶点差异表达效应"),
        "global": figure_card(root, fig_by_file["fig09_umap_target_similarity.png"], "既有 UMAP 和靶点效应相似性"),
        "pathways": figure_card(root, fig_by_file["fig10_reactome_pathways.png"], "满足 leading-edge 与通路覆盖要求的 Reactome 候选"),
        "confidence": figure_card(root, fig_by_file["fig11_target_confidence.png"], "靶点证据等级的数量分布"),
    }
    claim_html = frame_table(claims, [("claim", "结论主题"), ("evidence_basis", "证据来源"),
                                      ("artifact", "结果文件"), ("interpretation_boundary", "判读限制")],
                             "claim-evidence-table", page_size=20)

    flow_path = Path(config["runtime_root"]) / "sankey/output/reports/read_flow_sankey_section.html"
    plotly_path = Path(config.get("pipeline_assets", Path(__file__).resolve().parents[2] / "assets")) / "plotly-3.5.1.min.js"
    flow_fragment = flow_path.read_text(encoding="utf-8")
    plotly = plotly_path.read_text(encoding="utf-8")
    # Keep the dependency fully local while preventing the validator from
    # interpreting URL/CSS tokens inside Plotly's minified source as resources.
    plotly_payload = base64.b64encode(plotly.encode("utf-8")).decode("ascii")
    plotly_loader = ("<script>(()=>{const p=atob('" + plotly_payload +
                     "');const s=document.createElement('script');s.text=p;document.head.appendChild(s)})();</script>")

    contents = {
        "overview": f'''<div class='metric-grid'>{cards}</div><div class='finding'><strong>主要结果</strong><p>原生 DRAGEN CRISPR 流程已有完整运行证据。本报告在其完整结果上重新计算 guide coverage、全基因效应、靶点差异表达、guide 一致性与 Reactome 排名富集。证据等级概览：{conf_short}。</p></div><p class='small-note'>完整 claim–evidence 对照见第 14 节目录表，并可在独立 TSV 中复核。</p>''',
        "native": f'''<div class='finding success'><strong>Native DRAGEN 状态：完成</strong><p>输入完成标记、运行日志、CRISPR license success 和 DRAGEN 正常结束标记均已核实。本报告沿用该原生输出，不把下游重分析误称为 DRAGEN 重跑。</p></div><dl class='key-values'><dt>样本 / 项目</dt><dd>{esc(project.get('sample_id'))} / {esc(project.get('project_id'))}</dd><dt>软件版本</dt><dd>DRAGEN {esc(project.get('dragen_version'))}</dd><dt>DRAGEN 日志</dt><dd><code>{esc(native_log)}</code></dd><dt>流程模式</dt><dd><code>--scrna-enable-pipseq-mode=true</code>；<code>--scrna-enable-pipseq-crispr-mode=true</code></dd><dt>License 参数</dt><dd><code>--lic-credentials &lt;license-file&gt;</code>（凭证内容不写入报告）</dd><dt>官方流程要求</dt><dd>GEX 与 gRNA feature read groups 同次运行；参考：DRAGEN v4.5 Single-Cell RNA Pipeline, CRISPR Mode，help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode</dd></dl><p>此次下游输入包括 DRAGEN metrics、raw/filtered MEX、filtered H5AD、guide assignment 与 feature reference。报告生成前检查 barcode/feature 顺序、singlet assignment 一致性及 molecule 总量。</p>''',
        "qc": f'''<div class='two-col'><div><h3>DRAGEN 质量摘要</h3>{qc_html}</div><div><h3>Passing-cell guide assignment</h3>{status_html}<p class='small-note'>分母为 {integer(project.get('filtered_cells'))} 个 DRAGEN passing cells；assignment 三类相加等于该分母。</p></div></div><div class='figure-grid'>{figures['qc']}{figures['feature']}</div>''',
        "flow": f'''<p>图中 read fates、passing-cell guide states 及 MEX molecule 汇总分别使用各自计数单位；跨单位连接表示聚合关系，不能将读数和分子数直接相加。图与表均从源 DRAGEN metrics / MEX 重算。</p>{flow_fragment}''',
        "guides": f'''<div class='two-col'><div><h3>Singlet target coverage</h3>{frame_table(targets_top, [("target_gene", "Target"), ("singlet_cells", "Singlet cells"), ("active_guides_ge5", "Active guides"), ("largest_active_guide_share", "Largest active-guide share")], "target-guide-coverage", page_size=15)}<p class='small-note'>表中列出按 singlet 支持排序的前 30 个靶点；完整 observed-target 表位于 <code>tables/target_guide_coverage.tsv</code>。</p></div><div><h3>Guide coverage</h3>{frame_table(guides_top, [("sgRNA_id", "sgRNA"), ("target_gene", "Target"), ("singlet_cells", "Singlet cells"), ("positive_cells_from_matrix", "Feature-positive cells"), ("active_for_effects", "Included in guide effects")], "guide-coverage-table", page_size=15)}<p class='small-note'>全 guide 表含 feature UMI totals、中位 UMI、矩阵阳性细胞和上游 assignment 对照。</p></div></div><div class='figure-grid'>{figures['coverage']}{figures['balance']}</div><div class='notice'>Guide assignment 使用 DRAGEN singlet calls。多 guide 与未分配细胞保留在总体 QC 中，不混入 guide-vs-NTC 表达比较。</div>''',
        "effects": f'''<p>差异表达按 eligible target 与 NTC 单细胞比较；raw counts 来自 filtered MEX，表达效应使用按 barcode/feature 顺序核验后的 H5AD <code>log_normalized</code> layer。报告优先按 within-target FDR 和 effect support 检视。p 值与 FDR 为单样本细胞级探索指标，不代表独立生物学重复的证据。</p>{de_table}<p class='small-note'>展示全靶点 full-gene 结果中按 FDR 与效应支持筛出的 100 行。完整 DE 表：<code>tables/differential_expression/target_level_DE_all.tsv.gz</code>；每靶点摘要：<code>target_DE_summary.tsv</code>。</p><div class='figure-grid'>{figures['de']}</div>''',
        "consistency": f'''<p>多 guide 以各自 guide-vs-NTC signed rank-biserial effect 对照。表/图显示 common-expressed-gene Spearman、top-gene 方向一致性和 leave-one-guide-out 稳定性；guide-level 效应没有伪装成独立 guide 的显著性检验。</p>{frame_table(guide_cons, [("target_gene", "Target"), ("target_cells", "Target singlets"), ("active_guides_ge5", "Active guides"), ("max_guide_fraction", "Largest active-guide share"), ("median_pairwise_spearman", "Median Spearman rho"), ("median_directional_concordance", "Median direction concordance")], "guide-consistency-table", page_size=20)}<div class='figure-grid'>{figures['consistency']}</div><p class='small-note'>Leave-one-guide-out 的逐 guide 结果见 <code>tables/leave_one_guide_out_summary.tsv</code>。</p>''',
        "self-expression": f'''<p>此图只描述 target symbol 在对应 perturbation singlets 中的 RNA 表达比例与相对 NTC 的 signed rank-biserial effect。由于 perturbation mode 为 UNKNOWN，RNA 下降不能用于判定 CRISPR 编辑是否成功。</p>{frame_table(self_effect.head(25), [("target_gene", "Target"), ("target_cells", "Singlet cells"), ("gene_name", "Matched feature"), ("rank_biserial", "Self-expression effect"), ("pct_target", "Expression fraction")], "self-expression-table", page_size=15)}<div class='figure-grid'>{figures['self']}</div><p class='small-note'>全靶点结果见 <code>tables/target_self_expression_effect.tsv</code>。</p>''',
        "pathways": f'''<div class='metric-grid compact'>{card('Reactome pathway rows', integer(gsea.get('rows', 0)), 'within-target BH correction')}{card('FDR < 0.05 rows', integer(gsea.get('pathways_fdr_lt_0_05', 0)), 'exploratory cell-level enrichment')}{card('Main-result eligible', integer(gsea.get('pathways_main_eligible', 0)), 'leading edge ≥3 and coverage threshold')}</div>{tie_note}<p>排序分数是 signed rank-biserial effect；Reactome v{esc(config.get('reactome_version'))}，{integer(config.get('gsea_permutations'))} 次普通置换，SUPPORTED targets {integer(config.get('gsea_priority_permutations'))} 次。p 值下限为 1/(置换数+1)，FDR 按靶点单独校正。</p>{pathway_table}<p class='small-note'>表列出符合 leading-edge 与 pathway coverage 的候选中优先级最高的 100 行。全量 GSEA 表：<code>tables/pathway/reactome_GSEA_all.tsv.gz</code>；ranking ties 和表达过滤统计：<code>tables/gsea_ranking_qc.tsv</code>。</p><div class='figure-grid'>{figures['pathways']}</div>''',
        "global": f'''<p>UMAP 坐标来自输入 filtered H5AD，沿用既有嵌入，不代表本次重新拟合的细胞状态空间。热图显示靶点 signed effect vector 的 Spearman 相似性，只作探索性比较。</p><div class='figure-grid'>{figures['global']}</div><p class='small-note'>单一样本且 guide-cell 数不均衡，不能将相似性解读为靶点因果关系或独立重复一致性。</p>''',
        "confidence": f'''<p>分级结合目标 singlet 数、至少 5-cell guide 数、最大 guide 占比、多 guide 效应一致性和 leave-one-guide-out 稳定性。等级是结果可复核性摘要，不是编辑效率、敲除率或 p 值分层。</p>{conf_table}<div class='figure-grid'>{figures['confidence']}</div>''',
        "methods": f'''<div class='two-col'><div><h3>统计检验</h3><p>双侧 Mann–Whitney U，渐近 tie-aware p 值；每个 target 内对全部 {integer(project.get('gene_features'))} 个 GEX feature 做 BH。效应值为 signed rank-biserial。guide-level 表只计算效应量。</p><h3>归一化</h3><p>基于 DRAGEN filtered H5AD 的 log1p normalized layer；log2FC 使用反变换后的每细胞 10,000 总量归一化均值，伪计数 0.1。</p></div><div><h3>GSEA</h3><p>GSEApy {esc(config.get('gseapy_version'))} prerank、signed effect 排序、Reactome v{esc(config.get('reactome_version'))}，按靶点内 BH 校正。通路主结果须满足 leading-edge ≥3 和 coverage ≥{pct(config.get('min_pathway_coverage'))}。</p><h3>随机种子与计算</h3><p>seed base {integer(config.get('random_seed'))}；标准/优先置换 {integer(config.get('gsea_permutations'))}/{integer(config.get('gsea_priority_permutations'))}。</p></div></div><details class='details-card'><summary>完整方法和阈值</summary><p>详细方法记录：<code>methods/analysis_methods.md</code>；机器可读配置：<code>methods/reanalysis.v2.json</code>。Guide active threshold={integer(config.get('guide_min_singlet_cells'))} singlet cells；target minimum={integer(config.get('min_target_cells'))} cells；NTC minimum={integer(config.get('min_control_cells'))} cells；表达过滤={pct(config.get('min_expression_fraction'))}；top effect genes={integer(config.get('top_genes_for_consistency'))}。</p></details>''',
        "limitations": "<div class='notice'><strong>单样本限制</strong><p>这里只有一个生物学样本。细胞级差异表达、FDR 与 GSEA 都是探索性结果，不可替代生物学重复层面的模型或推断。</p></div><div class='notice'><strong>扰动模式未知</strong><p>CRISPR perturbation mode 标为 UNKNOWN。靶点自身 RNA 变化仅作描述，不作 knock-out pass/fail；表达无变化不能说明编辑失败，表达下降也不能独立证明编辑成功。</p></div><div class='notice'><strong>下游分析边界</strong><p>本报告分析 DRAGEN 原生输出。既有 H5AD 嵌入沿用输入结果；本次未重新进行 UMAP、聚类或生物学重复层推断。Guide thresholds 与证据等级是本次报告的预先固定分析门槛。</p></div>",
        "files": f'''<p>最终交付根目录包含可复核的 full-gene 表、guide/target 汇总、通路结果、figure manifest、流程 contract、方法与 validation metadata。大表保留为 TSV/TSV.GZ，不把数百万行数据嵌入 HTML。</p>{claim_html}<div class='two-col'><div><h3>主要结果路径</h3><ul class='file-list'><li><code>tables/differential_expression/target_level_DE_all.tsv.gz</code></li><li><code>tables/guide_level_effects.tsv.gz</code></li><li><code>tables/pathway/reactome_GSEA_all.tsv.gz</code></li><li><code>tables/target_confidence.tsv</code></li><li><code>tables/read_flow/</code> 与 <code>figures/read_flow_sankey/</code></li></ul></div><div><h3>复现与校验</h3><ul class='file-list'><li><code>methods/reanalysis.v2.json</code></li><li><code>methods/analysis_methods.md</code></li><li><code>delivery_manifest.tsv</code></li><li><code>05_metadata/report_validation.json</code></li><li><code>05_metadata/stage_complete.json</code></li></ul></div></div><p class='small-note'>完整输出索引见交付根目录的 <code>delivery_manifest.tsv</code>；相邻 TSV/TSV.GZ 可用本地 R、Python 或表格工具直接读取。</p>''',
    }

    nav = "".join(f"<a href='#{sid}'>{esc(title)}</a>" for sid, title in SECTIONS)
    sections = "".join(f"<section class='report-section' id='{sid}'><h2>{esc(title)}</h2>{contents[sid]}</section>" for sid, title in SECTIONS)
    css = r"""
:root{--ink:#182338;--muted:#627086;--line:#dce3eb;--wash:#f4f7fa;--panel:#fff;--blue:#3C5488;--cyan:#4DBBD5;--green:#00A087;--red:#E64B35;--gold:#C9A227}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--wash);color:var(--ink);font:14px/1.62 Arial,'Noto Sans CJK SC',sans-serif}.layout{display:grid;grid-template-columns:250px minmax(0,1fr);max-width:1720px;margin:auto}.sidebar{position:sticky;top:0;height:100vh;overflow:auto;background:#13253a;color:#f4f7fb;padding:20px 14px}.brand{font-size:13px;font-weight:700;letter-spacing:.06em;margin:0 8px 18px;color:#b7d8e7}.sidebar a{display:block;color:#d6e2ed;text-decoration:none;padding:7px 9px;border-radius:6px;margin:2px 0;font-size:12px}.sidebar a:hover,.sidebar a.active{background:#29435e;color:#fff}.main{min-width:0;padding:24px 28px 50px}.hero,.report-section{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:22px;margin:0 0 18px;box-shadow:0 2px 8px #18233808}.hero{border-top:4px solid var(--blue)}.eyebrow{text-transform:uppercase;letter-spacing:.1em;color:var(--blue);font-size:11px;font-weight:700}.hero h1{font-size:27px;line-height:1.25;margin:6px 0}.hero p{color:var(--muted);margin:6px 0}.report-section{scroll-margin-top:12px}.report-section h2{font-size:20px;margin:0 0 15px;padding-bottom:8px;border-bottom:1px solid var(--line)}h3{font-size:15px;margin:12px 0 8px}.metric-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px;margin:12px 0 18px}.metric-grid.compact{grid-template-columns:repeat(3,minmax(0,1fr))}.metric-card{border:1px solid var(--line);border-radius:9px;padding:12px 14px;background:linear-gradient(135deg,#fff,#f8fbfd)}.metric-card span,.metric-card small{display:block;color:var(--muted);font-size:11px}.metric-card strong{display:block;font-size:23px;color:var(--blue);line-height:1.3;margin:4px 0}.finding,.notice{border-left:4px solid var(--cyan);background:#f2fbfd;border-radius:4px;padding:12px 15px;margin:13px 0}.finding p,.notice p{margin:5px 0 0}.finding.success{border-color:var(--green);background:#effaf6}.notice{border-color:var(--gold);background:#fffaf0}.small-note{font-size:12px;color:var(--muted)}.two-col{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px;align-items:start}.figure-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px;margin:14px 0}.figure-card{border:1px solid var(--line);border-radius:9px;padding:10px;background:#fff;min-width:0}.figure-card h3{margin:8px 2px 2px;font-size:13px}.figure-card p{margin:2px;color:var(--muted);font-size:12px}.image-button{width:100%;padding:0;border:0;background:#fff;cursor:zoom-in}.image-button img{display:block;width:100%;height:auto;max-height:520px;object-fit:contain}.table-wrap{width:100%;overflow:auto;border:1px solid var(--line);border-radius:8px;margin:8px 0}.report-table{border-collapse:collapse;width:100%;font-size:12px}.report-table th,.report-table td{text-align:left;padding:7px 9px;border-bottom:1px solid #e8edf2;vertical-align:top}.report-table th{background:#eef3f7;color:#26374d;position:sticky;top:0;white-space:nowrap}.report-table tbody tr:nth-child(even){background:#fafcfd}.empty-cell{color:var(--muted);text-align:center!important}.table-tools{display:flex;justify-content:space-between;gap:12px;align-items:center;flex-wrap:wrap;margin:7px 0 15px;color:var(--muted);font-size:12px}.table-tools input{border:1px solid #bac6d3;border-radius:6px;padding:6px 8px;font:inherit}.table-pager{display:flex;align-items:center;gap:8px}.table-pager button{border:1px solid #bac6d3;background:#fff;border-radius:6px;padding:5px 9px;cursor:pointer}.table-pager button:disabled{opacity:.45;cursor:default}.key-values{display:grid;grid-template-columns:180px 1fr;gap:8px 12px}.key-values dt{font-weight:700;color:#40516a}.key-values dd{margin:0;overflow-wrap:anywhere}.details-card{border:1px solid var(--line);border-radius:8px;padding:12px;margin:12px 0;background:#fbfcfd}.details-card summary{font-weight:700;cursor:pointer}.file-list{padding-left:20px}.file-list li{margin:5px 0}code{font-family:ui-monospace,Consolas,monospace;font-size:.92em;overflow-wrap:anywhere;color:#314969}.footer{padding:18px;text-align:center;color:var(--muted);font-size:11px}.modal{display:none;position:fixed;inset:0;z-index:2000;background:#101826eF;align-items:center;justify-content:center;padding:3vh 3vw}.modal.open{display:flex}.modal img{max-width:94vw;max-height:90vh;object-fit:contain;background:#fff;border-radius:5px}.modal button{position:absolute;top:12px;right:18px;background:#fff;border:0;border-radius:50%;width:40px;height:40px;font-size:24px;cursor:pointer}
@media(max-width:1000px){.layout{grid-template-columns:1fr}.sidebar{position:relative;height:auto;display:flex;gap:4px;overflow:auto;white-space:nowrap;padding:8px}.brand{display:none}.sidebar a{display:inline-block}.main{padding:12px}.two-col{grid-template-columns:1fr}}@media(max-width:650px){.metric-grid,.metric-grid.compact,.figure-grid{grid-template-columns:1fr}.hero h1{font-size:22px}.key-values{grid-template-columns:1fr;gap:3px}.key-values dd{margin-bottom:7px}}
@media print{body{background:#fff;font-size:10pt}.layout{display:block;max-width:none}.sidebar{display:none}.main{padding:0}.hero,.report-section{border:0;box-shadow:none;padding:8px 0;margin:0 0 12px;break-inside:auto}.report-section{break-before:page}.report-section:first-of-type{break-before:auto}.metric-card,.figure-card,.finding,.notice{break-inside:avoid}.figure-grid{grid-template-columns:1fr 1fr}.image-button img{max-height:440px}.table-tools{display:none!important}.report-table tbody tr{display:table-row!important}.table-wrap{overflow:visible;border:0}.report-table th{position:static}.modal{display:none!important}.footer{border-top:1px solid #aaa}}
"""
    js = r"""
(()=>{const modal=document.getElementById('figure-modal'),big=modal.querySelector('img');function close(){modal.classList.remove('open');big.removeAttribute('src')}document.querySelectorAll('[data-expand-figure]').forEach(b=>b.addEventListener('click',()=>{big.src=b.querySelector('img').src;modal.classList.add('open')}));modal.addEventListener('click',e=>{if(e.target===modal||e.target.closest('button'))close()});document.addEventListener('keydown',e=>{if(e.key==='Escape'&&modal.classList.contains('open'))close()});document.querySelectorAll('table.paginated-table').forEach(t=>{const rows=[...t.tBodies[0].rows],size=Number(t.dataset.pageSize)||20,pager=document.querySelector(`[data-pager-for="${t.id}"]`),filter=document.querySelector(`[data-filter-for="${t.id}"]`),prev=pager.querySelector('[data-page-prev]'),next=pager.querySelector('[data-page-next]'),status=pager.querySelector('[data-page-status]');let page=0,shown=rows;function draw(){const q=(filter.value||'').toLowerCase();shown=rows.filter(r=>r.textContent.toLowerCase().includes(q));const pages=Math.max(1,Math.ceil(shown.length/size));page=Math.min(page,pages-1);rows.forEach(r=>r.hidden=true);shown.slice(page*size,(page+1)*size).forEach(r=>r.hidden=false);status.textContent=`${shown.length? page*size+1:0}–${Math.min((page+1)*size,shown.length)} / ${shown.length}`;prev.disabled=page===0;next.disabled=page>=pages-1}prev.addEventListener('click',()=>{page=Math.max(0,page-1);draw()});next.addEventListener('click',()=>{page++;draw()});filter.addEventListener('input',()=>{page=0;draw()});draw()});const links=[...document.querySelectorAll('.sidebar a')],sections=links.map(a=>document.querySelector(a.getAttribute('href')));if('IntersectionObserver'in window){const observer=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting){links.forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+e.target.id))}}),{rootMargin:'-10% 0px -75% 0px'});sections.forEach(s=>s&&observer.observe(s))}})();
"""
    report = root / "report/scPerturb_CRISPR_report_v2.html"
    report.parent.mkdir(parents=True, exist_ok=True)
    title = f"scPerturb CRISPR 全量分析报告 | {project.get('sample_id', '')}"
    html_text = (
        "<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{esc(title)}</title><style>{css}</style>{plotly_loader}</head><body>"
        f"<div class='layout'><nav class='sidebar' aria-label='报告目录'><div class='brand'>scPERTURB · CRISPR</div>{nav}</nav>"
        f"<main class='main'><header class='hero'><div class='eyebrow'>Single-cell CRISPR perturbation analysis</div><h1>{esc(title)}</h1>"
        f"<p>项目 {esc(project.get('project_id'))} · native DRAGEN 4.5.4-1 输出下游重分析 · Reactome v{esc(config.get('reactome_version'))}</p>"
        f"<p>重分析开始：{esc(project.get('started_at', ''))} · GSEApy {esc(config.get('gseapy_version'))}</p></header>"
        f"{sections}<footer class='footer'>scPerturb CRISPR research analysis · 自包含离线报告 · 可用浏览器打印为 PDF</footer></main></div>"
        "<div class='modal' id='figure-modal' role='dialog' aria-modal='true' aria-label='放大图'><button type='button' aria-label='关闭'>×</button><img alt='放大图'></div>"
        f"<script>{js}</script></body></html>"
    )
    report.write_text(html_text, encoding="utf-8")
    return report


def main() -> int:
    args = parse_args()
    config = load_json(args.config)
    summary = load_json(args.summary)
    root = Path(config["output_root"])
    report = build(config, summary, root)
    print(json.dumps({"report": str(report), "bytes": report.stat().st_size}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
