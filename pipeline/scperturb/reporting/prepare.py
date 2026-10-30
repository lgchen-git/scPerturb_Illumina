#!/usr/bin/env python3
"""Prepare the evidence-led scPerturb CRISPR V3 result package from audited V2 outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import mannwhitneyu, spearmanr


SAMPLE_ID = ""
H5AD = Path()
FEATURE_REFERENCE = Path()
FILTERED_FEATURE_MATRIX = Path()
FILTERED_FEATURES = Path()
FILTERED_BARCODES = Path()
REACTOME_GMT = Path()
CONTROL = "Non-Targeting Control"

from scperturb.reporting.style import NPG_CATEGORICAL, NPG_EVIDENCE, apply_matplotlib_style

NPG = list(NPG_CATEGORICAL)
CLASS_COLORS = dict(NPG_EVIDENCE)
CLASS_ORDER = ["SUPPORTED", "PARTIAL_SUPPORT", "DISCORDANT", "GUIDE_DOMINATED", "LOW_GUIDE_SUPPORT", "LOW_CELL", "INCONCLUSIVE"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--control-gene", default="Non-Targeting Control")
    parser.add_argument("--v2-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--runtime-root", type=Path, required=True)
    parser.add_argument("--reactome-gmt", type=Path, required=True)
    parser.add_argument("--skip-pathway-sensitivity", action="store_true", help="Only for smoke tests; never use for delivery.")
    parser.add_argument("--reuse-pathway-sensitivity", action="store_true", help="Resume final packaging from a complete secondary-ranking output already in output-root.")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False, default=json_scalar) + "\n", encoding="utf-8")


def json_scalar(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_frame(path: Path, frame: pd.DataFrame, gzip_output: bool | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    use_gzip = path.suffix == ".gz" if gzip_output is None else gzip_output
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", compression="gzip" if use_gzip else None)


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle]


def as_float(value) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return math.nan
    return value if math.isfinite(value) else math.nan


def fmt_count(value) -> str:
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError, OverflowError):
        return "—"


def pct(n: int | float, d: int | float) -> float:
    return float(n) / float(d) if d else math.nan


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.clip(q, 0, 1)
    return result


def gmt_maps(path: Path) -> tuple[dict[str, set[str]], dict[str, str], dict[str, str]]:
    sets: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    name_to_id: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            # Reactome v97 GMT: pathway name, Reactome stable identifier, genes...
            name, stable_id = fields[0], fields[1]
            names[stable_id] = name
            name_to_id[name] = stable_id
            sets[stable_id] = set(fields[2:])
    if not sets:
        raise RuntimeError(f"No pathway sets found in {path}")
    return sets, names, name_to_id


def update_primary_reactome_labels(path: Path, name_to_id: dict[str, str], names: dict[str, str]) -> int:
    tmp = path.with_suffix(path.suffix + ".v3tmp")
    corrected = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as src, gzip.open(tmp, "wt", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        columns = list(reader.fieldnames or [])
        writer = csv.DictWriter(dst, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            old_id, old_name = row.get("pathway_id", ""), row.get("pathway_name", "")
            # V2 stored GMT field 1 as pathway_id and field 2 as pathway_name; reverse only labels.
            candidate_id = old_name if old_name in names else name_to_id.get(old_id, "")
            candidate_name = names.get(candidate_id, old_id)
            if candidate_id:
                row["pathway_id"] = candidate_id
                row["pathway_name"] = candidate_name
                corrected += int(candidate_id != old_id or candidate_name != old_name)
            writer.writerow(row)
    tmp.replace(path)
    return corrected


def classify_row(row: dict, thresholds: dict) -> str:
    n = int(float(row["target_cells"]))
    n_guides = int(float(row["active_guides_ge5"]))
    max_share = as_float(row["max_guide_fraction"])
    rho = as_float(row["median_pairwise_spearman"])
    direction = as_float(row["median_directional_concordance"])
    loo_rho = as_float(row["minimum_loogo_spearman"])
    loo_direction = as_float(row["minimum_loogo_direction_retention"])
    if n < thresholds["min_target_cells"]:
        return "LOW_CELL"
    if n_guides < thresholds["min_active_guides"]:
        return "LOW_GUIDE_SUPPORT"
    if math.isfinite(max_share) and max_share > thresholds["max_guide_fraction"]:
        return "GUIDE_DOMINATED"
    if (math.isfinite(rho) and rho < 0) or (math.isfinite(direction) and direction < thresholds["direction_discordant_below"]):
        return "DISCORDANT"
    loo_stable = (math.isfinite(loo_rho) and math.isfinite(loo_direction)
                  and loo_rho >= thresholds["loogo_spearman_min"]
                  and loo_direction >= thresholds["loogo_direction_retention_min"])
    if (math.isfinite(rho) and rho >= thresholds["spearman_supported_at_least"]
            and math.isfinite(direction) and direction >= thresholds["direction_supported_at_least"] and loo_stable):
        return "SUPPORTED"
    if (math.isfinite(rho) and rho >= 0 and math.isfinite(direction)
            and direction >= thresholds["direction_discordant_below"]):
        return "PARTIAL_SUPPORT"
    return "INCONCLUSIVE"


def threshold_sensitivity(conf: pd.DataFrame, primary_classes: dict[str, str], baseline: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    grids = {
        "min_target_cells": [10, 20, 30],
        "min_active_guides": [1, 2, 3],
        "max_guide_fraction": [0.4, 0.5, 0.6],
        "spearman_supported_at_least": [0.2, 0.3, 0.4],
        "direction_discordant_below": [0.4, 0.5, 0.6],
        "direction_supported_at_least": [0.7, 0.8, 0.9],
        "loogo_spearman_min": [0.7, 0.8, 0.9],
        "loogo_direction_retention_min": [0.7, 0.8, 0.9],
    }
    rows, changes = [], []
    for parameter, values in grids.items():
        for value in values:
            thresholds = dict(baseline)
            thresholds[parameter] = value
            scenario = f"{parameter}={value:g}"
            results = {}
            for record in conf.to_dict("records"):
                target = record["target_gene"]
                cls = classify_row(record, thresholds)
                results[target] = cls
                if cls != primary_classes[target]:
                    changes.append((scenario, parameter, value, target, primary_classes[target], cls))
            counts = Counter(results.values())
            rows.append({"scenario": scenario, "parameter": parameter, "threshold": value,
                         "targets": len(results), "changed_from_primary": sum(results[t] != primary_classes[t] for t in results),
                         **{cls: counts.get(cls, 0) for cls in CLASS_ORDER}})
    return pd.DataFrame(rows), pd.DataFrame(changes, columns=["scenario", "parameter", "threshold", "target_gene", "primary_class", "sensitivity_class"])


def iter_target_groups(path: Path, chunksize: int = 180000):
    usecols = ["target_gene", "target_cells", "control_cells", "gene_id", "gene_name", "mean_target_normalized",
               "mean_control_normalized", "delta_pct", "rank_biserial", "pct_target", "pct_control", "p_value", "fdr"]
    carry = None
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, chunksize=chunksize,
                             dtype={"target_gene": str, "gene_id": str, "gene_name": str}):
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        if chunk.empty:
            continue
        last_target = chunk.iloc[-1]["target_gene"]
        for target, group in chunk.groupby("target_gene", sort=False):
            if target == last_target:
                carry = group.copy()
            else:
                yield str(target), group
        if last_target not in set(chunk["target_gene"].iloc[:-1].tolist()):
            continue
    if carry is not None and not carry.empty:
        yield str(carry.iloc[0]["target_gene"]), carry


def run_secondary_gsea(de_path: Path, out_path: Path, qc_path: Path, config: dict, confidence: pd.DataFrame,
                       targets: pd.DataFrame, pathways: dict[str, set[str]], pathway_names: dict[str, str]) -> int:
    expected = str(config.get("gseapy_version", ""))
    if expected and gseapy.__version__ != expected:
        raise RuntimeError(f"GSEApy version mismatch: expected {expected}, found {gseapy.__version__}")
    min_pct = float(config["min_expression_fraction"])
    min_size, max_size = int(config["min_gene_set_size"]), int(config["max_gene_set_size"])
    permutations = int(config["gsea_permutations"])
    cell_counts = dict(zip(targets["target_gene"], pd.to_numeric(targets["singlet_cells"], errors="coerce").astype(int)))
    eligible = set(targets.loc[pd.to_numeric(targets["singlet_cells"], errors="coerce") >= int(config["min_target_cells"]), "target_gene"])
    qc_primary = read_tsv(Path(config["v2_root"]) / "tables/gsea_ranking_qc.tsv").set_index("target_gene")
    write_rows = []
    qc_rows = []
    count = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["target_gene", "target_cells", "control_cells", "pathway_id", "pathway_name", "ES", "NES", "p_value", "fdr",
                         "fwer_p_value", "leading_edge", "leading_edge_size", "ranked_unique_symbols", "rank_unique_scores",
                         "pathway_gene_count", "genes_in_ranked_universe", "coverage_fraction", "permutations", "priority_class",
                         "main_result_eligible"])
        for i, (target, frame) in enumerate(iter_target_groups(de_path), start=1):
            if target not in eligible:
                continue
            expressed = ((pd.to_numeric(frame["pct_target"], errors="coerce") >= min_pct)
                         | (pd.to_numeric(frame["pct_control"], errors="coerce") >= min_pct))
            score = (pd.to_numeric(frame["mean_target_normalized"], errors="coerce")
                     - pd.to_numeric(frame["mean_control_normalized"], errors="coerce"))
            ranked = frame.loc[expressed & score.notna(), ["gene_name"]].copy()
            ranked["score"] = score.loc[ranked.index].astype(float)
            ranked["gene"] = ranked["gene_name"].astype(str).str.strip()
            ranked = ranked[ranked["gene"].ne("") & ranked["gene"].str.lower().ne("nan")]
            ranked = ranked.sort_values(["gene", "score"], key=lambda s: s.abs() if s.name == "score" else s,
                                        ascending=[True, False], kind="mergesort").drop_duplicates("gene")
            ranked = ranked.sort_values(["score", "gene"], ascending=[False, True], kind="mergesort")
            if ranked.empty:
                raise RuntimeError(f"No expression-filtered gene symbols for secondary ranking: {target}")
            scores = ranked["score"].to_numpy(dtype=float)
            unique_scores = int(np.unique(scores).size)
            primary_qc = qc_primary.loc[target]
            qc_rows.append((target, int(cell_counts[target]), int(len(ranked)), unique_scores,
                            1.0 - unique_scores / len(ranked), float(np.mean(scores == 0)),
                            primary_qc.get("tie_fraction", "NA"), primary_qc.get("zero_effect_fraction", "NA"),
                            primary_qc.get("expressed_genes_before_symbol_collapse", "NA")))
            result = gseapy.prerank(
                rnk=ranked[["gene", "score"]], gene_sets=pathways, outdir=None,
                min_size=min_size, max_size=max_size, permutation_num=permutations,
                seed=int(config["random_seed"]) + i, threads=int(config["analysis_threads"]),
                no_plot=True, verbose=False,
            )
            res = result.res2d.reset_index(drop=True)
            local = []
            universe = set(ranked["gene"].astype(str))
            for item in res.to_dict("records"):
                term = str(item.get("Term", ""))
                pid = term if term in pathways else next((k for k, v in pathway_names.items() if v == term), "")
                pname = pathway_names.get(pid, term)
                lead = str(item.get("Lead_genes", "") or "")
                lead_genes = [x for x in __import__("re").split(r"[;/, ]+", lead) if x]
                try:
                    pval = float(item.get("NOM p-val", 1.0))
                except (TypeError, ValueError):
                    pval = 1.0
                if not math.isfinite(pval):
                    pval = 1.0
                pval = max(pval, 1.0 / (permutations + 1))
                set_genes = pathways.get(pid, set())
                overlap = len(universe & set_genes)
                coverage = overlap / len(set_genes) if set_genes else 0.0
                try:
                    es, nes = float(item.get("ES", math.nan)), float(item.get("NES", math.nan))
                except (TypeError, ValueError):
                    es, nes = math.nan, math.nan
                try:
                    fwer = float(item.get("FWER p-val", math.nan))
                except (TypeError, ValueError):
                    fwer = math.nan
                local.append({"target_gene": target, "target_cells": cell_counts[target], "control_cells": int(config.get("ntc_cells", config["min_control_cells"])),
                              "pathway_id": pid, "pathway_name": pname, "ES": es, "NES": nes, "p_value": pval,
                              "fwer_p_value": fwer, "leading_edge": ";".join(lead_genes), "leading_edge_size": len(lead_genes),
                              "ranked_unique_symbols": len(ranked), "rank_unique_scores": unique_scores,
                              "pathway_gene_count": len(set_genes), "genes_in_ranked_universe": overlap,
                              "coverage_fraction": coverage, "permutations": permutations, "priority_class": "secondary_mean_difference",
                              "main_result_eligible": len(lead_genes) >= 3 and coverage >= float(config["min_pathway_coverage"])})
            if local:
                qvals = bh_adjust(np.asarray([r["p_value"] for r in local], dtype=float))
                for row, q in zip(local, qvals):
                    row["fdr"] = float(q)
                    writer.writerow([row[k] for k in ["target_gene", "target_cells", "control_cells", "pathway_id", "pathway_name", "ES", "NES", "p_value", "fdr",
                                                             "fwer_p_value", "leading_edge", "leading_edge_size", "ranked_unique_symbols", "rank_unique_scores",
                                                             "pathway_gene_count", "genes_in_ranked_universe", "coverage_fraction", "permutations", "priority_class", "main_result_eligible"]])
                    count += 1
            print(f"SECONDARY_GSEA\t{i}\t{target}\tgenes={len(ranked)}\tpathways={len(local)}", flush=True)
    write_frame(qc_path, pd.DataFrame(qc_rows, columns=["target_gene", "target_cells", "secondary_ranked_unique_symbols", "secondary_unique_scores",
                                                        "secondary_tie_fraction", "secondary_zero_effect_fraction", "primary_tie_fraction",
                                                        "primary_zero_effect_fraction", "primary_expressed_genes"]))
    return count


def calc_target_summary(tables: Path, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    cov = read_tsv(tables / "target_guide_coverage.tsv")
    conf = read_tsv(tables / "target_confidence.tsv")
    de_sum = read_tsv(tables / "differential_expression/target_DE_summary.tsv")
    de_path = tables / "differential_expression/target_level_DE_all.tsv.gz"
    effect_records: dict[str, dict] = {}
    for target, frame in iter_target_groups(de_path):
        effect = pd.to_numeric(frame["rank_biserial"], errors="coerce")
        delta = pd.to_numeric(frame["delta_pct"], errors="coerce")
        support = pd.to_numeric(frame["pct_target"], errors="coerce").fillna(0) + pd.to_numeric(frame["pct_control"], errors="coerce").fillna(0)
        # Keep an effect-size maximum and an expression-supported representative gene separately.
        max_idx = effect.abs().idxmax() if effect.notna().any() else None
        if max_idx is None:
            continue
        rows = frame.loc[frame.index]
        maximum = frame.loc[max_idx]
        valid = frame[(support >= 0.10) & effect.notna()].copy()
        if valid.empty:
            valid = frame[effect.notna()].copy()
        v_effect = pd.to_numeric(valid["rank_biserial"], errors="coerce").abs()
        v_support = (pd.to_numeric(valid["pct_target"], errors="coerce").fillna(0)
                     + pd.to_numeric(valid["pct_control"], errors="coerce").fillna(0))
        valid = valid.assign(_abs_effect=v_effect, _support=v_support,
                             _abs_delta=pd.to_numeric(valid["delta_pct"], errors="coerce").abs(),
                             _fdr=pd.to_numeric(valid["fdr"], errors="coerce"))
        rep = valid.sort_values(["_abs_effect", "_support", "_abs_delta", "_fdr"],
                                ascending=[False, False, False, True], kind="mergesort").iloc[0]
        effect_records[target] = {
            "strongest_effect_abs_rank_biserial": abs(float(maximum["rank_biserial"])),
            "strongest_effect_gene": str(maximum["gene_name"]),
            "strongest_effect_delta_pct": as_float(maximum["delta_pct"]),
            "representative_gene_id": str(rep["gene_id"]),
            "representative_gene": str(rep["gene_name"]),
            "representative_rank_biserial": as_float(rep["rank_biserial"]),
            "representative_pct_target": as_float(rep["pct_target"]),
            "representative_pct_control": as_float(rep["pct_control"]),
            "representative_fdr": as_float(rep["fdr"]),
        }
    conf_index = conf.set_index("target_gene")
    de_index = de_sum.set_index("target_gene")
    rows = []
    for r in cov.to_dict("records"):
        target = r["target_gene"]
        n_cells = int(float(r["singlet_cells"]))
        if n_cells < int(config["min_target_cells"]):
            continue
        c = conf_index.loc[target].to_dict()
        d = de_index.loc[target].to_dict() if target in de_index.index else {}
        e = effect_records.get(target, {})
        rows.append({"target_gene": target, "target_singlet_cells": n_cells,
                     "active_guides": int(float(r["active_guides_ge5"])),
                     "largest_guide_share": as_float(r["max_active_guide_fraction"]),
                     "de_genes_fdr_lt_0_05": int(float(d.get("genes_fdr_lt_0_05", 0) or 0)),
                     "genes_p_lt_0_05": int(float(d.get("genes_p_lt_0_05", 0) or 0)),
                     **e,
                     "median_pairwise_spearman": as_float(c.get("median_pairwise_spearman")),
                     "median_directional_concordance": as_float(c.get("median_directional_concordance")),
                     "minimum_loogo_spearman": as_float(c.get("minimum_loogo_spearman")),
                     "minimum_loogo_direction_retention": as_float(c.get("minimum_loogo_direction_retention")),
                     "loogo_stable": str(c.get("loogo_stable", "False")).lower() in {"true", "1"},
                     "evidence_class": c["confidence_class"]})
    summary = pd.DataFrame(rows).sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True], kind="mergesort")
    return summary, pd.DataFrame(effect_records).T.rename_axis("target_gene").reset_index()


def choose_representatives(summary: pd.DataFrame) -> pd.DataFrame:
    records = summary.set_index("target_gene").to_dict("index")
    chosen: list[tuple[str, str]] = []
    def add(target: str | None, reason: str):
        if target and target not in {x[0] for x in chosen}:
            chosen.append((target, reason))
    partial = summary.loc[summary.evidence_class.eq("PARTIAL_SUPPORT"), "target_gene"].tolist()
    if partial:
        add(partial[0], "唯一达到部分 guide 一致性门槛的 target")
    if not summary.empty:
        add(str(summary.sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True]).iloc[0]["target_gene"]),
            "singlet coverage 最高")
    multi = summary[pd.to_numeric(summary["active_guides"], errors="coerce") >= 2].copy()
    balanced = multi[pd.to_numeric(multi["largest_guide_share"], errors="coerce") <= 0.5].copy()
    if not balanced.empty:
        balanced = balanced.sort_values(["median_pairwise_spearman", "median_directional_concordance", "loogo_stable", "target_singlet_cells", "target_gene"],
                                        ascending=[False, False, False, False, True], kind="mergesort", na_position="last")
        add(str(balanced.iloc[0]["target_gene"]), "在 guide-balanced eligible targets 中 guide consistency 指标相对靠前")
    discordant = summary[summary.evidence_class.eq("DISCORDANT")].sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True])
    if not discordant.empty:
        add(str(discordant.iloc[0]["target_gene"]), "展示高覆盖 DISCORDANT 结果及其限制")
    if len(chosen) < 4 and not multi.empty:
        effect_order = multi.sort_values(["strongest_effect_abs_rank_biserial", "target_singlet_cells", "target_gene"],
                                         ascending=[False, False, True], kind="mergesort", na_position="last")
        for target in effect_order["target_gene"]:
            if target not in {x[0] for x in chosen}:
                add(str(target), "展示表达支持下较大的效应量；非 FDR-only 选择")
                break
    output = []
    for target, reason in chosen[:4]:
        row = records[target].copy()
        row["target_gene"] = target
        row["selection_reason"] = reason
        output.append(row)
    return pd.DataFrame(output)


def prepare_cell_state(tables: Path, source_root: Path, h5ad_path: Path, sample_id: str) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray, pd.DataFrame]:
    assignments = read_tsv(tables / "guide_assignment_all.tsv.gz")
    dragen_dir = source_root / "02_dragen_scrna_crispr/04_tool_output/dragen"
    barcodes_mex = read_gzip_lines(dragen_dir / f"{sample_id}.scRNA.filtered.barcodes.tsv.gz")
    with gzip.open(dragen_dir / f"{sample_id}.scRNA.filtered.features.tsv.gz", "rt", encoding="utf-8") as handle:
        gex_features = [line.split("\t", 1)[0] for line in handle if line.rstrip("\r\n").endswith("Gene Expression")]
    adata = ad.read_h5ad(h5ad_path, backed="r")
    h5_barcodes = list(map(str, adata.obs_names))
    h5_features = list(map(str, adata.var_names))
    order = {
        "h5ad_barcode_order_matches_filtered_mex": h5_barcodes == barcodes_mex,
        "h5ad_barcode_order_matches_assignment": h5_barcodes == assignments["barcode"].astype(str).tolist(),
        "h5ad_gene_feature_order_matches_filtered_mex": h5_features == gex_features,
        "h5ad_cells": len(h5_barcodes), "h5ad_gene_features": len(h5_features),
    }
    if not all([order["h5ad_barcode_order_matches_filtered_mex"], order["h5ad_barcode_order_matches_assignment"], order["h5ad_gene_feature_order_matches_filtered_mex"]]):
        adata.file.close()
        raise RuntimeError(f"H5AD/MEX barcode or feature order mismatch: {order}")
    obs = adata.obs.copy()
    obs["barcode"] = h5_barcodes
    for col in ["guide_status", "target_genes", "leiden"]:
        obs[col] = obs[col].astype(str)
    if not (obs["guide_status"].to_numpy() == assignments["guide_status"].to_numpy()).all():
        adata.file.close()
        raise RuntimeError("H5AD guide_status order/values differ from assignment table")
    if not (obs["target_genes"].to_numpy() == assignments["target_genes"].to_numpy()).all():
        adata.file.close()
        raise RuntimeError("H5AD target assignment order/values differ from assignment table")
    singlet = obs["guide_status"].eq("singlet")
    ntc = singlet & obs["target_genes"].eq(CONTROL)
    target_names = sorted(set(obs.loc[singlet & ~obs["target_genes"].isin(["", CONTROL]), "target_genes"]))
    cluster_order = sorted(obs["leiden"].unique(), key=lambda x: (not str(x).isdigit(), int(x) if str(x).isdigit() else str(x)))
    ntc_counts = obs.loc[ntc, "leiden"].value_counts().to_dict()
    state_rows = []
    js_rows = []
    group_sets = [(CONTROL, obs.loc[ntc])]
    for target in target_names:
        subset = obs.loc[singlet & obs["target_genes"].eq(target)]
        group_sets.append((target, subset))
        target_dist = np.asarray([int((subset["leiden"] == c).sum()) for c in cluster_order], dtype=float)
        ntc_dist = np.asarray([int((obs.loc[ntc, "leiden"] == c).sum()) for c in cluster_order], dtype=float)
        distance = float(jensenshannon(target_dist, ntc_dist, base=2)) if target_dist.sum() and ntc_dist.sum() else math.nan
        js_rows.append((target, len(subset), distance))
    for group, frame in group_sets:
        counts = frame["leiden"].value_counts().to_dict()
        total = len(frame)
        for cluster in cluster_order:
            n = int(counts.get(cluster, 0))
            state_rows.append((group, cluster, n, pct(n, total)))
    state_comp = pd.DataFrame(state_rows, columns=["group", "leiden_cluster", "singlet_cells", "fraction_within_group"])
    state_js = pd.DataFrame(js_rows, columns=["target_gene", "target_singlet_cells", "jensen_shannon_distance_to_ntc"])
    overlap_rows = []
    for target in target_names:
        subset = obs.loc[singlet & obs["target_genes"].eq(target)]
        target_counts = subset["leiden"].value_counts().to_dict()
        supported_clusters = [c for c in cluster_order if int(ntc_counts.get(c, 0)) >= 5 and int(target_counts.get(c, 0)) > 0]
        shared_cells = sum(int(target_counts.get(c, 0)) for c in supported_clusters)
        testable = shared_cells >= 20 and len(supported_clusters) >= 2
        overlap_rows.append((target, len(subset), shared_cells, len(supported_clusters), ";".join(map(str, supported_clusters)),
                             sum(int(ntc_counts.get(c, 0)) for c in supported_clusters), "TESTABLE" if testable else "NOT_TESTABLE",
                             "" if testable else "<20 target singlets or fewer than 2 shared Leiden clusters with >=5 NTC singlets"))
    overlap = pd.DataFrame(overlap_rows, columns=["target_gene", "target_singlet_cells", "target_cells_in_shared_states", "shared_leiden_clusters",
                                                 "shared_cluster_ids", "ntc_cells_in_shared_states", "status", "reason"])
    umap = np.asarray(adata.obsm["X_umap"])
    top_targets = set(sorted([(t, int((singlet & obs["target_genes"].eq(t)).sum())) for t in target_names], key=lambda x: (-x[1], x[0]))[:5][i][0] for i in range(min(5, len(target_names))))
    display = np.full(len(obs), "Other passing cells", dtype=object)
    display[ntc.to_numpy()] = "NTC singlets"
    for target in top_targets:
        display[(singlet & obs["target_genes"].eq(target)).to_numpy()] = target
    umap_df = pd.DataFrame({"UMAP1": umap[:, 0], "UMAP2": umap[:, 1], "display_group": display})
    # Summarize target-vs-NTC state composition before closing the backed object.
    state_summary = {"leiden_clusters": len(cluster_order), "ntc_singlet_cells": int(ntc.sum()),
                     "ntc_cells_by_cluster": {str(k): int(v) for k, v in ntc_counts.items()},
                     "ntc_clusters_with_at_least_5_cells": sum(int(v) >= 5 for v in ntc_counts.values()),
                     "targets_with_testable_state_overlap": int((overlap["status"] == "TESTABLE").sum()),
                     "eligible_target_count": int((overlap["target_singlet_cells"].astype(int) >= 20).sum())}
    adata.file.close()
    return order | state_summary, state_comp, state_js, overlap, umap, umap_df


def compute_cell_state_effects(tables: Path, h5ad_path: Path, overlap: pd.DataFrame, top100_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    assignment = read_tsv(tables / "guide_assignment_all.tsv.gz")
    obs = assignment[["barcode", "guide_status", "target_genes"]].copy()
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs["leiden"] = adata.obs["leiden"].astype(str).to_numpy()
    cell_index = {str(barcode): i for i, barcode in enumerate(adata.obs_names)}
    control_mask = obs["guide_status"].eq("singlet") & obs["target_genes"].eq(CONTROL)
    controls_by_cluster: dict[str, np.ndarray] = {}
    for cluster, frame in obs.loc[control_mask].groupby("leiden", sort=False):
        controls_by_cluster[str(cluster)] = np.asarray([cell_index[b] for b in frame["barcode"]], dtype=np.int32)
    var_index = {str(gene): i for i, gene in enumerate(adata.var_names)}
    top = read_tsv(top100_path)
    top["rank_biserial"] = pd.to_numeric(top["rank_biserial"], errors="coerce")
    gene_outputs, summary_rows = [], []
    overlap_map = overlap.set_index("target_gene").to_dict("index")
    layer = adata.layers["log_normalized"]
    for target, genes_frame in top.groupby("target_gene", sort=False):
        state = overlap_map.get(target, {})
        if state.get("status") != "TESTABLE":
            summary_rows.append((target, state.get("status", "NOT_TESTABLE"), state.get("target_singlet_cells", 0),
                                 state.get("target_cells_in_shared_states", 0), state.get("shared_leiden_clusters", 0), 0,
                                 math.nan, math.nan, math.nan))
            continue
        gene_frame = genes_frame.drop_duplicates("gene_id").copy()
        gene_frame["var_idx"] = gene_frame["gene_id"].map(var_index)
        gene_frame = gene_frame[gene_frame["var_idx"].notna()].copy()
        if gene_frame.empty:
            summary_rows.append((target, "NOT_TESTABLE", state["target_singlet_cells"], state["target_cells_in_shared_states"], state["shared_leiden_clusters"], 0,
                                 math.nan, math.nan, math.nan))
            continue
        gene_frame["var_idx"] = gene_frame["var_idx"].astype(int)
        gene_indices = gene_frame["var_idx"].to_numpy(dtype=np.int32)
        primary_effect = gene_frame["rank_biserial"].to_numpy(dtype=float)
        target_mask = obs["guide_status"].eq("singlet") & obs["target_genes"].eq(target)
        target_by_cluster = {}
        for cluster, frame in obs.loc[target_mask].groupby("leiden", sort=False):
            target_by_cluster[str(cluster)] = np.asarray([cell_index[b] for b in frame["barcode"]], dtype=np.int32)
        shared = [c for c in sorted(controls_by_cluster) if len(controls_by_cluster[c]) >= 5 and c in target_by_cluster]
        shared_n = sum(len(target_by_cluster[c]) for c in shared)
        if shared_n < 20 or len(shared) < 2:
            summary_rows.append((target, "NOT_TESTABLE", int(target_mask.sum()), shared_n, len(shared), 0, math.nan, math.nan, math.nan))
            continue
        adjusted = np.zeros(len(gene_indices), dtype=float)
        total_weight = 0.0
        for cluster in shared:
            t_idx = target_by_cluster[cluster]
            c_idx = controls_by_cluster[cluster]
            xt = layer[t_idx, :][:, gene_indices].toarray()
            xc = layer[c_idx, :][:, gene_indices].toarray()
            u = mannwhitneyu(xt, xc, axis=0, alternative="two-sided", method="asymptotic").statistic
            effect = 2.0 * np.asarray(u, dtype=float) / (len(t_idx) * len(c_idx)) - 1.0
            weight = len(t_idx) / shared_n
            adjusted += weight * effect
            total_weight += weight
        if total_weight:
            adjusted /= total_weight
        direction_match = np.sign(primary_effect) == np.sign(adjusted)
        rho = spearmanr(primary_effect, adjusted).statistic if len(primary_effect) > 2 else math.nan
        median_delta = float(np.median(np.abs(adjusted - primary_effect)))
        concordance = float(np.mean(direction_match))
        summary_rows.append((target, "TESTED_EXPLORATORY", int(target_mask.sum()), shared_n, len(shared), len(gene_indices),
                             concordance, rho, median_delta))
        for (_, row), adj in zip(gene_frame.iterrows(), adjusted):
            gene_outputs.append((target, row["gene_id"], row["gene_name"], int(target_mask.sum()), shared_n, len(shared),
                                 float(row["rank_biserial"]), float(adj), bool(np.sign(row["rank_biserial"]) == np.sign(adj))))
        print(f"CELL_STATE_SENSITIVITY\t{target}\tshared_cells={shared_n}\tclusters={len(shared)}\tgenes={len(gene_indices)}", flush=True)
    adata.file.close()
    summary = pd.DataFrame(summary_rows, columns=["target_gene", "status", "target_singlet_cells", "target_cells_in_shared_states", "shared_leiden_clusters",
                                                  "genes_compared", "direction_concordance", "spearman_primary_vs_state_adjusted", "median_abs_effect_change"])
    detailed = pd.DataFrame(gene_outputs, columns=["target_gene", "gene_id", "gene_name", "target_singlet_cells", "target_cells_in_shared_states", "shared_leiden_clusters",
                                                  "primary_rank_biserial", "state_adjusted_rank_biserial", "direction_concordant"])
    return summary, detailed


def save_figures(out: Path, summary: dict, funnel: pd.DataFrame, coverage_bins: pd.DataFrame, confidence: pd.DataFrame,
                 target_summary: pd.DataFrame, representatives: pd.DataFrame, primary_qc: pd.DataFrame,
                 secondary_qc: pd.DataFrame, pathway_sensitivity: pd.DataFrame, state_comp: pd.DataFrame,
                 umap_df: pd.DataFrame, cell_state_summary: pd.DataFrame) -> list[dict]:
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "regular", "axes.labelcolor": "#333333",
                         "text.color": "#222222", "axes.edgecolor": "#555555", "figure.facecolor": "white", "axes.facecolor": "white"})
    manifest = []
    def save(fig, filename, title, section, caption):
        path = fig_dir / filename
        fig.savefig(path, dpi=450, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        manifest.append({"file": filename, "title": title, "section": section, "caption": caption, "source_analysis_version": "v3"})

    # Cell partition + target funnel; cell fates are not represented as a serial loss chain.
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0), gridspec_kw={"width_ratios": [1.05, 1.55]}, constrained_layout=True)
    ax = axes[0]
    cell_metrics = [("Singlet", int(summary["singlet_cells"]), NPG[2]), ("Multi-guide", int(summary["multi_guide_cells"]), NPG[1]),
                    ("Unassigned", int(summary["unassigned_cells"]), "#B8BDC6")]
    total = int(summary["filtered_cells"])
    fractions = [100 * value / total for _, value, _ in cell_metrics]
    ax.barh(np.arange(len(cell_metrics)), fractions, color=[c for _, _, c in cell_metrics], height=.62)
    ax.set_yticks(np.arange(len(cell_metrics)), [name for name, _, _ in cell_metrics], fontsize=9)
    ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    ax.set_title(f"Assignment composition (n={total:,})", fontsize=12)
    ax.set_xlabel("Share of passing cells")
    for y, ((_, value, _), fraction) in enumerate(zip(cell_metrics, fractions)):
        ax.text(fraction + 1.2, y, f"{value:,} ({fraction:.2f}%)", ha="left", va="center", fontsize=8.5, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax = axes[1]
    vals = funnel["count"].astype(int).to_numpy()
    labels = funnel["stage"].tolist()
    colors = ["#8491B4", "#4DBBD5", "#00A087", "#3C5488", "#E64B35"]
    y = np.arange(len(vals))
    ax.barh(y, vals, color=colors, height=.62)
    ax.invert_yaxis(); ax.set_yticks(y, labels, fontsize=8.5)
    ax.set_xlim(0, max(vals) * 1.12)
    for yi, value in zip(y, vals):
        ax.text(value + max(vals)*.015, yi, f"{value:,}", va="center", fontsize=9, fontweight="bold")
    ax.set_xlabel("Observed targets")
    ax.set_title("Target analysis funnel", fontsize=12)
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "fig01_screen_usability_funnel.png", "Screen usability funnel", "overview",
         "Passing-cell guide classes are mutually exclusive. Target stages use ≥20 singlet cells, ≥2 active guides, largest-guide share ≤0.5, then median rho ≥0 and direction concordance ≥0.5; final evidence classes remain the all-target partition.")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2), constrained_layout=True)
    ax = axes[0]
    read_values = [summary["fraction_reads_in_passing_cells"]*100, summary["feature_fraction_reads_in_passing_cells"]*100]
    ax.bar(["GEX reads", "Feature reads"], read_values, color=[NPG[1], NPG[2]], width=.58)
    ax.set_ylim(0, 100); ax.set_ylabel("Reads in passing cells (%)")
    for i, value in enumerate(read_values): ax.text(i, value + 2, f"{value:.2f}%", ha="center", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax = axes[1]
    vals = [summary["sequencing_saturation"]*100, summary["mapq40_rate_gex"]*100, summary["q30_rate_gex"]*100]
    ax.bar(["Saturation", "MAPQ ≥40", "Q30 bases"], vals, color=[NPG[5], NPG[3], NPG[0]], width=.58)
    ax.set_ylim(0, 100); ax.set_ylabel("Percent")
    for i, value in enumerate(vals): ax.text(i, value + 2, f"{value:.2f}%", ha="center", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax = axes[2]
    ax.bar(["Genes/cell", "Molecules/cell"], [summary["median_genes_per_cell"], summary["median_molecules_per_cell"]], color=[NPG[4], NPG[1]], width=.58)
    for i, value in enumerate([summary["median_genes_per_cell"], summary["median_molecules_per_cell"]]): ax.text(i, value+20, f"{value:,}", ha="center", fontsize=9)
    ax.set_ylabel("Median per passing cell"); ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig02_gex_crispr_qc.png", "GEX and CRISPR capture QC", "qc",
         "Passing-read fractions, sequencing saturation, mapping/Q30 and per-cell medians use the DRAGEN denominators shown in the QC tables. Saturation is descriptive; no pass/fail threshold is applied.")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), constrained_layout=True)
    ax = axes[0]
    bins = coverage_bins
    bar_colors = ["#B8BDC6", NPG[5], NPG[1], NPG[2], NPG[0]]
    bars = ax.bar(bins["bin"], bins["target_count"].astype(int), color=bar_colors, width=.68)
    ax.bar_label(bars, fmt="%d", padding=3, fontsize=8)
    ax.set_ylabel("Observed targets"); ax.set_xlabel("Singlet cells per target")
    ax.spines[["top", "right"]].set_visible(False)
    ax = axes[1]
    eligible = target_summary.sort_values(["target_singlet_cells", "target_gene"], ascending=[False, True]).head(20).sort_values("target_singlet_cells")
    colors = [CLASS_COLORS.get(x, "#888888") for x in eligible["evidence_class"]]
    ax.barh(eligible["target_gene"], eligible["target_singlet_cells"], color=colors)
    ax.set_xlabel("Singlet guide cells"); ax.set_ylabel("")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "fig03_target_coverage.png", "Target coverage distribution and top coverage", "coverage",
         "The left panel shows target counts in fixed singlet-cell bins. The right panel shows the highest-coverage eligible targets, colored by the shared evidence-class legend.")

    fig, ax = plt.subplots(figsize=(8.4, 6.0), constrained_layout=True)
    plotted = confidence.copy()
    plotted["x"] = pd.to_numeric(plotted["median_pairwise_spearman"], errors="coerce")
    plotted["y"] = pd.to_numeric(plotted["minimum_loogo_spearman"], errors="coerce")
    plotted["size"] = pd.to_numeric(plotted["target_cells"], errors="coerce").fillna(1)
    plotted["guides"] = pd.to_numeric(plotted["active_guides_ge5"], errors="coerce").fillna(0)
    plotted = plotted[(plotted["target_cells"].astype(float) >= 20) & plotted["x"].notna() & plotted["y"].notna()]
    plotted["guide_bucket"] = np.where(plotted["guides"] >= 3, "≥3 active guides", "2 active guides")
    guide_markers = {"2 active guides": "o", "≥3 active guides": "^"}
    for cls in CLASS_ORDER:
        group = plotted[plotted["confidence_class"] == cls]
        if group.empty: continue
        for guide_bucket, marker in guide_markers.items():
            subset = group[group["guide_bucket"] == guide_bucket]
            if subset.empty: continue
            ax.scatter(subset["x"], subset["y"], s=22 + 2.1*np.sqrt(subset["size"]), c=CLASS_COLORS[cls], alpha=.82,
                       marker=marker, edgecolors="white", linewidths=.4)
    ax.axvline(0, color="#666666", linewidth=.8)
    ax.axhline(.8, color="#777777", linestyle="--", linewidth=.8)
    ax.set_xlabel("Median guide-guide Spearman rho")
    ax.set_ylabel("Minimum leave-one-guide-out Spearman rho")
    ax.set_title("Multi-guide consistency and LOOGO stability")
    class_handles = [ax.scatter([], [], s=42, c=CLASS_COLORS[cls], marker="o", edgecolors="white", linewidths=.4, label=cls)
                     for cls in CLASS_ORDER if cls in set(plotted["confidence_class"])]
    class_legend = ax.legend(handles=class_handles, title="Evidence class", frameon=False, fontsize=7, title_fontsize=7,
                             ncol=2, loc="lower left")
    ax.add_artist(class_legend)
    guide_handles = [ax.scatter([], [], s=48, c="#566574", marker=marker, edgecolors="white", linewidths=.4, label=name)
                     for name, marker in guide_markers.items() if name in set(plotted["guide_bucket"])]
    ax.legend(handles=guide_handles, title="Active guides", frameon=False, fontsize=7, title_fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    label_targets = []
    for cls in ["PARTIAL_SUPPORT", "GUIDE_DOMINATED", "DISCORDANT"]:
        candidate = plotted[plotted["confidence_class"] == cls].sort_values(["size", "target_gene"], ascending=[False, True])
        if not candidate.empty:
            label_targets.append(candidate.iloc[0])
    for row in label_targets:
        ax.annotate(row["target_gene"], (row["x"], row["y"]), xytext=(5, 4), textcoords="offset points", fontsize=7.5)
    save(fig, "fig04_multi_guide_robustness.png", "Multi-guide robustness map", "robustness",
         "Eligible targets: x=median guide-guide rho, y=minimum LOOGO rho, point size=target singlet cells, color=evidence class, shape=active guides (2 or ≥3). Selected labels identify the highest-coverage target per PARTIAL_SUPPORT, GUIDE_DOMINATED and DISCORDANT class; direction retention is reported separately and no composite score is used.")

    rep = representatives.copy()
    fig, ax = plt.subplots(figsize=(8.8, max(3.8, .75*max(1, len(rep)))), constrained_layout=True)
    if not rep.empty:
        y = np.arange(len(rep))
        effect = pd.to_numeric(rep["representative_rank_biserial"], errors="coerce").fillna(0).to_numpy()
        support = (pd.to_numeric(rep["representative_pct_target"], errors="coerce").fillna(0)
                   + pd.to_numeric(rep["representative_pct_control"], errors="coerce").fillna(0)).to_numpy()
        sizes = 36 + 110*np.clip(support/2, 0, 1)
        colors = [CLASS_COLORS.get(x, "#777777") for x in rep["evidence_class"]]
        ax.scatter(effect, y, s=sizes, c=colors, alpha=.9, edgecolors="white", linewidths=.6)
        labels = [f"{r.target_gene} · {r.representative_gene}" for r in rep.itertuples()]
        ax.set_yticks(y, labels, fontsize=8)
    ax.axvline(0, color="#555555", linewidth=.8)
    ax.set_xlabel("Signed rank-biserial effect (target vs NTC)")
    ax.set_title("Illustrative expression effects; selection reasons shown in table")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "fig05_representative_perturbation_effects.png", "Representative perturbation expression effects", "effects",
         "Examples are selected using coverage, expression support and multi-guide evidence; they are not ranked as validated perturbations. Point size reflects expression prevalence support.")

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), constrained_layout=True)
    for ax, frame, col, title, color in [
        (axes[0], primary_qc, "tie_fraction", "Primary ranking: signed rank-biserial ties", NPG[3]),
        (axes[1], secondary_qc, "secondary_tie_fraction", "Sensitivity ranking: normalized mean difference ties", NPG[1]),
    ]:
        vals = pd.to_numeric(frame[col], errors="coerce").dropna().to_numpy()
        ax.hist(vals*100, bins=14, color=color, edgecolor="white")
        med = float(np.median(vals))*100 if len(vals) else math.nan
        ax.axvline(med, color="#222222", linestyle="--", linewidth=1, label=f"Median {med:.1f}%")
        ax.set_xlabel("Tied rank scores (%)"); ax.set_ylabel("Targets"); ax.set_title(title, fontsize=10); ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig06_pathway_ranking_qc.png", "Reactome ranking quality control", "pathways",
         "Tie-fraction distributions are shown for the primary and secondary rankings. Zero-effect fraction and retained gene counts are given in the accompanying ranking-QC table.")

    ps = pathway_sensitivity.copy()
    ps["primary_NES"] = pd.to_numeric(ps.get("primary_NES"), errors="coerce")
    ps["secondary_NES"] = pd.to_numeric(ps.get("secondary_NES"), errors="coerce")
    fig, ax = plt.subplots(figsize=(6.8, 5.6), constrained_layout=True)
    if not ps.empty:
        c = [NPG[2] if bool(x) else NPG[0] for x in ps["direction_concordant"]]
        ax.scatter(ps["primary_NES"], ps["secondary_NES"], c=c, s=24, alpha=.72, linewidths=0)
        ax.scatter([], [], c=NPG[2], s=34, label="Direction concordant")
        ax.scatter([], [], c=NPG[0], s=34, label="Direction discordant")
    ax.axhline(0, color="#777777", linewidth=.7); ax.axvline(0, color="#777777", linewidth=.7)
    ax.set_xlabel("Primary NES (signed rank-biserial)"); ax.set_ylabel("Sensitivity NES (mean-expression difference)")
    ax.set_title("Pathway direction sensitivity")
    if not ps.empty:
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig07_reactome_ranking_sensitivity.png", "Reactome direction sensitivity", "pathways",
         "Each point is a target-pathway pair tested under both rankings; the diagonal quadrants show whether enrichment direction is preserved. All pathway signals remain exploratory.")

    display_groups = [CONTROL] + [t for t in target_summary.sort_values("target_singlet_cells", ascending=False)["target_gene"].tolist() if t != CONTROL][:7]
    sub = state_comp[state_comp["group"].isin(display_groups)].copy()
    clusters = list(dict.fromkeys(state_comp["leiden_cluster"].tolist()))
    fig, ax = plt.subplots(figsize=(9.6, 4.8), constrained_layout=True)
    bottoms = np.zeros(len(display_groups))
    for j, cluster in enumerate(clusters):
        values = []
        for group in display_groups:
            row = sub[(sub["group"] == group) & (sub["leiden_cluster"] == cluster)]
            values.append(float(row.iloc[0]["fraction_within_group"]) if not row.empty else 0.0)
        ax.bar(display_groups, np.asarray(values)*100, bottom=bottoms, color=NPG[j % len(NPG)], label=f"Leiden cluster {cluster}", width=.72)
        bottoms += np.asarray(values)*100
    ax.set_ylabel("Singlet cells per group (%)"); ax.tick_params(axis="x", labelrotation=30)
    ax.legend(frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig08_cell_state_composition.png", "Leiden-cluster composition of NTC and high-coverage targets", "exploratory",
         "Uses existing, unannotated Leiden cluster labels in the input H5AD. Bars are descriptive within-group fractions; no cell-level significance or mechanism is inferred.")

    fig, ax = plt.subplots(figsize=(7.8, 6.2), constrained_layout=True)
    colors = {"Other passing cells": "#D9DDE2", "NTC singlets": NPG[0]}
    for i, g in enumerate([x for x in umap_df["display_group"].unique() if x not in colors]): colors[g] = NPG[(i+2) % len(NPG)]
    for group in ["Other passing cells"] + [x for x in umap_df["display_group"].unique() if x != "Other passing cells"]:
        f = umap_df[umap_df["display_group"] == group]
        if f.empty: continue
        ax.scatter(f["UMAP1"], f["UMAP2"], s=1.6 if group == "Other passing cells" else 5.0,
                   c=colors[group], alpha=.38 if group == "Other passing cells" else .86, label=group, linewidths=0, rasterized=True)
    ax.set_aspect("equal", adjustable="datalim"); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title("Existing embedding from input H5AD")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig09_existing_umap_exploratory.png", "Existing H5AD UMAP (not recomputed)", "exploratory",
         "Existing embedding from input H5AD. Not recomputed in this analysis. Proximity in this descriptive projection is not a perturbation mechanism.")
    return manifest


def main() -> int:
    args = parse_args()
    global SAMPLE_ID, H5AD, FEATURE_REFERENCE, FILTERED_FEATURE_MATRIX, FILTERED_FEATURES, FILTERED_BARCODES, REACTOME_GMT, CONTROL
    SAMPLE_ID = args.sample_id
    CONTROL = args.control_gene
    v2, source, out, runtime = args.v2_root.resolve(), args.source_root.resolve(), args.output_root.resolve(), args.runtime_root.resolve()
    dragen_dir = source / "02_dragen_scrna_crispr/04_tool_output/dragen"
    H5AD = source / "03_cell_guide_qc/04_tool_output" / f"{SAMPLE_ID}_gex_filtered.h5ad"
    FEATURE_REFERENCE = source / "01_input_preflight/04_tool_output/dragen_feature_reference.csv"
    FILTERED_FEATURE_MATRIX = dragen_dir / f"{SAMPLE_ID}.scRNA.filtered.features.matrix.tsv.gz"
    FILTERED_FEATURES = dragen_dir / f"{SAMPLE_ID}.scRNA.filtered.features.tsv.gz"
    FILTERED_BARCODES = dragen_dir / f"{SAMPLE_ID}.scRNA.filtered.barcodes.tsv.gz"
    REACTOME_GMT = args.reactome_gmt.resolve()
    runtime.mkdir(parents=True, exist_ok=True)
    (runtime / "mpl").mkdir(parents=True, exist_ok=True)
    if not v2.is_dir() or not source.is_dir():
        raise FileNotFoundError("V2 result root or native DRAGEN source root is missing")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(v2 / "tables", out / "tables", dirs_exist_ok=True)
    (out / "figures/read_flow_sankey").mkdir(parents=True, exist_ok=True)
    shutil.copytree(v2 / "figures/read_flow_sankey", out / "figures/read_flow_sankey", dirs_exist_ok=True)
    config = read_json(v2 / "methods/reanalysis.v2.json")
    summary = read_json(v2 / "tables/reanalysis_summary.json")
    project = read_json(v2 / "tables/project_summary.json")
    config.update({"report_version": 3, "analysis_output_version": "v3-package-with-reused-primary-DE",
                   "source_analysis_version": "v2-primary-methods-preserved", "v2_root": str(v2),
                   "source_root": str(source), "output_root": str(out), "runtime_root": str(runtime),
                   "reactome_gmt": str(REACTOME_GMT), "secondary_ranking": "mean_target_normalized - mean_control_normalized",
                   "min_target_cells_sensitivity": [10, 20, 30], "min_active_guides_per_target_sensitivity": [1, 2, 3],
                   "max_guide_fraction_sensitivity": [0.4, 0.5, 0.6], "spearman_supported_sensitivity": [0.2, 0.3, 0.4],
                   "direction_discordant_sensitivity": [0.4, 0.5, 0.6], "direction_supported_sensitivity": [0.7, 0.8, 0.9],
                   "loogo_spearman_sensitivity": [0.7, 0.8, 0.9], "loogo_direction_retention_sensitivity": [0.7, 0.8, 0.9],
                   "cell_state_method": "Leiden-stratified rank-biserial direction sensitivity on per-target top 100 genes; no p/FDR"})
    baseline = {"min_target_cells": 20, "min_active_guides": 2, "max_guide_fraction": 0.5,
                "direction_discordant_below": 0.5, "spearman_supported_at_least": 0.3,
                "direction_supported_at_least": 0.8, "loogo_spearman_min": 0.8,
                "loogo_direction_retention_min": 0.8}
    gsea_config = summary.get("analysis_config", {})
    config.update(gsea_config)
    config.update({"report_version": 3, "analysis_output_version": "v3-package-with-reused-primary-DE",
                   "source_analysis_version": "v2-primary-methods-preserved", "v2_root": str(v2),
                   "source_root": str(source), "output_root": str(out), "runtime_root": str(runtime),
                   "reactome_gmt": str(REACTOME_GMT), "secondary_ranking": "mean_target_normalized - mean_control_normalized",
                   **{f"sensitivity_{k}": v for k, v in config.items() if k.endswith("_sensitivity")}})
    config.update({"min_active_guides_per_target": 2, "ntc_cells": int(project["ntc_cells"]), **baseline})
    tables = out / "tables"
    # Keep the old conserved read-flow contract and its validated assets as an appendix input.
    runtime_source = Path(read_json(v2 / "methods/reanalysis.v2.json")["runtime_root"])
    flow_fragment = runtime_source / "sankey/output/reports/read_flow_sankey_section.html"
    if flow_fragment.is_file():
        shutil.copy2(flow_fragment, tables / "read_flow/read_flow_sankey_section.html")

    conf = read_tsv(tables / "target_confidence.tsv")
    cov = read_tsv(tables / "target_guide_coverage.tsv")
    guide = read_tsv(tables / "guide_coverage.tsv")
    assignment = read_tsv(tables / "guide_assignment_all.tsv.gz")
    observed_target_count = int(project["observed_targets"])
    if len(conf) != observed_target_count or len(cov) != observed_target_count or set(conf.target_gene) != set(cov.target_gene):
        raise RuntimeError("Target coverage and confidence table do not match the observed target count")
    primary_classes = dict(zip(conf["target_gene"], conf["confidence_class"]))
    class_counts = Counter(primary_classes.values())
    if sum(class_counts.values()) != int(project["observed_targets"]):
        raise RuntimeError("Evidence class counts do not sum to observed targets")
    base_conf_thresholds = {"min_target_cells": int(config["min_target_cells"]), "min_active_guides": 2,
                            "max_guide_fraction": float(config["max_guide_fraction"]),
                            "direction_discordant_below": float(config["direction_discordant_below"]),
                            "spearman_supported_at_least": float(config["spearman_supported_at_least"]),
                            "direction_supported_at_least": float(config["direction_supported_at_least"]),
                            "loogo_spearman_min": float(config["loogo_spearman_min"]),
                            "loogo_direction_retention_min": float(config["loogo_direction_retention_min"])}
    recomputed_primary = {t: classify_row(r, base_conf_thresholds) for t, r in zip(conf.target_gene, conf.to_dict("records"))}
    if recomputed_primary != primary_classes:
        raise RuntimeError("Recomputed baseline evidence classes differ from V2 classes")

    # Cell and target funnels, coverage bins and class table.
    cov_num = cov.copy()
    for c in ["singlet_cells", "active_guides_ge5", "max_active_guide_fraction"]:
        cov_num[c] = pd.to_numeric(cov_num[c], errors="coerce")
    eligible20 = cov_num[cov_num.singlet_cells >= 20]
    two_guides = eligible20[eligible20.active_guides_ge5 >= 2]
    balanced = two_guides[two_guides.max_active_guide_fraction <= 0.5]
    conf_index = conf.set_index("target_gene")
    consistent = []
    for r in balanced.to_dict("records"):
        c = conf_index.loc[r["target_gene"]]
        rho, direction = as_float(c["median_pairwise_spearman"]), as_float(c["median_directional_concordance"])
        if math.isfinite(rho) and rho >= 0 and math.isfinite(direction) and direction >= 0.5:
            consistent.append(r["target_gene"])
    funnel = pd.DataFrame([
        ("Observed targets", len(cov_num), "All targets with at least one assigned singlet"),
        ("≥20 singlet cells", len(eligible20), "Primary target-level analysis minimum"),
        ("≥2 active guides", len(two_guides), "Each active guide has ≥5 singlet cells"),
        ("Guide-balanced", len(balanced), "Largest active-guide share ≤0.5"),
        ("Guide-consistent", len(consistent), "Median rho ≥0 and direction concordance ≥0.5; includes partial support"),
    ], columns=["stage", "count", "definition"])
    write_frame(tables / "screen_target_funnel.tsv", funnel)
    bins = [("<5", 0, 4), ("5–9", 5, 9), ("10–19", 10, 19), ("20–29", 20, 29), ("≥30", 30, math.inf)]
    cov_dist = pd.DataFrame([{"bin": name, "target_count": int(((cov_num.singlet_cells >= lo) & (cov_num.singlet_cells <= hi)).sum()),
                              "fraction_of_observed_targets": float(((cov_num.singlet_cells >= lo) & (cov_num.singlet_cells <= hi)).sum()) / len(cov_num)}
                             for name, lo, hi in bins])
    write_frame(tables / "target_coverage_distribution.tsv", cov_dist)
    class_rows = []
    definitions = {
        "LOW_CELL": "<20 target singlet cells; insufficient target-level coverage",
        "LOW_GUIDE_SUPPORT": "<2 active guides, each active at ≥5 singlet cells",
        "GUIDE_DOMINATED": "largest active-guide share >0.5",
        "DISCORDANT": "median guide rho <0 or direction concordance <0.5",
        "PARTIAL_SUPPORT": "rho ≥0 and direction concordance ≥0.5, without all SUPPORTED gates",
        "SUPPORTED": "rho ≥0.3, direction concordance ≥0.8, and both LOOGO minima ≥0.8 after earlier gates",
        "INCONCLUSIVE": "does not meet another evidence-class rule",
    }
    for cls in CLASS_ORDER:
        n = int(class_counts.get(cls, 0))
        class_rows.append((cls, n, pct(n, len(conf)), definitions[cls]))
    class_table = pd.DataFrame(class_rows, columns=["evidence_class", "target_count", "fraction", "rule"])
    write_frame(tables / "target_evidence_class_summary.tsv", class_table)
    target_full = cov.merge(conf[["target_gene", "median_pairwise_spearman", "median_directional_concordance",
                                 "minimum_loogo_spearman", "minimum_loogo_direction_retention", "loogo_stable", "confidence_class"]],
                            on="target_gene", how="left", validate="one_to_one")
    target_full = target_full.sort_values(["singlet_cells", "target_gene"], ascending=[False, True], kind="mergesort")
    write_frame(tables / "target_coverage_full.tsv", target_full)

    # NTC guide support is measured from both guide assignment and the filtered feature MEX.
    ntc_guides = guide[guide.target_gene.eq(CONTROL)].copy()
    if ntc_guides.empty:
        raise RuntimeError("No NTC guide features were found in guide coverage")
    ref = pd.read_csv(FEATURE_REFERENCE, dtype=str, keep_default_na=False)
    ntc_ids = set(ref.loc[ref["gene_name"].eq(CONTROL), "id"].astype(str))
    unique_ntc_barcodes, positive_per_guide = set(), Counter()
    with gzip.open(FILTERED_FEATURE_MATRIX, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 5 and fields[0] in ntc_ids and int(fields[4]) > 0:
                unique_ntc_barcodes.add(fields[3])
                positive_per_guide[fields[0]] += 1
    if set(ntc_guides.sgRNA_id) != ntc_ids:
        raise RuntimeError("NTC guide reference and guide-coverage rows differ")
    ntc_guides["matrix_positive_cells"] = ntc_guides.sgRNA_id.map(positive_per_guide).fillna(0).astype(int)
    ntc_guides["singlet_cells"] = pd.to_numeric(ntc_guides["singlet_cells"], errors="coerce").fillna(0).astype(int)
    ntc_guides["active_for_effects"] = ntc_guides["active_for_effects"].astype(str).str.lower().isin(["true", "1"])
    if not (ntc_guides["matrix_positive_cells"].to_numpy() == pd.to_numeric(ntc_guides["positive_cells_from_matrix"], errors="coerce").fillna(0).astype(int).to_numpy()).all():
        raise RuntimeError("NTC per-guide feature-positive cell counts differ between MEX and guide coverage")
    ntc_singlet_cells = int(ntc_guides.singlet_cells.sum())
    active_ntc = int(ntc_guides.active_for_effects.sum())
    ntc_assignment = assignment[assignment.target_genes.astype(str).str.split(";").apply(lambda x: CONTROL in x)]
    ntc_status_counts = Counter(ntc_assignment.guide_status)
    if ntc_singlet_cells != int(project["ntc_cells"]) or ntc_status_counts.get("singlet", 0) != ntc_singlet_cells:
        raise RuntimeError("NTC singlet cell counts differ across guide table, assignments, and summary")
    ntc_guides = ntc_guides[["sgRNA_id", "sequence", "matrix_positive_cells", "singlet_cells", "total_feature_umi",
                             "median_feature_umi_per_positive_cell", "active_for_effects"]].sort_values("sgRNA_id")
    write_frame(tables / "ntc_guide_qc.tsv", ntc_guides)
    ntc_summary = {"ntc_guide_features": int(len(ntc_guides)), "active_ntc_guides_ge5_singlets": active_ntc,
                   "ntc_singlet_cells": ntc_singlet_cells, "ntc_feature_positive_cells_unique": int(len(unique_ntc_barcodes)),
                   "ntc_per_guide_positive_cell_sum_not_unique": int(ntc_guides.matrix_positive_cells.sum()),
                   "ntc_cells_with_assignment_singlet": int(ntc_status_counts.get("singlet", 0)),
                   "ntc_cells_with_assignment_multi": int(ntc_status_counts.get("multi", 0)),
                   "control_reference_is_limited": True}

    sensitivity, transitions = threshold_sensitivity(conf, primary_classes, baseline)
    write_frame(tables / "target_threshold_sensitivity.tsv", sensitivity)
    write_frame(tables / "target_threshold_class_transitions.tsv", transitions)

    target_summary, effect_summary = calc_target_summary(tables, config)
    representatives = choose_representatives(target_summary)
    write_frame(tables / "target_level_summary.tsv", target_summary)
    write_frame(tables / "target_effect_gene_summary.tsv", effect_summary)
    write_frame(tables / "representative_target_effects.tsv", representatives)

    # Read-only H5AD validation, state-composition summaries and conditional stratified sensitivity.
    order_and_state, state_comp, state_js, overlap, umap, umap_df = prepare_cell_state(tables, source, H5AD, SAMPLE_ID)
    write_frame(tables / "cell_state_composition.tsv", state_comp)
    write_frame(tables / "cell_state_distance_to_ntc.tsv", state_js)
    write_frame(tables / "cell_state_overlap.tsv", overlap)
    state_effect_summary, state_effect_genes = compute_cell_state_effects(tables, H5AD, overlap, tables / "differential_expression/target_level_DE_top100.tsv")
    write_frame(tables / "cell_state_effect_sensitivity.tsv", state_effect_summary)
    write_frame(tables / "cell_state_effect_sensitivity_genes.tsv.gz", state_effect_genes)

    # Reactome v97 source labels are corrected using the GMT name/ID contract; primary statistics are unchanged.
    pathway_sets, pathway_names, name_to_id = gmt_maps(REACTOME_GMT)
    primary_path = tables / "pathway/reactome_GSEA_all.tsv.gz"
    corrected_rows = update_primary_reactome_labels(primary_path, name_to_id, pathway_names)
    primary_qc = read_tsv(tables / "gsea_ranking_qc.tsv")
    secondary_count = 0
    if not args.skip_pathway_sensitivity:
        secondary_output = tables / "pathway/reactome_GSEA_mean_difference.tsv.gz"
        secondary_qc_output = tables / "gsea_ranking_sensitivity_qc.tsv"
        if args.reuse_pathway_sensitivity and secondary_output.is_file() and secondary_qc_output.is_file():
            with gzip.open(secondary_output, "rt", encoding="utf-8") as handle:
                secondary_count = max(0, sum(1 for _ in handle) - 1)
            secondary_qc_resume = read_tsv(secondary_qc_output)
            eligible_resume = set(cov.loc[pd.to_numeric(cov["singlet_cells"], errors="coerce") >= int(config["min_target_cells"]), "target_gene"])
            if set(secondary_qc_resume["target_gene"]) != eligible_resume:
                raise RuntimeError("Cannot resume: secondary ranking QC target set is incomplete")
            print(f"REUSED_COMPLETE_SECONDARY_GSEA\trows={secondary_count}\ttargets={len(secondary_qc_resume)}", flush=True)
        else:
            secondary_count = run_secondary_gsea(tables / "differential_expression/target_level_DE_all.tsv.gz",
                                                 secondary_output,
                                                 secondary_qc_output, config, conf, cov,
                                                 pathway_sets, pathway_names)
        primary = pd.read_csv(primary_path, sep="\t", dtype={"target_gene": str, "pathway_id": str})
        secondary = pd.read_csv(tables / "pathway/reactome_GSEA_mean_difference.tsv.gz", sep="\t", compression="gzip",
                                dtype={"target_gene": str, "pathway_id": str})
        merge_cols = ["target_gene", "pathway_id"]
        merged = primary.merge(secondary, on=merge_cols, how="inner", suffixes=("_primary", "_secondary"), validate="one_to_one")
        merged["direction_concordant"] = (np.sign(pd.to_numeric(merged["NES_primary"], errors="coerce"))
                                           == np.sign(pd.to_numeric(merged["NES_secondary"], errors="coerce")))
        merged["pathway_name"] = merged["pathway_name_primary"]
        merged["target_cells"] = pd.to_numeric(merged["target_cells_primary"], errors="coerce")
        merged["leading_edge_size"] = pd.to_numeric(merged["leading_edge_size_primary"], errors="coerce")
        merged["coverage_fraction"] = pd.to_numeric(merged["coverage_fraction_primary"], errors="coerce")
        merged["primary_FDR"] = pd.to_numeric(merged["fdr_primary"], errors="coerce")
        merged["secondary_FDR"] = pd.to_numeric(merged["fdr_secondary"], errors="coerce")
        merged["primary_NES"] = pd.to_numeric(merged["NES_primary"], errors="coerce")
        merged["secondary_NES"] = pd.to_numeric(merged["NES_secondary"], errors="coerce")
        merged["main_candidate"] = (merged.direction_concordant & (merged.primary_FDR < 0.05) & (merged.secondary_FDR < 0.05)
                                     & (merged.leading_edge_size >= 3) & (merged.coverage_fraction >= float(config["min_pathway_coverage"]))
                                     & (merged.target_cells >= 30))
        sensitivity_path = tables / "pathway/reactome_ranking_sensitivity.tsv.gz"
        keep_cols = ["target_gene", "pathway_id", "pathway_name", "target_cells", "primary_NES", "primary_FDR", "secondary_NES",
                     "secondary_FDR", "direction_concordant", "leading_edge_size", "coverage_fraction", "main_candidate"]
        write_frame(sensitivity_path, merged[keep_cols])
        candidates = merged[merged.main_candidate].sort_values(["primary_FDR", "coverage_fraction", "leading_edge_size"],
                                                               ascending=[True, False, False], kind="mergesort")
        write_frame(tables / "pathway/reactome_exploratory_candidates.tsv", candidates[keep_cols].head(50))
        direction_concordance = float(merged.direction_concordant.mean()) if len(merged) else math.nan
        pathway_summary = {"primary_method": "signed rank-biserial effect size", "sensitivity_method": "mean_target_normalized - mean_control_normalized",
                           "primary_rows": int(len(primary)), "secondary_rows": int(len(secondary)), "shared_target_pathway_pairs": int(len(merged)),
                           "shared_pairs_direction_concordant": int(merged.direction_concordant.sum()),
                           "direction_concordant_fraction": direction_concordance,
                           "jointly_filtered_exploratory_candidates": int(merged.main_candidate.sum()),
                           "display_filter": "FDR<0.05 in both rankings; same NES sign; leading edge>=3; pathway coverage>=10%; target cells>=30",
                           "pathway_label_corrections": corrected_rows}
    else:
        pathway_summary = {"secondary_ranking_status": "SKIPPED_SMOKE_TEST_ONLY", "pathway_label_corrections": corrected_rows}

    project.update({"report_version": 3, "ntc_qc": ntc_summary,
                    "target_funnel": {"observed": int(len(cov)), "minimum_20_cells": int(len(eligible20)),
                                      "minimum_2_active_guides": int(len(two_guides)), "guide_balanced": int(len(balanced)),
                                      "guide_consistent_partial_or_better": int(len(consistent)),
                                      "supported_targets": int(class_counts.get("SUPPORTED", 0)),
                                      "partial_support_targets": int(class_counts.get("PARTIAL_SUPPORT", 0))},
                    "confidence_counts": {k: int(class_counts.get(k, 0)) for k in CLASS_ORDER},
                    "active_guides_de": int(summary.get("active_guides_de", project.get("active_guides_ge5_singlets", 0))),
                    "gsea_primary_tie_median": float(pd.to_numeric(primary_qc.tie_fraction, errors="coerce").median()),
                    "pathway_sensitivity": pathway_summary,
                    "cell_state_sensitivity": order_and_state,
                    "source_analysis_version": "V2 primary DE and guide effect tables revalidated; V3 summaries, figures, state sensitivity and alternate pathway ranking generated"})
    summary.update({"report_version": 3, "source_analysis_version": project["source_analysis_version"],
                    "confidence_counts": project["confidence_counts"], "ntc_qc": ntc_summary,
                    "target_funnel": project["target_funnel"], "eligible_targets_de": int(len(target_summary)),
                    "supported_targets": [], "pathway_sensitivity": pathway_summary,
                    "cell_state_sensitivity": order_and_state, "analysis_config": config,
                    "source_root": str(source), "output_root": str(out)})
    write_json(tables / "project_summary.json", project)
    write_json(tables / "reanalysis_summary.json", summary)
    write_json(out / "methods/reanalysis.v3.json", config)
    analysis_methods = """# scPerturb CRISPR V3 分析方法与判读边界

## 数据来源与复用范围

本版本读取已完成的 DRAGEN 4.5.4-1 原生 PIPseq CRISPR 结果。原生数据、V2 主要差异表达/guide-level effect 表均保留不变；V3 在独立输出目录重算 target funnel、NTC QC、证据分类敏感性、目标摘要与图表，并新增 Leiden 分层效应敏感性和 normalized mean-expression difference Reactome 排名。V3 纠正 Reactome v97 GMT 的名称/稳定 ID 显示映射；primary NES、p 值和 FDR 数值不因标签修正而改变。

GEX feature 与 barcode 顺序逐项核对 filtered MEX 与输入 H5AD。表达展示和 Leiden-cluster 分层敏感性读取 H5AD 的 log_normalized layer；H5AD 的 X 不用作 raw counts。既有 UMAP 坐标沿用，不在 V3 重新拟合。

## Guide assignment、coverage 与 NTC

主统计对象为 passing cells 中 DRAGEN singlet assignment。Active guide 定义保持为每条 guide ≥5 singlet cells；target-level primary analysis minimum 保持为 target ≥20 singlet cells 且 NTC ≥20 cells。NTC singlet cell 数按当前样本实际 assignment 结果统计。NTC feature-positive cells 按 filtered feature MEX 中至少一个 NTC feature UMI>0 的独立 barcode 计数；per-guide positive counts 不相加当作独立细胞数。

Evidence class 保持原有门槛并按顺序判定：LOW_CELL（target cells<20）；LOW_GUIDE_SUPPORT（active guides<2）；GUIDE_DOMINATED（largest active-guide share>0.5）；DISCORDANT（median rho<0 或方向一致性<0.5）；SUPPORTED（rho≥0.3、方向一致性≥0.8 且最小 LOOGO rho 和方向保留均≥0.8）；PARTIAL_SUPPORT（rho≥0 且方向一致性≥0.5，但未满足 SUPPORTED）；其余为 INCONCLUSIVE。分类不是综合分数，也不代表细胞级统计显著性。

Threshold sensitivity 使用 one-factor-at-a-time 网格：target cells 10/20/30、target active-guide minimum 1/2/3、largest-guide share 0.4/0.5/0.6、Spearman supported cutoff 0.2/0.3/0.4、discordance direction boundary 0.4/0.5/0.6、supported direction cutoff 0.7/0.8/0.9，以及两项 LOOGO cutoff 各 0.7/0.8/0.9。每次仅改变一个门槛，primary classification 不随敏感性结果改写。

## 差异表达与多 guide 一致性

复用 V2 完成的逐 target 双侧 Mann–Whitney U 检验、渐近 tie-aware p 值、target 内 Benjamini–Hochberg 校正，以及 signed rank-biserial effect。Primary effect 方法和全基因表保持不变。Target summary 同时给出 FDR<0.05 gene 数、最大绝对 rank-biserial effect、表达支持、active guide 数、largest-guide share、guide-guide Spearman、方向一致性和 LOOGO 指标。代表性示例按 coverage、expression support、guide consistency 与证据类别选择，FDR 不作唯一排序依据。

## Leiden-cluster composition and stratified sensitivity

描述 target 与 NTC singlet cells 在输入 H5AD 的未注释 Leiden cluster labels 中的组成差异；这些 cluster 不直接代表已注释细胞类型或细胞状态。对至少有 20 个 target cells 落在共有 clusters、且至少 2 个 clusters 各有 ≥5 个 NTC singlets 的 target，基于 target cells 在共有 clusters 内的比例，对各 cluster 内 top 100 primary-effect gene 的 signed rank-biserial effects 加权平均；比较 primary 与 cluster-stratified 效应方向。该结果不报告 p/FDR，仅为单样本探索性敏感性分析。条件不足的 target 标为 NOT_TESTABLE。

## Reactome exploratory ranking

Primary ranking 继续使用 signed rank-biserial effect；secondary ranking 使用 target 与 NTC normalized mean expression difference。两者使用相同的表达过滤（任一组表达比例≥5%）、symbol collapse（保留绝对 ranking effect 较大的 feature）、Reactome v97 gene sets（15–500 genes）、每 target 1000 次置换和 target 内 BH。方向敏感性按 target/pathway 的 NES 符号比较。主文候选要求两种 ranking 均 FDR<0.05、NES 方向一致、leading edge≥3、pathway coverage≥10%、target cells≥30。Ranking ties、zero-effect fraction 与保留基因数随结果公开；全部 Reactome 信号均为探索性结果，低证据或非 SUPPORTED target 的 pathway enrichment 不构成已验证机制。

## 解释边界

本数据只有一个生物学样本。Cell-level p/FDR 不替代 biological replicate inference。Perturbation mode 为 UNKNOWN；target RNA 变化仅是表达观察，不能确定分子扰动模式或编辑结果。Guide imbalance、低 singlet coverage、有限 NTC reference 与 Reactome 排名 ties 都会限制 target-level robustness 和 pathway 解读。
"""
    (out / "methods/analysis_methods.md").write_text(analysis_methods, encoding="utf-8")
    changelog = f"""# REPORT V3 CHANGELOG

## Report restructuring

- DRAGEN run proof and detailed read-flow were moved to a collapsed reproducibility appendix.
- Target self-expression and the reused UMAP were moved to Exploratory analyses; large DE and Reactome tables remain companion data rather than the main narrative.
- Added a cell-assignment/target-usability funnel, explicit NTC QC, evidence-class definitions and class-threshold sensitivity, and a multi-guide robustness map.
- Added Reactome ranking-quality visualization and a secondary normalized mean-expression-difference ranking; primary signed rank-biserial ranking remains unchanged.
- Added descriptive Leiden-cluster composition review and conditional cluster-stratified effect sensitivity; no matched-cluster p/FDR is reported.

## Recomputed and preserved results

- Recomputed passing-cell partition, target coverage, target funnel, target-cell bins, guide balance, evidence-class counts, NTC counts, per-target summaries, and figures from source tables/MEX/H5AD.
- NTC: {ntc_summary['ntc_guide_features']} guides, {ntc_summary['active_ntc_guides_ge5_singlets']} active guides, {ntc_summary['ntc_singlet_cells']} singlet cells, {ntc_summary['ntc_feature_positive_cells_unique']} distinct feature-positive cells.
- Reactome v97 GMT name/ID fields were corrected in the V3 table/report labels; primary enrichment statistics were reused unchanged. Secondary ranking outcomes are provided side by side.
- Primary Mann–Whitney/tie-aware p, target-wise BH, rank-biserial effects, guide threshold, target threshold, and guide consistency rules were preserved.

## Conclusion changes

- The current evidence classes sum exactly to {len(cov)} observed targets: {", ".join(f"{k}={v}" for k, v in class_counts.items())}. SUPPORTED targets: {class_counts.get('SUPPORTED', 0)}; PARTIAL_SUPPORT: {class_counts.get('PARTIAL_SUPPORT', 0)}.
- V3 moves the interpretation to screen usability and evidence strength. Reactome and reused-embedding results are explicitly exploratory; target RNA is not used to infer editing outcome.
"""
    (out / "REPORT_V3_CHANGELOG.md").write_text(changelog, encoding="utf-8")
    report_summary = [
        "# V3 Analysis Summary", "", f"Sample: {project['sample_id']} / {project['project_id']}",
        f"Passing cells: {project['filtered_cells']:,}; singlet: {project['singlet_cells']:,}; multi-guide: {project['multi_guide_cells']:,}; unassigned: {project['unassigned_cells']:,}.",
        f"Observed targets: {len(cov):,}; target cells ≥20: {len(eligible20):,}; ≥2 active guides: {len(two_guides):,}; balanced: {len(balanced):,}; partial-or-better consistency: {len(consistent):,}.",
        f"Evidence classes: {dict(class_counts)}.", f"NTC: {ntc_summary}.",
        "Primary DE and guide-level effects reuse the verified existing calculations; V3 figures and summary tables are rebuilt."
    ]
    (out / "REANALYSIS_SUMMARY.md").write_text("\n".join(report_summary) + "\n", encoding="utf-8")

    # GSEA primary QC is reused; secondary QC is newly generated. Keep pathways labelled by stable Reactome ID.
    secondary_qc = read_tsv(tables / "gsea_ranking_sensitivity_qc.tsv") if (tables / "gsea_ranking_sensitivity_qc.tsv").exists() else pd.DataFrame(columns=["target_gene", "secondary_tie_fraction"])
    primary_qc["tie_fraction"] = pd.to_numeric(primary_qc["tie_fraction"], errors="coerce")
    secondary_qc["secondary_tie_fraction"] = pd.to_numeric(secondary_qc.get("secondary_tie_fraction"), errors="coerce")
    pathway_sensitivity = read_tsv(tables / "pathway/reactome_ranking_sensitivity.tsv.gz") if (tables / "pathway/reactome_ranking_sensitivity.tsv.gz").exists() else pd.DataFrame(columns=["primary_NES", "secondary_NES", "direction_concordant"])
    pathway_sensitivity["direction_concordant"] = pathway_sensitivity.get("direction_concordant", pd.Series(dtype=str)).astype(str).str.lower().isin(["true", "1"])
    figure_manifest = save_figures(out, project, funnel, cov_dist, conf, target_summary, representatives, primary_qc,
                                   secondary_qc, pathway_sensitivity, state_comp, umap_df, state_effect_summary)
    write_frame(tables / "figure_manifest.tsv", pd.DataFrame(figure_manifest))
    write_frame(tables / "claim_evidence_map.tsv", pd.DataFrame([
        ("Screen usability", "Passing-cell assignment and target coverage", "tables/screen_target_funnel.tsv; figures/fig01_screen_usability_funnel.png", "assignment classes partition passing cells"),
        ("NTC reference", "Filtered feature MEX, guide coverage and assignments", "tables/ntc_guide_qc.tsv; tables/project_summary.json", f"{ntc_summary['ntc_singlet_cells']} NTC singlets; target/control effects remain cell-level"),
        ("Target evidence class", "Target cells, guide count/balance, consistency and LOOGO", "tables/target_confidence.tsv; tables/target_threshold_sensitivity.tsv", "No combined score; no class is editing outcome"),
        ("Expression effects", "Reused full target-vs-NTC DE table and target summary", "tables/differential_expression/target_level_DE_all.tsv.gz; tables/target_level_summary.tsv", "One biological sample; cell-level statistics are exploratory"),
        ("Guide robustness", "Guide pairs and leave-one-guide-out source tables", "tables/guide_consistency_summary.tsv; tables/leave_one_guide_out_summary.tsv", "Guide imbalance is reported separately"),
        ("Reactome", "Primary signed rank-biserial and secondary mean-difference ranking", "tables/pathway/reactome_ranking_sensitivity.tsv.gz; tables/gsea_ranking_qc.tsv", "Exploratory; ranking ties and single-sample design limit interpretation"),
        ("Leiden cluster composition", "Unannotated input H5AD Leiden clusters and existing embedding", "tables/cell_state_composition.tsv; tables/cell_state_effect_sensitivity.tsv", "Descriptive/exploratory grouping labels; do not directly identify cell states or establish mechanism"),
    ], columns=["claim", "evidence_basis", "artifact", "interpretation_boundary"]))
    write_json(out / "05_metadata/v3_preparation.json", {"report_version": 3, "prepared": True,
               "primary_pathway_labels_corrected": corrected_rows, "secondary_gsea_rows": secondary_count,
               "h5ad_order_and_state": order_and_state, "ntc_qc": ntc_summary,
               "class_counts": {k: int(v) for k, v in class_counts.items()},
               "target_funnel": {r.stage: int(r.count) for r in funnel.itertuples()},
               "cell_state_effect_status_counts": dict(Counter(state_effect_summary.status)),
               "figure_count": len(figure_manifest), "analysis_version": "v3"})
    print(json.dumps({"output_root": str(out), "targets": len(cov), "eligible": len(target_summary),
                      "classes": dict(class_counts), "ntc": ntc_summary,
                      "target_funnel": {r.stage: int(r.count) for r in funnel.itertuples()},
                      "state_sensitivity": dict(Counter(state_effect_summary.status)),
                      "secondary_gsea_rows": secondary_count, "figures": len(figure_manifest)}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
