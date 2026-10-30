#!/usr/bin/env python3
"""Reanalyse a completed native DRAGEN PIPseq CRISPR run into v2 tables and figures."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.io import mmread


NPG = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
    "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
]
NPG_SIGNED = ["#3C5488", "#F7F7F7", "#E64B35"]


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=json_scalar) + "\n", encoding="utf-8")


def json_scalar(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_tsv(path: Path, columns: list[str] | tuple[str, ...], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_rows_gzip(path: Path, columns: list[str] | tuple[str, ...]):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=6)
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    return handle, writer


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required non-empty file missing: {path}")


def load_metric_index(path: Path) -> dict[tuple[str, str, str], tuple[str, str]]:
    result = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 4 and row[0] and row[2]:
                result[(row[0].strip(), row[1].strip(), row[2].strip())] = (
                    row[3].strip(), row[4].strip() if len(row) > 4 else ""
                )
    return result


def mvalue(index, section: str, name: str, sample: str) -> str:
    value = index.get((section, sample, name))
    if value is None:
        raise KeyError(f"DRAGEN metric absent: {section} / {sample} / {name}")
    return value[0]


def metric_int(index, section: str, name: str, sample: str) -> int:
    value = mvalue(index, section, name, sample)
    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"DRAGEN metric is not a count: {section}/{name}={value!r}")
    return int(value)


def metric_float(index, section: str, name: str, sample: str) -> float:
    return float(mvalue(index, section, name, sample))


def metric_percent(index, section: str, name: str, sample: str) -> float:
    value = index.get((section, sample, name))
    if value is None or len(value) < 2 or not value[1] or value[1].upper() == "NA":
        raise KeyError(f"DRAGEN metric percentage absent: {section} / {sample} / {name}")
    return float(value[1]) / 100.0


def pct(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else math.nan


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=np.float64)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    idx = np.flatnonzero(valid)
    order = np.argsort(p[idx], kind="mergesort")
    ranked = p[idx][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    dest = idx[order]
    out[dest] = adjusted
    return out


def finite_float(value) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def sparse_dense(matrix, row_indices: np.ndarray, col_start: int, col_end: int) -> np.ndarray:
    return matrix[col_start:col_end, :][:, row_indices].toarray().astype(np.float64, copy=False)


def rank_effect_block(norm_csc, raw_csc, target_idx: np.ndarray, control_idx: np.ndarray,
                      start: int, end: int, need_tests: bool = True) -> dict[str, np.ndarray]:
    xt = sparse_dense(norm_csc, target_idx, start, end)
    xc = sparse_dense(norm_csc, control_idx, start, end)
    ct = sparse_dense(raw_csc, target_idx, start, end)
    cc = sparse_dense(raw_csc, control_idx, start, end)
    nt, nc = len(target_idx), len(control_idx)
    if need_tests:
        test = stats.mannwhitneyu(
            xt, xc, axis=1, alternative="two-sided", method="asymptotic", use_continuity=True
        )
        u_statistic = np.asarray(test.statistic)
    else:
        # Guide concordance and leave-one-guide-out need only U, not a p-value.
        ranks = stats.rankdata(np.concatenate((xt, xc), axis=1), axis=1, method="average")
        u_statistic = ranks[:, :nt].sum(axis=1) - nt * (nt + 1) / 2.0
    result = {
        "rank_biserial": (2.0 * u_statistic / (nt * nc) - 1.0).astype(np.float32),
        "pct_target": (np.count_nonzero(ct, axis=1) / nt).astype(np.float32),
        "pct_control": (np.count_nonzero(cc, axis=1) / nc).astype(np.float32),
        "mean_target_log1p": xt.mean(axis=1),
        "mean_control_log1p": xc.mean(axis=1),
        "mean_target_norm": np.expm1(xt).mean(axis=1),
        "mean_control_norm": np.expm1(xc).mean(axis=1),
    }
    if need_tests:
        result["p_value"] = np.asarray(test.pvalue, dtype=np.float64)
    return result


def read_matrix(path: Path) -> sparse.csr_matrix:
    require_file(path)
    with gzip.open(path, "rb") as handle:
        matrix = mmread(handle)
    return matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)


def read_lines_gzip(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def stable_spearman(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3 or np.ptp(a[mask]) == 0 or np.ptp(b[mask]) == 0:
        return math.nan
    return float(stats.spearmanr(a[mask], b[mask]).statistic)


def gene_metadata(adata) -> tuple[list[str], list[str]]:
    ids = adata.var_names.astype(str).tolist()
    names = adata.var["gene_name"].astype(str).tolist()
    names = [n if n and n.lower() != "nan" else gid for gid, n in zip(ids, names)]
    return ids, names


def build_flow_contract(config, source, root, metrics, passing, total_unique, status_counts,
                        raw_gex_molecules, filtered_gene_molecules,
                        raw_feature_molecules, filtered_feature_molecules):
    sample = str(config["sample_id"])
    gex_reads = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Total input reads", sample)
    exact = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Reads with exactly matching barcodes", sample)
    corrected = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Reads with corrected barcodes", sample)
    nonmatching = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Reads with non-matching barcodes", sample)
    missing = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Reads missing barcodes", sample)
    barcoded = metric_int(metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Total barcoded reads", sample)
    if exact + corrected + nonmatching + missing != gex_reads or exact + corrected != barcoded:
        raise RuntimeError("GEX barcode read fates do not conserve their DRAGEN denominators")
    gex_fate_names = (
        "Unique exon matching reads", "Unique intron matching reads", "Filtered antisense reads",
        "Filtered ambiguously matching reads", "Filtered non-matching reads", "Mitochondrial reads",
    )
    gex_fates = [(name, metric_int(metrics, "SINGLE-CELL RNA READ MATCHING METRICS", name, sample))
                 for name in gex_fate_names]
    if sum(value for _, value in gex_fates) != barcoded:
        raise RuntimeError("GEX read-matching fates do not sum to total barcoded reads")

    feature_reads = metric_int(metrics, "SINGLE-CELL RNA FEATURE METRICS", "Total input feature reads", sample)
    feature_barcoded = metric_int(metrics, "SINGLE-CELL RNA FEATURE METRICS", "Total barcoded feature reads", sample)
    feature_invalid = metric_int(metrics, "SINGLE-CELL RNA FEATURE METRICS", "Feature reads with invalid barcodes", sample)
    feature_matched = metric_int(metrics, "SINGLE-CELL RNA FEATURE METRICS", "Feature matching reads", sample)
    feature_nonmatched = metric_int(metrics, "SINGLE-CELL RNA FEATURE METRICS", "Feature non-matching reads", sample)
    if feature_barcoded + feature_invalid != feature_reads or feature_matched + feature_nonmatched != feature_barcoded:
        raise RuntimeError("CRISPR feature read fates do not conserve their DRAGEN denominators")

    nodes, links, transforms = [], [], []

    def add_node(view, node, stage, label, value, unit, scope, kind):
        nodes.append({"view_id": view, "node_id": node, "stage": stage, "label": label,
                      "value": int(value), "unit": unit, "scope": scope, "node_kind": kind})

    def add_link(view, link, source_id, target_id, value, reason):
        links.append({"view_id": view, "link_id": link, "source": source_id,
                      "target": target_id, "value": int(value), "reason": reason})

    # GEX read fate chain. DRAGEN counts transcript R2 reads as single-ended reads.
    v = "gex_reads"
    add_node(v, "gex_input", 0, "GEX total input reads", gex_reads, "reads", "GEX library", "retained")
    for node, label, value, kind in (
        ("gex_exact", "Exact barcode match", exact, "branch"),
        ("gex_corrected", "Corrected barcode", corrected, "branch"),
        ("gex_no_barcode", "Missing barcode", missing, "loss"),
        ("gex_unmatched_barcode", "Non-matching barcode", nonmatching, "loss"),
    ):
        add_node(v, node, 1, label, value, "reads", "barcode classification", kind)
        add_link(v, f"input_to_{node}", "gex_input", node, value, "DRAGEN barcode read classification")
    add_node(v, "gex_barcoded", 2, "Reads with valid exact or corrected barcode", barcoded,
             "reads", "barcoded GEX reads", "retained")
    add_link(v, "exact_to_barcoded", "gex_exact", "gex_barcoded", exact, "exact barcode reads retained")
    add_link(v, "corrected_to_barcoded", "gex_corrected", "gex_barcoded", corrected, "corrected barcode reads retained")
    for i, (name, value) in enumerate(gex_fates):
        node = f"gex_map_{i}"
        retained = name in {"Unique exon matching reads", "Unique intron matching reads", "Mitochondrial reads"}
        add_node(v, node, 3, name.replace(" reads", ""), value, "reads", "DRAGEN primary read-matching fate",
                 "retained" if retained else "loss")
        add_link(v, f"barcoded_to_map_{i}", "gex_barcoded", node, value,
                 "mutually exclusive DRAGEN read-matching category")

    # CRISPR feature reads use a separate view because they are a separate library/read group.
    v = "crispr_reads"
    add_node(v, "crispr_input", 0, "CRISPR feature input reads", feature_reads, "reads", "CRISPR library", "retained")
    add_node(v, "crispr_barcoded", 1, "Valid CRISPR barcode and linker", feature_barcoded,
             "reads", "CRISPR barcode/linker classification", "retained")
    add_node(v, "crispr_invalid", 1, "Invalid barcode or CRISPR linker", feature_invalid,
             "reads", "CRISPR barcode/linker classification", "loss")
    add_link(v, "input_to_barcoded", "crispr_input", "crispr_barcoded", feature_barcoded,
             "DRAGEN feature reads with valid barcode/linker")
    add_link(v, "input_to_invalid", "crispr_input", "crispr_invalid", feature_invalid,
             "DRAGEN feature reads with invalid barcode or linker")
    for node, label, value, kind, reason in (
        ("crispr_matched", "Matched guide feature", feature_matched, "retained", "matched feature reference within DRAGEN allowance"),
        ("crispr_nonmatched", "No matching guide feature", feature_nonmatched, "loss", "no feature-reference match"),
    ):
        add_node(v, node, 2, label, value, "reads", "guide-reference matching", kind)
        add_link(v, f"barcoded_to_{node}", "crispr_barcoded", node, value, reason)

    # Cell fates: the complement is defined by DRAGEN's binary cell-calling result.
    v = "cells"
    not_passing = total_unique - passing
    if not_passing < 0:
        raise RuntimeError("passing cells exceed total unique barcodes")
    add_node(v, "all_barcodes", 0, "All unique GEX barcodes", total_unique, "cells", "DRAGEN barcode cell calling", "retained")
    add_node(v, "cells_passing", 1, "DRAGEN passing cells", passing, "cells", "DRAGEN cell-calling threshold", "retained")
    add_node(v, "barcodes_not_passing", 1, "Barcodes below cell-calling threshold", not_passing,
             "cells", "complement of DRAGEN passing cells", "loss")
    add_link(v, "barcodes_to_passing", "all_barcodes", "cells_passing", passing, "barcodes called as passing cells")
    add_link(v, "barcodes_to_not_passing", "all_barcodes", "barcodes_not_passing", not_passing,
             "all unique barcodes not called as passing cells")
    for node, label, key in (
        ("cells_singlet", "Single-guide assigned", "singlet"),
        ("cells_multi", "Multi-guide assigned", "multi"),
        ("cells_unassigned", "Unassigned guide", "none"),
    ):
        value = int(status_counts.get(key, 0))
        add_node(v, node, 2, label, value, "cells", "cell-level guide assignment", "branch")
        add_link(v, f"passing_to_{node}", "cells_passing", node, value, "partition by upstream guide-assignment status")

    # Independent molecule views come from raw and filtered MEX sums, not from read-count proxies.
    for view, label, total, filtered, retained_label, background_label, raw_scope in (
        ("gex_molecules", "GEX molecules", raw_gex_molecules, filtered_gene_molecules,
         "Molecules in passing cells", "Molecules in barcodes below threshold", "barcodes represented in raw GEX MEX"),
        ("feature_molecules", "CRISPR feature molecules", raw_feature_molecules, filtered_feature_molecules,
         "Feature molecules in passing cells", "Feature molecules in barcodes below threshold", "barcodes represented in raw feature MEX"),
    ):
        remainder = total - filtered
        if remainder < 0:
            raise RuntimeError(f"negative unfiltered molecule remainder for {view}")
        add_node(view, f"{view}_all", 0, f"All {label.lower()} in raw MEX", total, "molecules", raw_scope, "retained")
        add_node(view, f"{view}_filtered", 1, retained_label, filtered, "molecules", "DRAGEN filtered barcodes", "branch")
        add_node(view, f"{view}_background", 1, background_label, remainder, "molecules", "barcodes excluded by DRAGEN cell calling", "loss")
        add_link(view, f"all_to_filtered_{view}", f"{view}_all", f"{view}_filtered", filtered,
                 "sum of filtered MEX counts")
        add_link(view, f"all_to_background_{view}", f"{view}_all", f"{view}_background", remainder,
                 "all-barcode MEX sum minus filtered-barcode MEX sum")

    transforms.extend([
        {"transform_id": "gex_reads_to_molecules", "source_view_id": "gex_reads", "source_node_id": "gex_barcoded",
         "target_view_id": "gex_molecules", "target_node_id": "gex_molecules_all", "relation": "aggregate",
         "label": "GEX reads aggregated to molecules", "reason": "DRAGEN UMI/IMI molecule counting; units differ"},
        {"transform_id": "crispr_reads_to_molecules", "source_view_id": "crispr_reads", "source_node_id": "crispr_barcoded",
         "target_view_id": "feature_molecules", "target_node_id": "feature_molecules_all", "relation": "aggregate",
         "label": "Feature reads aggregated to molecules", "reason": "DRAGEN feature matching and molecule counting; units differ"},
        {"transform_id": "gex_molecules_to_cells", "source_view_id": "gex_molecules", "source_node_id": "gex_molecules_filtered",
         "target_view_id": "cells", "target_node_id": "cells_passing", "relation": "aggregate",
         "label": "GEX molecules assigned across passing cells", "reason": "one molecule count is distributed among passing cells"},
        {"transform_id": "feature_molecules_to_cells", "source_view_id": "feature_molecules", "source_node_id": "feature_molecules_filtered",
         "target_view_id": "cells", "target_node_id": "cells_passing", "relation": "aggregate",
         "label": "Feature molecules assigned across passing cells", "reason": "guide feature molecules are aggregated per passing cell"},
    ])
    views = [
        {"view_id": "gex_reads", "display_name": "GEX reads", "order": 0, "is_default": "true"},
        {"view_id": "crispr_reads", "display_name": "CRISPR feature reads", "order": 1, "is_default": "false"},
        {"view_id": "cells", "display_name": "Cells", "order": 2, "is_default": "false"},
        {"view_id": "gex_molecules", "display_name": "GEX molecules", "order": 3, "is_default": "false"},
        {"view_id": "feature_molecules", "display_name": "CRISPR feature molecules", "order": 4, "is_default": "false"},
    ]
    return nodes, links, transforms, views


def run(config: dict, config_path: Path) -> dict:
    sample = str(config["sample_id"])
    started_at = time.strftime("%Y-%m-%d %H:%M:%S %Z")
    source = Path(config["source_root"])
    root = Path(config["output_root"])
    if (root / "05_metadata/stage_complete.json").exists():
        raise RuntimeError(f"refusing to overwrite a completed reanalysis: {root}")
    runtime = Path(config["runtime_root"])
    for name in ("tables", "figures", "report", "methods", "logs", "05_metadata"):
        (root / name).mkdir(parents=True, exist_ok=True)
    for name in ("runtime", "cache", "matplotlib", "numba", "sankey"):
        (runtime / name).mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, root / "methods/reanalysis.v2.json")
    write_methods(root, config)
    os.environ.update({
        "TMPDIR": str(runtime / "runtime"), "TMP": str(runtime / "runtime"),
        "TEMP": str(runtime / "runtime"), "XDG_CACHE_HOME": str(runtime / "cache"),
        "MPLCONFIGDIR": str(runtime / "matplotlib"), "NUMBA_CACHE_DIR": str(runtime / "numba"),
        "OMP_NUM_THREADS": str(config["analysis_threads"]),
        "OPENBLAS_NUM_THREADS": str(config["analysis_threads"]),
        "MKL_NUM_THREADS": str(config["analysis_threads"]),
        "NUMBA_NUM_THREADS": str(config["analysis_threads"]),
    })
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})

    native_marker = source / "02_dragen_scrna_crispr/05_metadata/stage_complete.json"
    require_file(native_marker)
    marker = read_json(native_marker)
    if marker.get("status") != "complete":
        raise RuntimeError("source DRAGEN stage is not marked complete")
    dragen_dir = source / "02_dragen_scrna_crispr/04_tool_output/dragen"
    run_log = next(dragen_dir.glob("dragen_run_*.log"), None)
    if run_log is None:
        raise RuntimeError("native DRAGEN run log is absent")
    log_text = run_log.read_text(encoding="utf-8", errors="replace")
    if "PipSeq Software License Authentication Success" not in log_text or "DRAGEN finished normally" not in log_text:
        raise RuntimeError("source does not document native CRISPR license success and normal DRAGEN completion")
    metrics_path = dragen_dir / f"{sample}.scRNA_metrics.csv"
    mapping_path = dragen_dir / f"{sample}.mapping_metrics.csv"
    h5ad_path = source / "03_cell_guide_qc/04_tool_output" / f"{sample}_gex_filtered.h5ad"
    assignment_path = source / "03_cell_guide_qc/01_tables/cell_guide_assignments.tsv"
    guide_summary_path = source / "03_cell_guide_qc/01_tables/guide_summary.tsv"
    singlet_assignment_path = dragen_dir / f"{sample}.scRNA.singlet_cell_guide_assignments.csv"
    feature_reference_path = source / "01_input_preflight/04_tool_output/dragen_feature_reference.csv"
    for path in (metrics_path, mapping_path, h5ad_path, assignment_path, guide_summary_path,
                 singlet_assignment_path, feature_reference_path):
        require_file(path)
    sc_metrics = load_metric_index(metrics_path)
    map_metrics = load_metric_index(mapping_path)

    adata = ad.read_h5ad(h5ad_path)
    if "log_normalized" not in adata.layers:
        raise RuntimeError("GEX H5AD has no log_normalized layer")
    norm_layer = adata.layers["log_normalized"]
    norm_cell_gene = norm_layer.tocsr(copy=True) if sparse.issparse(norm_layer) else sparse.csr_matrix(norm_layer)
    if norm_cell_gene.shape != adata.shape:
        raise RuntimeError("H5AD log_normalized layer differs from H5AD dimensions")

    mex_features_path = dragen_dir / f"{sample}.scRNA.filtered.features.tsv.gz"
    mex_barcodes_path = dragen_dir / f"{sample}.scRNA.filtered.barcodes.tsv.gz"
    filtered_matrix_path = dragen_dir / f"{sample}.scRNA.filtered.matrix.mtx.gz"
    raw_features_path = dragen_dir / f"{sample}.scRNA.features.tsv.gz"
    raw_barcodes_path = dragen_dir / f"{sample}.scRNA.barcodes.tsv.gz"
    raw_matrix_path = dragen_dir / f"{sample}.scRNA.matrix.mtx.gz"
    for path in (mex_features_path, mex_barcodes_path, filtered_matrix_path, raw_features_path,
                 raw_barcodes_path, raw_matrix_path):
        require_file(path)
    feature_rows = [line.split("\t") for line in read_lines_gzip(mex_features_path)]
    raw_feature_rows = [line.split("\t") for line in read_lines_gzip(raw_features_path)]
    barcodes = read_lines_gzip(mex_barcodes_path)
    if barcodes != adata.obs_names.astype(str).tolist():
        raise RuntimeError("filtered MEX barcode order does not match the GEX H5AD")
    if [row[0] for row in feature_rows if len(row) > 2 and row[2] == "Gene Expression"] != adata.var_names.astype(str).tolist():
        raise RuntimeError("filtered MEX gene feature order does not match the GEX H5AD")
    if len(feature_rows) != len(raw_feature_rows):
        raise RuntimeError("raw and filtered MEX feature counts differ")
    if [(row[0], row[2] if len(row) > 2 else "") for row in feature_rows] != [
        (row[0], row[2] if len(row) > 2 else "") for row in raw_feature_rows
    ]:
        raise RuntimeError("raw and filtered MEX feature identities or order differ")
    guide_feature_indices = [i for i, row in enumerate(feature_rows) if len(row) > 2 and row[2] == "CRISPR Direct Capture"]
    gene_feature_indices = [i for i, row in enumerate(feature_rows) if len(row) > 2 and row[2] == "Gene Expression"]
    guide_feature_ids = [feature_rows[i][0] for i in guide_feature_indices]
    matrix_feature_x_cell = read_matrix(filtered_matrix_path)
    if matrix_feature_x_cell.shape != (len(feature_rows), len(barcodes)):
        raise RuntimeError(f"filtered MEX dimensions unexpected: {matrix_feature_x_cell.shape}")
    raw_cell_gene = matrix_feature_x_cell[gene_feature_indices, :].T.tocsr()
    if raw_cell_gene.shape != norm_cell_gene.shape:
        raise RuntimeError("filtered MEX GEX counts and H5AD normalized layer differ in shape")
    sample_raw_totals = np.asarray(raw_cell_gene.sum(axis=1)).ravel().astype(np.int64)
    if not np.array_equal(sample_raw_totals, pd.to_numeric(adata.obs["total_counts"]).to_numpy(dtype=np.int64)):
        raise RuntimeError("filtered MEX raw gene counts do not match per-cell total_counts")
    raw_csc = raw_cell_gene.T.tocsc()
    norm_csc = norm_cell_gene.T.tocsc()
    n_cells, n_genes = raw_cell_gene.shape
    cell_x_feature = matrix_feature_x_cell.T.tocsr()
    guide_csr = cell_x_feature[:, guide_feature_indices].tocsr()
    filtered_gene_molecules = int(matrix_feature_x_cell[gene_feature_indices, :].sum())
    filtered_feature_molecules = int(matrix_feature_x_cell[guide_feature_indices, :].sum())
    del matrix_feature_x_cell, cell_x_feature

    guide_reference = pd.read_csv(feature_reference_path, dtype=str).fillna("")
    guide_reference = guide_reference.loc[guide_reference["feature_type"].eq("CRISPR Direct Capture")].copy()
    guide_reference = guide_reference.drop_duplicates("id", keep="first").set_index("id")
    missing_reference = sorted(set(guide_feature_ids) - set(guide_reference.index.astype(str)))
    if missing_reference:
        raise RuntimeError(f"MEX guide features absent from feature reference: {missing_reference[:5]}")
    gene_ids, gene_names = gene_metadata(adata)
    gene_names_array = np.asarray(gene_names, dtype=object)
    gene_ids_array = np.asarray(gene_ids, dtype=object)
    control_gene = str(config["control_gene"])

    assignment = read_tsv(assignment_path)
    if len(assignment) != n_cells or set(assignment["barcode"]) != set(barcodes):
        raise RuntimeError("cell guide assignment table does not cover the filtered barcode set")
    assignment = assignment.set_index("barcode").loc[barcodes].reset_index()
    singlet_rows = []
    with singlet_assignment_path.open(encoding="utf-8", newline="") as handle:
        singlet_rows = list(csv.DictReader(handle))
    singlet_map = {row["cell_barcode"]: row["feature_call"] for row in singlet_rows if row.get("cell_barcode") and row.get("feature_call")}
    assignment_guide = dict(zip(assignment["barcode"], assignment["guide_ids"].astype(str)))
    for barcode, guide in singlet_map.items():
        if barcode not in assignment_guide or assignment_guide[barcode] != guide:
            raise RuntimeError(f"native singlet guide call disagrees with cell assignment table for {barcode}")
    if int(assignment["guide_status"].eq("singlet").sum()) != len(singlet_map):
        raise RuntimeError("native singlet assignment count differs from cell assignment table")
    cell_index = {barcode: i for i, barcode in enumerate(barcodes)}
    guide_cells: dict[str, np.ndarray] = defaultdict(list)
    target_cells: dict[str, list[int]] = defaultdict(list)
    for barcode, guide in singlet_map.items():
        idx = cell_index[barcode]
        guide_cells[guide].append(idx)
        target = str(guide_reference.loc[guide, "gene_name"])
        target_cells[target].append(idx)
    guide_cells = {key: np.asarray(value, dtype=np.int32) for key, value in guide_cells.items()}
    target_cells = {key: np.asarray(value, dtype=np.int32) for key, value in target_cells.items()}
    control_mask = assignment["guide_status"].eq("singlet").to_numpy() & assignment["target_genes"].eq(control_gene).to_numpy()
    control_idx = np.flatnonzero(control_mask).astype(np.int32)
    if len(control_idx) < int(config["min_control_cells"]):
        raise RuntimeError(f"control cells below configured minimum: {len(control_idx)}")
    observed_target_counts = {
        target: len(indices) for target, indices in target_cells.items() if target and target != control_gene
    }

    # Exact raw/filtered molecule tallies keep molecules separate from reads and cells.
    raw_matrix = read_matrix(raw_matrix_path)
    raw_barcodes = read_lines_gzip(raw_barcodes_path)
    if raw_matrix.shape != (len(raw_feature_rows), len(raw_barcodes)):
        raise RuntimeError(f"raw MEX dimensions {raw_matrix.shape} do not match its features/barcodes ({len(raw_feature_rows)}, {len(raw_barcodes)})")
    if len(raw_barcodes) != len(set(raw_barcodes)):
        raise RuntimeError("raw MEX barcode list contains duplicates")
    raw_type = [row[2] if len(row) > 2 else "" for row in raw_feature_rows]
    raw_gene_rows = [i for i, value in enumerate(raw_type) if value == "Gene Expression"]
    raw_guide_rows = [i for i, value in enumerate(raw_type) if value == "CRISPR Direct Capture"]
    raw_gex_molecules = int(raw_matrix[raw_gene_rows, :].sum())
    raw_feature_molecules = int(raw_matrix[raw_guide_rows, :].sum())
    reported_gex_molecules = metric_int(sc_metrics, "SINGLE-CELL RNA MOLECULE COUNT METRICS", "Total molecules", sample)
    reported_feature_molecules = metric_int(sc_metrics, "SINGLE-CELL RNA FEATURE METRICS", "Total feature molecules", sample)
    if raw_gex_molecules != reported_gex_molecules:
        raise RuntimeError(f"raw MEX GEX molecules {raw_gex_molecules} != DRAGEN {reported_gex_molecules}")
    if raw_feature_molecules != reported_feature_molecules:
        raise RuntimeError(f"raw MEX feature molecules {raw_feature_molecules} != DRAGEN {reported_feature_molecules}")
    if filtered_gene_molecules > raw_gex_molecules or filtered_feature_molecules > raw_feature_molecules:
        raise RuntimeError("filtered molecule matrix exceeds all-barcode molecule matrix")
    del raw_matrix

    # Project summary, cell-level feature UMI measurements, and guide coverage.
    guide_index = {guide: i for i, guide in enumerate(guide_feature_ids)}
    guide_umi_totals = np.asarray(guide_csr.sum(axis=0)).ravel().astype(np.int64)
    guide_positive_cells = np.asarray(guide_csr.getnnz(axis=0)).ravel().astype(np.int64)
    guide_medians = []
    for j in range(guide_csr.shape[1]):
        values = guide_csr[:, j].data
        guide_medians.append(float(np.median(values)) if len(values) else 0.0)
    cell_feature_rows = []
    top1 = np.zeros(n_cells, dtype=np.int32)
    top2 = np.zeros(n_cells, dtype=np.int32)
    n_detected_guides = np.diff(guide_csr.indptr).astype(np.int32)
    feature_umi_per_cell = np.asarray(guide_csr.sum(axis=1)).ravel().astype(np.int64)
    for i in range(n_cells):
        values = guide_csr.data[guide_csr.indptr[i]:guide_csr.indptr[i + 1]]
        if len(values):
            ordered = np.partition(values, max(0, len(values) - 2))
            top1[i] = int(ordered[-1])
            top2[i] = int(ordered[-2]) if len(ordered) > 1 else 0
        top1_fraction = float(top1[i] / feature_umi_per_cell[i]) if feature_umi_per_cell[i] else 0.0
        top_ratio = float(top1[i] / top2[i]) if top2[i] else (math.inf if top1[i] else math.nan)
        cell_feature_rows.append((
            barcodes[i], int(n_detected_guides[i]), int(feature_umi_per_cell[i]), int(top1[i]),
            int(top2[i]), top1_fraction, top_ratio,
        ))
    write_tsv(root / "tables/cell_guide_feature_umi.tsv.gz",
              ("barcode", "detected_guides", "total_feature_umi", "top1_umi", "top2_umi", "top1_fraction", "top1_top2_ratio"),
              cell_feature_rows)
    guide_rows = []
    raw_guide_summary = read_tsv(guide_summary_path).set_index("sgRNA_id")
    for j, guide in enumerate(guide_feature_ids):
        ref = guide_reference.loc[guide]
        assigned = int(len(guide_cells.get(guide, [])))
        positive = int(guide_positive_cells[j])
        guide_rows.append((
            guide, str(ref["gene_name"]), str(ref["sequence"]), positive, assigned,
            int(guide_umi_totals[j]), float(guide_medians[j]),
            int(raw_guide_summary.loc[guide, "positive_cells"]) if guide in raw_guide_summary.index else 0,
            int(raw_guide_summary.loc[guide, "singlet_cells"]) if guide in raw_guide_summary.index else 0,
            assigned >= int(config["guide_min_singlet_cells"]),
        ))
    guide_columns = ("sgRNA_id", "target_gene", "sequence", "positive_cells_from_matrix", "singlet_cells",
                    "total_feature_umi", "median_feature_umi_per_positive_cell", "upstream_positive_cells",
                    "upstream_singlet_cells", "active_for_effects")
    write_tsv(root / "tables/guide_coverage.tsv", guide_columns, guide_rows)
    active_guides = {
        row[0]: guide_cells[row[0]] for row in guide_rows
        if row[0] in guide_cells and len(guide_cells[row[0]]) >= int(config["guide_min_singlet_cells"])
    }
    if len(active_guides) == 0:
        raise RuntimeError("no guide reached the configured singlet-cell threshold")
    cell_rows = []
    for row, feature_row in zip(assignment.itertuples(index=False), cell_feature_rows):
        cell_rows.append(tuple(row) + tuple(feature_row[1:]))
    write_tsv(root / "tables/guide_assignment_all.tsv.gz",
              tuple(assignment.columns) + ("detected_guides", "total_feature_umi", "top1_umi", "top2_umi", "top1_fraction", "top1_top2_ratio"),
              cell_rows)

    # Counts and explicit denominators are regenerated from source artifacts.
    total_unique = metric_int(sc_metrics, "SINGLE-CELL RNA CELL METRICS", "Total unique barcodes", sample)
    passing = metric_int(sc_metrics, "SINGLE-CELL RNA CELL METRICS", "Passing cells", sample)
    status_counts = assignment["guide_status"].value_counts().to_dict()
    if sum(int(status_counts.get(k, 0)) for k in ("singlet", "multi", "none")) != passing:
        raise RuntimeError(f"guide assignment states do not partition DRAGEN passing cells: {status_counts}")
    if len(barcodes) != passing:
        raise RuntimeError("filtered MEX barcode count differs from DRAGEN passing cells")
    project_summary = {
        "sample_id": sample, "project_id": config["project_id"],
        "dragen_native_crispr": True, "dragen_version": "4.5.4-1",
        "dragen_license_success_in_source_log": True, "dragen_finished_normally": True,
        "perturbation_mode": config["perturbation_mode"],
        "gex_reads": metric_int(sc_metrics, "SINGLE-CELL RNA BARCODE READ METRICS", "Total input reads", sample),
        "crispr_feature_reads": metric_int(sc_metrics, "SINGLE-CELL RNA FEATURE METRICS", "Total input feature reads", sample),
        "filtered_cells": passing, "total_unique_barcodes": total_unique,
        "singlet_cells": int(status_counts.get("singlet", 0)),
        "multi_guide_cells": int(status_counts.get("multi", 0)),
        "unassigned_cells": int(status_counts.get("none", 0)),
        "observed_targets": len(observed_target_counts),
        "eligible_targets": sum(v >= int(config["min_target_cells"]) for v in observed_target_counts.values()),
        "ntc_cells": int(len(control_idx)), "guide_features": len(guide_feature_ids),
        "active_guides_ge5_singlets": len(active_guides), "gene_features": n_genes,
        "gex_molecules_all_barcodes": raw_gex_molecules,
        "gex_molecules_passing_cells": filtered_gene_molecules,
        "feature_molecules_all_barcodes": raw_feature_molecules,
        "feature_molecules_passing_cells": filtered_feature_molecules,
        "singlet_rate_of_filtered_cells": pct(int(status_counts.get("singlet", 0)), passing),
        "multi_guide_rate_of_filtered_cells": pct(int(status_counts.get("multi", 0)), passing),
        "unassigned_rate_of_filtered_cells": pct(int(status_counts.get("none", 0)), passing),
        "mapq40_rate_gex": metric_percent(map_metrics, "MAPPING/ALIGNING PER RG", "Reads with MAPQ [40:inf)", "GEX_1"),
        "q30_rate_gex": metric_percent(map_metrics, "MAPPING/ALIGNING PER RG", "Q30 bases", "GEX_1"),
        "median_molecules_per_cell": metric_int(sc_metrics, "SINGLE-CELL RNA CELL METRICS", "Median molecules per passing cell", sample),
        "median_genes_per_cell": metric_int(sc_metrics, "SINGLE-CELL RNA CELL METRICS", "Median genes per passing cell", sample),
        "sequencing_saturation": metric_float(sc_metrics, "SINGLE-CELL RNA MOLECULE COUNT METRICS", "Sequencing saturation", sample),
        "fraction_reads_in_passing_cells": metric_float(sc_metrics, "SINGLE-CELL RNA CELL METRICS", "Fraction of reads in passing cells", sample),
        "feature_fraction_reads_in_passing_cells": metric_float(sc_metrics, "SINGLE-CELL RNA FEATURE METRICS", "Fraction of feature reads in passing cells", sample),
    }
    write_json(root / "tables/project_summary.json", project_summary)

    # All observed target counts and library/guide coverage thresholds.
    coverage_rows = []
    for threshold in (5, 10, 20, 30, 50):
        n = sum(value >= threshold for value in observed_target_counts.values())
        coverage_rows.append((threshold, n, len(observed_target_counts), pct(n, len(observed_target_counts))))
    write_tsv(root / "tables/target_coverage_thresholds.tsv",
              ("minimum_singlet_cells", "targets_meeting_threshold", "observed_targets", "fraction_targets"), coverage_rows)
    target_rows = []
    guides_by_target: dict[str, list[str]] = defaultdict(list)
    for guide, idx in active_guides.items():
        guides_by_target[str(guide_reference.loc[guide, "gene_name"])].append(guide)
    for target, count in sorted(observed_target_counts.items(), key=lambda item: (-item[1], item[0])):
        guide_ids = guides_by_target.get(target, [])
        guide_counts = [len(active_guides[g]) for g in guide_ids]
        target_rows.append((target, count, int(count >= int(config["min_target_cells"])),
                            len(guide_ids), sum(guide_counts),
                            max(guide_counts) / count if guide_counts and count else math.nan,
                            ";".join(sorted(guide_ids))))
    write_tsv(root / "tables/target_guide_coverage.tsv",
              ("target_gene", "singlet_cells", "eligible_for_DE", "active_guides_ge5", "active_guide_singlet_cells",
               "max_active_guide_fraction", "active_sgRNA_ids"), target_rows)

    # Build explicit data-flow contracts only after source counts have reconciled.
    flow_nodes, flow_links, flow_transforms, flow_views = build_flow_contract(
        config, source, root, sc_metrics, passing, total_unique, status_counts,
        raw_gex_molecules, filtered_gene_molecules, raw_feature_molecules, filtered_feature_molecules,
    )
    flow_dir = runtime / "sankey/output"
    flow_dir.mkdir(parents=True, exist_ok=True)
    node_columns = ("view_id", "node_id", "stage", "label", "value", "unit", "scope", "node_kind")
    link_columns = ("view_id", "link_id", "source", "target", "value", "reason")
    view_columns = ("view_id", "display_name", "order", "is_default")
    transform_columns = ("transform_id", "source_view_id", "source_node_id", "target_view_id", "target_node_id", "relation", "label", "reason")
    write_tsv(flow_dir / "read_flow_nodes.tsv", node_columns,
              [tuple(row[key] for key in node_columns) for row in flow_nodes])
    write_tsv(flow_dir / "read_flow_links.tsv", link_columns,
              [tuple(row[key] for key in link_columns) for row in flow_links])
    write_tsv(flow_dir / "read_flow_views.tsv", view_columns,
              [tuple(row[key] for key in view_columns) for row in flow_views])
    write_tsv(flow_dir / "read_flow_transformations.tsv",
              transform_columns, [tuple(row[key] for key in transform_columns) for row in flow_transforms])
    sankey_args = [sys.executable, "-m", "scperturb.reporting.read_flow", "--nodes", str(flow_dir / "read_flow_nodes.tsv"),
                   "--links", str(flow_dir / "read_flow_links.tsv"), "--views", str(flow_dir / "read_flow_views.tsv"),
                   "--transformations", str(flow_dir / "read_flow_transformations.tsv"), "--default-mode", "view",
                   "--node-layout", "main_top_losses", "--orientation", "horizontal", "--output-dir", str(flow_dir),
                   "--title", "scPerturb DRAGEN 原生流程：reads、molecules 与 cells"]
    plotly_asset = Path(config.get("pipeline_assets", Path(__file__).resolve().parents[2] / "assets")) / "plotly-3.5.1.min.js"
    if plotly_asset.is_file():
        sankey_args += ["--plotly-js", str(plotly_asset)]
    with (root / "logs/read_flow_sankey.log").open("w", encoding="utf-8") as log_handle:
        subprocess.run(sankey_args, check=True, stdout=log_handle, stderr=subprocess.STDOUT)
    subprocess.run([sys.executable, "-m", "scperturb.validation.read_flow", "--output-dir", str(flow_dir)], check=True)
    # Keep machine-readable flow inputs, payload, validation, and SVGs; the standalone
    # module HTML is only an intermediate because the final report is the sole report.
    shutil.copytree(flow_dir / "tables", root / "tables/read_flow", dirs_exist_ok=True)
    shutil.copytree(flow_dir / "figures/read_flow_sankey", root / "figures/read_flow_sankey", dirs_exist_ok=True)
    shutil.copy2(flow_dir / "reports/sankey_validation.json", root / "05_metadata/sankey_validation.json")
    write_sankey_delivery_manifest(flow_dir / "reports/sankey_manifest.tsv", root)

    # Differential expression and guide-level rank effects follow below.
    result = run_expression_analyses(config, root, adata, raw_csc, norm_csc, gene_ids, gene_names,
                                     control_idx, target_cells, guide_cells,
                                     active_guides, guide_reference, control_gene)
    result.update(project_summary)
    result["flow_contract"] = str(root / "tables/read_flow/read_flow_nodes.tsv")
    result["gseapy_version"] = gseapy.__version__
    result["started_at"] = started_at
    write_json(root / "tables/reanalysis_summary.json", result)
    write_figures(config, root, adata, result)
    write_claim_evidence_map(root, result)
    build_report_cmd = [sys.executable, "-m", "scperturb.reporting.intermediate",
                        "--config", str(config_path), "--summary", str(root / "tables/reanalysis_summary.json")]
    subprocess.run(build_report_cmd, check=True)
    write_json(root / "05_metadata/report_validation.json", {"status": "PENDING"})
    write_summary_docs(config, root, result, report_validation="PENDING")
    write_delivery_manifest(root)
    validate_report(config, root)
    write_summary_docs(config, root, result, report_validation="PASS")
    completion = {
        "status": "complete", "sample_id": sample, "stage": "scperturb_reanalysis_v2",
        "completed_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "source_native_dragen_stage": str(native_marker),
        "outputs": [str(p) for p in sorted(root.rglob("*")) if p.is_file()],
        "validation": {"flow": "PASS", "report": "PASS", "tables": "PASS"},
    }
    write_json(root / "05_metadata/stage_complete.json", completion)
    return result


DE_COLUMNS = (
    "target_gene", "target_cells", "control_cells", "gene_id", "gene_name",
    "mean_target_log1p", "mean_control_log1p", "mean_target_normalized", "mean_control_normalized",
    "log2fc", "pct_target", "pct_control", "delta_pct", "rank_biserial", "p_value", "fdr",
)
GUIDE_EFFECT_COLUMNS = (
    "sgRNA_id", "target_gene", "guide_cells", "control_cells", "gene_id", "gene_name",
    "mean_guide_log1p", "mean_control_log1p", "mean_guide_normalized", "mean_control_normalized",
    "log2fc", "pct_guide", "pct_control", "delta_pct", "rank_biserial",
)


def _stats_for_comparison(norm_csc, raw_csc, target_idx, control_idx, n_genes, block_size,
                          target_label, gene_ids, gene_names, need_p=True):
    parts = defaultdict(list)
    for start in range(0, n_genes, block_size):
        end = min(start + block_size, n_genes)
        values = rank_effect_block(norm_csc, raw_csc, target_idx, control_idx, start, end, need_tests=need_p)
        for key, value in values.items():
            parts[key].append(np.asarray(value))
    combined = {key: np.concatenate(value) for key, value in parts.items()}
    if need_p:
        combined["p_value"] = np.maximum(combined["p_value"], 1e-300)
        combined["fdr"] = bh_adjust(combined["p_value"])
    combined["gene_id"] = np.asarray(gene_ids, dtype=object)
    combined["gene_name"] = np.asarray(gene_names, dtype=object)
    return combined


def _load_effect_cache(path: Path, group_column: str, gene_column: str, effect_column: str,
                       pct_columns: tuple[str, ...], expected_groups, gene_ids):
    """Load a complete per-group effect table without retaining its wide source table."""
    n_genes = len(gene_ids)
    columns = [group_column, gene_column, effect_column, *pct_columns]
    cached: dict[str, tuple[np.ndarray, ...]] = {}
    row_count = 0
    chunksize = n_genes * 4
    for chunk in pd.read_csv(path, sep="\t", usecols=columns, dtype=str,
                             chunksize=chunksize, compression="infer"):
        for group, part in chunk.groupby(group_column, sort=False, dropna=False):
            group = str(group)
            if group in cached:
                raise ValueError(f"effect cache group is repeated across chunks: {group}")
            if len(part) != n_genes or not np.array_equal(part[gene_column].astype(str).to_numpy(),
                                                          np.asarray(gene_ids, dtype=str)):
                raise ValueError(f"effect cache gene rows are incomplete or out of order: {group}")
            arrays = [pd.to_numeric(part[col], errors="coerce").to_numpy(dtype=np.float32)
                      for col in (effect_column, *pct_columns)]
            cached[group] = tuple(arrays)
            row_count += len(part)
    if set(cached) != set(map(str, expected_groups)):
        missing = sorted(set(map(str, expected_groups)) - set(cached))
        extra = sorted(set(cached) - set(map(str, expected_groups)))
        raise ValueError(f"effect cache groups differ from expected groups; missing={missing[:5]}, extra={extra[:5]}")
    return cached, row_count


def _top_indices(effect: np.ndarray, pct_a: np.ndarray, pct_b: np.ndarray, threshold: float, top_n: int) -> np.ndarray:
    universe = np.flatnonzero((pct_a >= threshold) | (pct_b >= threshold))
    if not len(universe):
        universe = np.arange(len(effect))
    order = np.argsort(-np.abs(effect[universe]), kind="mergesort")
    return universe[order[:min(top_n, len(order))]]


def run_expression_analyses(config, root, adata, raw_csc, norm_csc, gene_ids, gene_names,
                            control_idx, target_cell_map, guide_cell_map, active_guides,
                            guide_reference, control_gene):
    sample = str(config["sample_id"])
    n_genes = len(gene_ids)
    block = int(config["gene_block_size"])
    min_target = int(config["min_target_cells"])
    min_guide = int(config["guide_min_singlet_cells"])
    min_pct = float(config["min_expression_fraction"])
    top_n = int(config["top_genes_for_consistency"])
    observed_target_counts = {
        target: len(indices) for target, indices in target_cell_map.items()
        if target and target != control_gene
    }
    gene_names_array = np.asarray(gene_names, dtype=object)
    eligible = sorted(
        (target for target, indices in target_cell_map.items()
         if target != control_gene and len(indices) >= min_target and len(control_idx) >= int(config["min_control_cells"])),
        key=lambda target: (-len(target_cell_map[target]), target),
    )
    if not eligible:
        raise RuntimeError("no target meets the configured target and NTC cell thresholds")

    de_path = root / "tables/differential_expression/target_level_DE_all.tsv.gz"
    de_effects: dict[str, np.ndarray] = {}
    de_pct: dict[str, np.ndarray] = {}
    de_pct_control: dict[str, np.ndarray] = {}
    de_rows = 0
    de_top_rows = []
    de_summary_rows = []
    try:
        de_cache, de_rows = _load_effect_cache(
            de_path, "target_gene", "gene_id", "rank_biserial", ("pct_target", "pct_control"),
            eligible, gene_ids,
        )
        for target, values in de_cache.items():
            de_effects[target], de_pct[target], de_pct_control[target] = values
        require_file(root / "tables/differential_expression/target_level_DE_top100.tsv")
        require_file(root / "tables/differential_expression/target_DE_summary.tsv")
        print(f"CACHE_HIT\ttarget_DE\tgroups={len(de_cache)}\trows={de_rows}", flush=True)
        reuse_de = True
    except Exception as exc:
        print(f"CACHE_MISS\ttarget_DE\t{type(exc).__name__}: {exc}", flush=True)
        reuse_de = False
    if not reuse_de:
        de_handle, de_writer = write_rows_gzip(de_path, DE_COLUMNS)
        try:
            for k, target in enumerate(eligible, start=1):
                indices = target_cell_map[target]
                result = _stats_for_comparison(norm_csc, raw_csc, indices, control_idx, n_genes, block,
                                               target, gene_ids, gene_names, need_p=True)
                effect = result["rank_biserial"].astype(np.float32)
                pct_target = result["pct_target"].astype(np.float32)
                pct_control = result["pct_control"].astype(np.float32)
                log2fc = np.log2((result["mean_target_norm"] + 0.1) / (result["mean_control_norm"] + 0.1))
                delta = pct_target - pct_control
                for i in range(n_genes):
                    de_writer.writerow((target, len(indices), len(control_idx), gene_ids[i], gene_names[i],
                                        result["mean_target_log1p"][i], result["mean_control_log1p"][i],
                                        result["mean_target_norm"][i], result["mean_control_norm"][i],
                                        log2fc[i], pct_target[i], pct_control[i], delta[i], effect[i],
                                        result["p_value"][i], result["fdr"][i]))
                de_rows += n_genes
                de_effects[target] = effect
                de_pct[target] = pct_target
                de_pct_control[target] = pct_control
                rank_order = np.argsort(-np.abs(effect), kind="mergesort")[:min(100, n_genes)]
                for i in rank_order:
                    de_top_rows.append((target, len(indices), len(control_idx), gene_ids[i], gene_names[i],
                                        result["mean_target_log1p"][i], result["mean_control_log1p"][i],
                                        result["mean_target_norm"][i], result["mean_control_norm"][i],
                                        log2fc[i], pct_target[i], pct_control[i], delta[i], effect[i],
                                        result["p_value"][i], result["fdr"][i]))
                de_summary_rows.append((target, len(indices), len(control_idx), n_genes,
                                        int(np.sum(result["p_value"] < 0.05)),
                                        int(np.sum(result["fdr"] < 0.05)),
                                        float(np.median(pct_target)),
                                        float(np.median(np.abs(effect))),
                                        float(np.median(result["fdr"]))))
                print(f"TARGET_DE\t{k}/{len(eligible)}\t{target}\tcells={len(indices)}\tgenes={n_genes}", flush=True)
        finally:
            de_handle.close()
        write_tsv(root / "tables/differential_expression/target_level_DE_top100.tsv", DE_COLUMNS, de_top_rows)
        write_tsv(root / "tables/differential_expression/target_DE_summary.tsv",
                  ("target_gene", "target_cells", "control_cells", "genes_tested", "genes_p_lt_0_05",
                   "genes_fdr_lt_0_05", "median_pct_target", "median_abs_rank_biserial", "median_fdr"), de_summary_rows)

    # Per-guide effect sizes are full-gene but deliberately have no per-guide p-values/FDR.
    guide_path = root / "tables/guide_level_effects.tsv.gz"
    guide_effects: dict[str, np.ndarray] = {}
    guide_pcts: dict[str, np.ndarray] = {}
    active_order = sorted(active_guides)
    try:
        guide_cache, guide_rows = _load_effect_cache(
            guide_path, "sgRNA_id", "gene_id", "rank_biserial", ("pct_guide",), active_order, gene_ids,
        )
        for guide, values in guide_cache.items():
            guide_effects[guide], guide_pcts[guide] = values
        if guide_rows != len(active_order) * n_genes:
            raise ValueError(f"guide effect row count differs: {guide_rows}")
        print(f"CACHE_HIT\tguide_effects\tgroups={len(guide_cache)}\trows={guide_rows}", flush=True)
        reuse_guides = True
    except Exception as exc:
        print(f"CACHE_MISS\tguide_effects\t{type(exc).__name__}: {exc}", flush=True)
        reuse_guides = False
    if not reuse_guides:
        guide_handle, guide_writer = write_rows_gzip(guide_path, GUIDE_EFFECT_COLUMNS)
        try:
            for k, guide in enumerate(active_order, start=1):
                indices = guide_cell_map[guide]
                result = _stats_for_comparison(norm_csc, raw_csc, indices, control_idx, n_genes, block,
                                               guide, gene_ids, gene_names, need_p=False)
                effect = result["rank_biserial"].astype(np.float32)
                pct_guide = result["pct_target"].astype(np.float32)
                pct_control = result["pct_control"].astype(np.float32)
                log2fc = np.log2((result["mean_target_norm"] + 0.1) / (result["mean_control_norm"] + 0.1))
                target = str(guide_reference.loc[guide, "gene_name"])
                for start in range(0, n_genes, block):
                    end = min(start + block, n_genes)
                    for i in range(start, end):
                        guide_writer.writerow((guide, target, len(indices), len(control_idx), gene_ids[i], gene_names[i],
                                               result["mean_target_log1p"][i], result["mean_control_log1p"][i],
                                               result["mean_target_norm"][i], result["mean_control_norm"][i],
                                               log2fc[i], pct_guide[i], pct_control[i], pct_guide[i] - pct_control[i], effect[i]))
                guide_effects[guide] = effect
                guide_pcts[guide] = pct_guide
                print(f"GUIDE_EFFECTS\t{k}/{len(active_order)}\t{guide}\tcells={len(indices)}", flush=True)
        finally:
            guide_handle.close()

    guides_by_target: dict[str, list[str]] = defaultdict(list)
    for guide in active_order:
        guides_by_target[str(guide_reference.loc[guide, "gene_name"])].append(guide)
    pair_rows = []
    consistency_rows = []
    top_sets: dict[str, np.ndarray] = {}
    for target in sorted(target_cell_map):
        if target == control_gene:
            continue
        guides = sorted(guides_by_target.get(target, []))
        pair_metrics = []
        for i in range(len(guides)):
            for j in range(i + 1, len(guides)):
                a, b = guides[i], guides[j]
                ea, eb = guide_effects[a], guide_effects[b]
                pa, pb = guide_pcts[a], guide_pcts[b]
                common = (pa >= min_pct) | (pb >= min_pct)
                rho = stable_spearman(ea[common], eb[common])
                top_a = _top_indices(ea, pa, np.zeros(n_genes, dtype=np.float32), min_pct, top_n)
                top_b = _top_indices(eb, pb, np.zeros(n_genes, dtype=np.float32), min_pct, top_n)
                set_a, set_b = set(map(int, top_a)), set(map(int, top_b))
                union, intersection = set_a | set_b, set_a & set_b
                jaccard = len(intersection) / len(union) if union else math.nan
                union_idx = np.fromiter(union, dtype=np.int64)
                signs_ok = np.sign(ea[union_idx]) == np.sign(eb[union_idx]) if len(union_idx) else np.array([], dtype=bool)
                nonzero = (ea[union_idx] != 0) | (eb[union_idx] != 0) if len(union_idx) else np.array([], dtype=bool)
                concordance = float(signs_ok[nonzero].mean()) if nonzero.any() else math.nan
                pair_rows.append((target, a, b, len(guide_cell_map[a]), len(guide_cell_map[b]), rho,
                                  len(intersection), len(union), jaccard, concordance))
                pair_metrics.append((rho, jaccard, concordance))
        active_counts = [len(guide_cell_map[g]) for g in guides]
        total_target_cells = len(target_cell_map[target])
        max_fraction = max(active_counts) / total_target_cells if active_counts and total_target_cells else math.nan
        rhos = [x[0] for x in pair_metrics if math.isfinite(x[0])]
        overlaps = [x[1] for x in pair_metrics if math.isfinite(x[1])]
        concordances = [x[2] for x in pair_metrics if math.isfinite(x[2])]
        consistency_rows.append((target, total_target_cells, len(guides), max_fraction,
                                 len(pair_metrics), float(np.median(rhos)) if rhos else math.nan,
                                 float(np.median(overlaps)) if overlaps else math.nan,
                                 float(np.median(concordances)) if concordances else math.nan,
                                 "BALANCED" if math.isfinite(max_fraction) and max_fraction <= float(config["max_guide_fraction"]) else "GUIDE_DOMINATED"))
    write_tsv(root / "tables/guide_consistency_pairs.tsv",
              ("target_gene", "guide_a", "guide_b", "guide_a_cells", "guide_b_cells", "spearman_rho",
               "top_gene_overlap", "top_gene_union", "top100_jaccard", "directional_concordance"), pair_rows)
    write_tsv(root / "tables/guide_consistency_summary.tsv",
              ("target_gene", "target_cells", "active_guides_ge5", "max_guide_fraction", "guide_pairs",
               "median_pairwise_spearman", "median_top100_jaccard", "median_directional_concordance", "guide_balance_status"),
              consistency_rows)
    consistency_by_target = {row[0]: row for row in consistency_rows}

    # Leave one active guide out from the pooled target cells and compare its effect vector to the full target.
    loo_rows = []
    loo_by_target: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for target, target_idx in sorted(target_cell_map.items()):
        if target == control_gene:
            continue
        guides = sorted(guides_by_target.get(target, []))
        if len(guides) < 2 or len(target_idx) < 3:
            continue
        pooled = de_effects.get(target)
        pooled_pct = de_pct.get(target)
        if pooled is None:
            pooled_result = _stats_for_comparison(norm_csc, raw_csc, target_idx, control_idx, n_genes, block,
                                                   target, gene_ids, gene_names, need_p=False)
            pooled = pooled_result["rank_biserial"].astype(np.float32)
            pooled_pct = pooled_result["pct_target"].astype(np.float32)
        top_idx = _top_indices(pooled, pooled_pct, np.zeros(n_genes, dtype=np.float32), min_pct, top_n)
        for omit in guides:
            omit_idx = guide_cell_map[omit]
            keep_mask = ~np.isin(target_idx, omit_idx, assume_unique=False)
            kept = target_idx[keep_mask]
            if len(kept) < 3:
                rho, retention = math.nan, math.nan
            else:
                loo = _stats_for_comparison(norm_csc, raw_csc, kept, control_idx, n_genes, block,
                                            target, gene_ids, gene_names, need_p=False)["rank_biserial"].astype(np.float32)
                rho = stable_spearman(pooled, loo)
                valid = (pooled[top_idx] != 0) | (loo[top_idx] != 0)
                retention = float((np.sign(pooled[top_idx][valid]) == np.sign(loo[top_idx][valid])).mean()) if valid.any() else math.nan
            stable = bool(math.isfinite(rho) and rho >= float(config["loogo_spearman_min"])
                          and math.isfinite(retention) and retention >= float(config["loogo_direction_retention_min"]))
            loo_rows.append((target, omit, len(target_idx), len(kept), rho, retention, stable))
            loo_by_target[target].append((rho, retention))
        print(f"LOOGO\t{target}\tguides={len(guides)}", flush=True)
    write_tsv(root / "tables/leave_one_guide_out_summary.tsv",
              ("target_gene", "omitted_sgRNA_id", "pooled_target_cells", "remaining_target_cells",
               "spearman_vs_pooled", "top100_direction_retention", "stable"), loo_rows)

    # Target self-expression is an evidence field, never a knockout pass/fail when modality is unknown.
    self_rows = []
    gene_indices_by_name: dict[str, list[int]] = defaultdict(list)
    for i, name in enumerate(gene_names):
        gene_indices_by_name[str(name)].append(i)
    for target in sorted(target_cell_map):
        if target == control_gene:
            continue
        candidates = gene_indices_by_name.get(target, [])
        if not candidates:
            self_rows.append((target, len(target_cell_map[target]), math.nan, math.nan, math.nan, math.nan, "target symbol absent from feature matrix"))
            continue
        best = max(candidates, key=lambda i: abs(float(de_effects.get(target, np.zeros(n_genes))[i])) if target in de_effects else abs(float(target_cell_map[target].size)))
        if target in de_effects:
            self_rows.append((target, len(target_cell_map[target]), gene_ids[best], gene_names[best],
                              float(de_effects[target][best]), float(de_pct[target][best]),
                              "target-level signed rank-biserial effect; exploratory"))
        else:
            self_rows.append((target, len(target_cell_map[target]), gene_ids[best], gene_names[best],
                              math.nan, math.nan, "below target DE eligibility threshold"))
    write_tsv(root / "tables/target_self_expression_effect.tsv",
              ("target_gene", "target_cells", "gene_id", "gene_name", "rank_biserial", "pct_target", "interpretation"), self_rows)

    # Confidence uses evidence gates, not a combined numerical score.
    confidence_rows = []
    for target, n_target in sorted(observed_target_counts.items()):
        c = consistency_by_target.get(target)
        n_active = int(c[2]) if c else 0
        max_fraction = float(c[3]) if c else math.nan
        rho = float(c[5]) if c else math.nan
        direction = float(c[7]) if c else math.nan
        loo_metrics = loo_by_target.get(target, [])
        loo_rhos = [x[0] for x in loo_metrics if math.isfinite(x[0])]
        loo_directions = [x[1] for x in loo_metrics if math.isfinite(x[1])]
        loo_stable = bool(loo_metrics and len(loo_rhos) == len(loo_metrics) and len(loo_directions) == len(loo_metrics)
                          and min(loo_rhos) >= float(config["loogo_spearman_min"])
                          and min(loo_directions) >= float(config["loogo_direction_retention_min"]))
        if n_target < min_target:
            status = "LOW_CELL"
        elif n_active < 2:
            status = "LOW_GUIDE_SUPPORT"
        elif math.isfinite(max_fraction) and max_fraction > float(config["max_guide_fraction"]):
            status = "GUIDE_DOMINATED"
        elif (math.isfinite(rho) and rho < 0) or (math.isfinite(direction) and direction < float(config["direction_discordant_below"])):
            status = "DISCORDANT"
        elif (math.isfinite(rho) and rho >= float(config["spearman_supported_at_least"])
              and math.isfinite(direction) and direction >= float(config["direction_supported_at_least"])
              and loo_stable):
            status = "SUPPORTED"
        elif math.isfinite(rho) and rho >= 0 and math.isfinite(direction) and direction >= float(config["direction_discordant_below"]):
            status = "PARTIAL_SUPPORT"
        else:
            status = "INCONCLUSIVE"
        confidence_rows.append((target, n_target, n_active, max_fraction, rho, direction,
                                min(loo_rhos) if loo_rhos else math.nan,
                                min(loo_directions) if loo_directions else math.nan,
                                loo_stable, status,
                                "target self-expression is reported separately; mode UNKNOWN"))
    confidence_columns = ("target_gene", "target_cells", "active_guides_ge5", "max_guide_fraction",
                          "median_pairwise_spearman", "median_directional_concordance",
                          "minimum_loogo_spearman", "minimum_loogo_direction_retention", "loogo_stable",
                          "confidence_class", "self_effect_interpretation")
    write_tsv(root / "tables/target_confidence.tsv", confidence_columns, confidence_rows)

    similarity_targets = eligible
    if len(similarity_targets) >= 2:
        common = np.zeros(n_genes, dtype=bool)
        for target in similarity_targets:
            common |= (de_pct[target] >= min_pct)
        common |= np.asarray((raw_csc[:, control_idx].getnnz(axis=1)).ravel() / len(control_idx) >= min_pct)
        sub = np.vstack([de_effects[t][common] for t in similarity_targets]).astype(np.float64)
        ranked = np.vstack([stats.rankdata(row, method="average") for row in sub])
        corr = np.corrcoef(ranked)
        corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
        np.fill_diagonal(corr, 1.0)
        try:
            order = leaves_list(linkage(1.0 - np.clip(corr, -1, 1), method="average", optimal_ordering=True))
        except Exception:
            order = np.arange(len(similarity_targets))
        sim_path = root / "tables/target_effect_similarity.tsv"
        write_tsv(sim_path, ("target_gene", *similarity_targets),
                  [(similarity_targets[i], *corr[i].tolist()) for i in range(len(similarity_targets))])
        write_json(root / "tables/target_effect_similarity_order.json",
                   {"targets": similarity_targets, "cluster_order": [similarity_targets[int(i)] for i in order],
                    "correlation": "Spearman correlation of target-level rank-biserial effects over genes expressed in any eligible target/control at >=5%"})
    else:
        corr = np.eye(len(similarity_targets))
        order = np.arange(len(similarity_targets))
        write_tsv(root / "tables/target_effect_similarity.tsv", ("target_gene", *similarity_targets),
                  [(target, 1.0) for target in similarity_targets])

    gsea_result = run_reactome(config, root, eligible, target_cell_map, control_idx,
                               de_effects, de_pct, de_pct_control, gene_ids, gene_names, confidence_rows)
    return {
        "eligible_targets_de": len(eligible), "target_de_rows": de_rows,
        "active_guides_de": len(active_order), "guide_effect_rows": len(active_order) * n_genes,
        "guide_pairs": len(pair_rows), "loogo_rows": len(loo_rows),
        "confidence_counts": dict(Counter(row[-2] for row in confidence_rows)),
        "supported_targets": [row[0] for row in confidence_rows if row[-2] == "SUPPORTED"],
        "median_guide_pair_spearman": float(np.nanmedian([r[5] for r in consistency_rows])) if consistency_rows else math.nan,
        "median_guide_direction_concordance": float(np.nanmedian([r[7] for r in consistency_rows])) if consistency_rows else math.nan,
        "target_similarity_size": len(similarity_targets),
        "gsea": gsea_result,
        "analysis_config": config,
        "statistical_method": "two-sided Mann-Whitney U; asymptotic tie-aware p; BH within target; effect=2U/(n_target*n_control)-1",
        "log2fc_formula": "log2((mean normalized linear expression + 0.1)/(control mean normalized linear expression + 0.1)); normalized counts per 10,000 cells",
    }


def run_reactome(config, root, eligible, target_cell_map, control_idx,
                 de_effects, de_pct, de_pct_control, gene_ids, gene_names, confidence_rows):
    gmt = Path(config["reactome_gmt"])
    require_file(gmt)
    pathways: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    for line in gmt.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            pid = fields[0]
            names[pid] = fields[1]
            pathways[pid] = set(fields[2:])
    if not pathways:
        raise RuntimeError("Reactome GMT contains no pathways")
    confidence = {row[0]: row[-2] for row in confidence_rows}
    gsea_rows = []
    qc_rows = []
    permutation_by_target = {}
    expected_version = str(config.get("gseapy_version", gseapy.__version__))
    if expected_version and gseapy.__version__ != expected_version:
        raise RuntimeError(f"GSEApy version mismatch: configured {expected_version}, found {gseapy.__version__}")
    min_pct = float(config["min_expression_fraction"])
    min_coverage = float(config.get("min_pathway_coverage", 0.1))
    min_size, max_size = int(config["min_gene_set_size"]), int(config["max_gene_set_size"])
    for k, target in enumerate(eligible, start=1):
        effects = de_effects[target]
        expressed = (de_pct[target] >= min_pct) | (de_pct_control[target] >= min_pct)
        candidate_idx = np.flatnonzero(expressed & np.isfinite(effects))
        by_symbol: dict[str, int] = {}
        for i in candidate_idx:
            symbol = str(gene_names[i]).strip()
            if not symbol or symbol.lower() == "nan":
                continue
            current = by_symbol.get(symbol)
            if current is None or abs(float(effects[i])) > abs(float(effects[current])):
                by_symbol[symbol] = int(i)
        ranked_idx = sorted(by_symbol.values(), key=lambda i: (-float(effects[i]), str(gene_names[i]), str(gene_ids[i])))
        if not ranked_idx:
            raise RuntimeError(f"no expressed unique gene symbols for Reactome ranking: {target}")
        rank_symbols = [str(gene_names[i]) for i in ranked_idx]
        rank_scores = np.asarray([float(effects[i]) for i in ranked_idx], dtype=np.float64)
        unique_scores = int(np.unique(rank_scores).size)
        tie_fraction = 1.0 - unique_scores / len(rank_scores)
        zero_fraction = float(np.mean(rank_scores == 0))
        low_expression_fraction = float(1.0 - expressed.mean())
        is_priority = confidence.get(target) == str(config["gsea_priority_confidence"])
        permutations = int(config["gsea_priority_permutations"] if is_priority else config["gsea_permutations"])
        permutation_by_target[target] = permutations
        qc_rows.append((target, len(gene_ids), len(candidate_idx), len(rank_symbols), unique_scores,
                        tie_fraction, zero_fraction, low_expression_fraction, permutations,
                        "SUPPORTED priority" if is_priority else "standard"))
        ranked = pd.DataFrame({"gene": rank_symbols, "score": rank_scores})
        result = gseapy.prerank(
            rnk=ranked, gene_sets=pathways, outdir=None,
            min_size=min_size, max_size=max_size,
            permutation_num=permutations, seed=int(config["random_seed"]) + k,
            threads=int(config["analysis_threads"]), no_plot=True, verbose=False,
        )
        frame = result.res2d.reset_index(drop=True)
        target_rows = []
        rank_universe = set(rank_symbols)
        for _, item in frame.iterrows():
            term = str(item.get("Term", ""))
            pid = term if term in pathways else next((key for key, name in names.items() if name == term), "")
            pname = names.get(pid, term)
            leading = str(item.get("Lead_genes", "") or "")
            leading_genes = [x for x in re.split(r"[;/, ]+", leading) if x]
            pvalue = finite_float(item.get("NOM p-val", math.nan))
            if not math.isfinite(pvalue):
                pvalue = 1.0
            pvalue = max(pvalue, 1.0 / (permutations + 1))
            pathway_gene_count = len(pathways.get(pid, set()))
            overlap = len(rank_universe & pathways.get(pid, set()))
            coverage = overlap / pathway_gene_count if pathway_gene_count else 0.0
            target_rows.append({
                "target_gene": target, "target_cells": len(target_cell_map[target]),
                "control_cells": len(control_idx), "pathway_id": pid, "pathway_name": pname,
                "ES": finite_float(item.get("ES", math.nan)), "NES": finite_float(item.get("NES", math.nan)),
                "p_value": pvalue, "fdr": math.nan, "fwer_p_value": finite_float(item.get("FWER p-val", math.nan)),
                "leading_edge": ";".join(leading_genes), "leading_edge_size": len(leading_genes),
                "ranked_unique_symbols": len(rank_symbols), "rank_unique_scores": unique_scores,
                "pathway_gene_count": pathway_gene_count, "genes_in_ranked_universe": overlap,
                "coverage_fraction": coverage, "permutations": permutations,
                "priority_class": "SUPPORTED" if is_priority else "standard",
                "main_result_eligible": len(leading_genes) >= 3 and coverage >= min_coverage,
            })
        if target_rows:
            q = bh_adjust(np.asarray([row["p_value"] for row in target_rows], dtype=np.float64))
            for row, fdr in zip(target_rows, q):
                row["fdr"] = float(fdr)
        gsea_rows.extend(target_rows)
        print(f"REACTOME\t{k}/{len(eligible)}\t{target}\tgenes={len(rank_symbols)}\tpermutations={permutations}\tpathways={len(target_rows)}", flush=True)
    qc_columns = ("target_gene", "all_gex_features", "expressed_genes_before_symbol_collapse",
                  "ranked_unique_symbols", "unique_rank_scores", "tie_fraction", "zero_effect_fraction",
                  "low_expression_fraction_excluded", "permutations", "priority_class")
    write_tsv(root / "tables/gsea_ranking_qc.tsv", qc_columns, qc_rows)
    gsea_columns = ("target_gene", "target_cells", "control_cells", "pathway_id", "pathway_name", "ES", "NES",
                    "p_value", "fdr", "fwer_p_value", "leading_edge", "leading_edge_size",
                    "ranked_unique_symbols", "rank_unique_scores", "pathway_gene_count",
                    "genes_in_ranked_universe", "coverage_fraction", "permutations", "priority_class",
                    "main_result_eligible")
    write_tsv(root / "tables/pathway/reactome_GSEA_all.tsv.gz", gsea_columns,
              [tuple(row[col] for col in gsea_columns) for row in gsea_rows])
    if gsea_rows:
        gsea_df = pd.DataFrame(gsea_rows)
        fdr_values = pd.to_numeric(gsea_df["fdr"], errors="coerce")
        summary = {
            "database": config["reactome_database"], "version": config["reactome_version"],
            "gmt": str(gmt), "method": "GSEApy prerank with signed rank-biserial effect size",
            "gseapy_version": gseapy.__version__, "min_gene_set_size": min_size,
            "max_gene_set_size": max_size, "min_expression_fraction": min_pct,
            "min_main_pathway_coverage": min_coverage,
            "permutations_standard": int(config["gsea_permutations"]),
            "permutations_priority": int(config["gsea_priority_permutations"]),
            "priority_confidence_class": config["gsea_priority_confidence"],
            "random_seed_base": int(config["random_seed"]),
            "targets": len(eligible), "rows": len(gsea_rows),
            "targets_with_fdr_lt_0_05": int(gsea_df.loc[fdr_values < 0.05, "target_gene"].nunique()),
            "pathways_fdr_lt_0_05": int((fdr_values < 0.05).sum()),
            "pathways_main_eligible": int(gsea_df["main_result_eligible"].astype(bool).sum()),
            "p_value_floor": "1/(permutations+1); no empirical p-value is reported as zero",
            "fdr_scope": "Benjamini-Hochberg within each target across tested Reactome sets",
            "interpretation": "Exploratory cell-level pathway analysis from one biological sample.",
        }
    else:
        summary = {"database": config["reactome_database"], "version": config["reactome_version"],
                   "method": "GSEApy prerank with signed rank-biserial effect size", "targets": len(eligible), "rows": 0}
    write_json(root / "tables/pathway/pathway_metrics.json", summary)
    return summary


def write_figures(config, root, adata, summary):
    from matplotlib.colors import LinearSegmentedColormap

    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    manifest = []

    def save(fig, name, title, caption):
        path = figures / name
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        manifest.append((name, title, caption, str(path), "PASS"))

    obs = adata.obs.copy()
    total = pd.to_numeric(obs["total_counts"], errors="coerce").fillna(0).to_numpy()
    genes = pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").fillna(0).to_numpy()
    mito = pd.to_numeric(obs["pct_counts_mt"], errors="coerce").fillna(0).to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    for ax, values, label, color in zip(axes, (total, genes, mito),
                                        ("UMIs per cell", "Detected genes per cell", "Mitochondrial fraction (%)"),
                                        (NPG[3], NPG[2], NPG[0])):
        clean = values[np.isfinite(values) & (values >= 0)]
        ax.hist(clean, bins=60, color=color, alpha=0.88, edgecolor="white", linewidth=0.4)
        if label != "Mitochondrial fraction (%)":
            ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("Cells")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Gene-expression QC across DRAGEN passing cells", fontsize=14, fontweight="bold")
    save(fig, "fig02_gex_qc.png", "GEX 单细胞质量", "passing cells 的 nCount、nFeature 与线粒体比例分布；只展示已完成 cell calling 的细胞。")

    feature_path = root / "tables/cell_guide_feature_umi.tsv.gz"
    feature_df = pd.read_csv(feature_path, sep="\t")
    assignment_path = root / "tables/guide_assignment_all.tsv.gz"
    assign_df = pd.read_csv(assignment_path, sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    for ax, col, title, color in (
        (axes[0], "total_feature_umi", "Feature UMI per cell", NPG[1]),
        (axes[1], "top1_umi", "Top-guide UMI per cell", NPG[2]),
        (axes[2], "top1_fraction", "Top-guide fraction", NPG[0]),
    ):
        values = pd.to_numeric(feature_df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        ax.hist(values, bins=50, color=color, alpha=0.88, edgecolor="white", linewidth=0.4)
        ax.set_xlabel(title)
        ax.set_ylabel("Cells")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Guide-feature UMI support per filtered cell", fontsize=14, fontweight="bold")
    save(fig, "fig03_crispr_feature_qc.png", "CRISPR feature UMI 质量", "CRISPR read counts、分子数和单细胞 guide UMI 分布；read capture 不等同于可靠 singlet guide assignment。")

    coverage = pd.read_csv(root / "tables/target_coverage_thresholds.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    vals = coverage["targets_meeting_threshold"].astype(int).to_numpy()
    labels = [f"≥{int(x)} cells" for x in coverage["minimum_singlet_cells"]]
    bars = ax.bar(labels, vals, color=[NPG[i] for i in range(len(vals))], width=0.72)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylabel("Observed targets")
    ax.set_title("Target cell coverage thresholds")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig04_target_coverage.png", "靶点细胞覆盖", "横轴为 singlet cell 数门槛；分母为从 native guide assignment 观察到的全部靶点。")

    target_cov = pd.read_csv(root / "tables/target_guide_coverage.tsv", sep="\t")
    confidence = pd.read_csv(root / "tables/target_confidence.tsv", sep="\t")
    confidence_map = dict(zip(confidence["target_gene"], confidence["confidence_class"]))
    target_cov["class"] = target_cov["target_gene"].map(confidence_map).fillna("INCONCLUSIVE")
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    levels = list(dict.fromkeys(target_cov["class"].tolist()))
    colors = {value: NPG[i % len(NPG)] for i, value in enumerate(levels)}
    for level in levels:
        f = target_cov[target_cov["class"].eq(level)]
        ax.scatter(f["active_guides_ge5"], f["max_active_guide_fraction"],
                   s=20 + 10 * np.log1p(f["singlet_cells"].astype(float)), alpha=0.78,
                   color=colors[level], label=level, linewidths=0)
    ax.axhline(float(config["max_guide_fraction"]), color="#555555", linestyle="--", linewidth=0.9)
    ax.set_xlabel("Active guides (≥5 singlet cells)")
    ax.set_ylabel("Largest active-guide share of target singlets")
    ax.set_title("Guide balance and target support")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig05_guide_balance.png", "Guide 数量与平衡性", "点大小随靶点 singlet cells 增加；虚线为 GUIDE_DOMINATED 的 0.5 阈值。")

    consistency = pd.read_csv(root / "tables/guide_consistency_summary.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True)
    good = consistency[consistency["active_guides_ge5"] >= 2].copy()
    x = pd.to_numeric(good["median_pairwise_spearman"], errors="coerce")
    y = pd.to_numeric(good["median_directional_concordance"], errors="coerce")
    cc = [NPG[2] if s == "BALANCED" else "#B6BDC5" for s in good["guide_balance_status"]]
    ax.scatter(x, y, c=cc, s=24 + 7 * np.log1p(good["target_cells"].astype(float)), alpha=0.82, linewidths=0)
    ax.axhline(float(config["direction_supported_at_least"]), color="#666", linestyle="--", linewidth=0.8)
    ax.axhline(float(config["direction_discordant_below"]), color="#999", linestyle=":", linewidth=0.8)
    ax.axvline(0, color="#777", linewidth=0.7)
    ax.set_xlabel("Median pairwise Spearman rho")
    ax.set_ylabel("Median top-gene directional concordance")
    ax.set_title("Independent guide-level effect consistency")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig06_guide_consistency.png", "Guide-level 效应一致性", "每点一个有至少两个活跃 guide 的靶点；一致性指标基于 guide-vs-NTC 的全表达效应和 top 100 基因。")

    self_effect = pd.read_csv(root / "tables/target_self_expression_effect.tsv", sep="\t")
    x = pd.to_numeric(self_effect["pct_target"], errors="coerce")
    y = pd.to_numeric(self_effect["rank_biserial"], errors="coerce")
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    mask = x.notna() & y.notna()
    ax.scatter(x[mask], y[mask], s=24, alpha=0.74, color=NPG[0], linewidths=0)
    ax.axhline(0, color="#555", linewidth=0.8)
    ax.set_xlabel("Target-gene expression fraction in perturbation cells")
    ax.set_ylabel("Target-gene signed rank-biserial effect vs NTC")
    ax.set_title("Perturbation target self-expression")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig07_target_self_effect.png", "靶点自身表达效应", "CRISPR modality is UNKNOWN; target-gene RNA change is descriptive and is not used as a knockout success/failure test.")

    de_top = pd.read_csv(root / "tables/differential_expression/target_level_DE_top100.tsv", sep="\t")
    de_top["rank_biserial"] = pd.to_numeric(de_top["rank_biserial"], errors="coerce")
    de_top["delta_pct"] = pd.to_numeric(de_top["delta_pct"], errors="coerce")
    de_top["fdr"] = pd.to_numeric(de_top["fdr"], errors="coerce")
    de_top["support"] = pd.to_numeric(de_top["pct_target"], errors="coerce").fillna(0) + pd.to_numeric(de_top["pct_control"], errors="coerce").fillna(0)
    de_top["review_rank"] = de_top["rank_biserial"].abs() * de_top["delta_pct"].abs() * np.sqrt(de_top["support"])
    representative = de_top.sort_values(["review_rank", "fdr"], ascending=[False, True], kind="mergesort").head(24).copy()
    fig, ax = plt.subplots(figsize=(8.2, 7.0), constrained_layout=True)
    rep_labels = (representative["target_gene"].astype(str) + " · " + representative["gene_name"].astype(str)).tolist()
    vals = representative["rank_biserial"].to_numpy(dtype=float)
    sizes = 28 + 180 * np.clip(representative["support"].to_numpy(dtype=float) / 2, 0, 1)
    ax.scatter(vals, np.arange(len(vals)), c=vals, cmap=LinearSegmentedColormap.from_list("npg_signed", NPG_SIGNED),
               vmin=-1, vmax=1, s=sizes, edgecolors="white", linewidths=0.5)
    ax.set_yticks(np.arange(len(vals)), rep_labels, fontsize=7.2)
    ax.axvline(0, color="#555", linewidth=0.8)
    ax.set_xlabel("Signed rank-biserial effect (target vs NTC)")
    ax.set_title("Representative target-level expression effects")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "fig08_representative_target_effects.png", "代表性靶点效应", "按效应量、表达比例差和表达支持排序；圆点大小代表两组表达比例支持，不以极端 logFC 单独排序。")

    similarity = pd.read_csv(root / "tables/target_effect_similarity.tsv", sep="\t", index_col=0)
    order_json = read_json(root / "tables/target_effect_similarity_order.json")
    ordered = [x for x in order_json.get("cluster_order", []) if x in similarity.index]
    if ordered:
        similarity = similarity.loc[ordered, ordered]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), constrained_layout=True, gridspec_kw={"width_ratios": [1.0, 1.35]})
    if "X_umap" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap"])
        labels = np.asarray(adata.obs["target_genes"].astype(str))
        status = np.asarray(adata.obs["guide_status"].astype(str))
        counts = pd.Series(labels[(status == "singlet") & (labels != str(config["control_gene"]))]).value_counts()
        selected = set(counts.head(6).index.tolist())
        group = np.full(len(labels), "Other cells", dtype=object)
        group[labels == str(config["control_gene"])] = "NTC"
        for target in selected:
            group[(labels == target) & (status == "singlet")] = target
        levels = list(dict.fromkeys(group.tolist()))
        palette = {level: NPG[i % len(NPG)] for i, level in enumerate(levels)}
        for level in levels:
            m = group == level
            axes[0].scatter(xy[m, 0], xy[m, 1], s=2.0, alpha=0.72, color=palette[level], label=level, linewidths=0,
                            rasterized=True)
        axes[0].set_xlabel("UMAP1")
        axes[0].set_ylabel("UMAP2")
        axes[0].set_title("Existing UMAP · exploratory")
        axes[0].set_aspect("equal", adjustable="datalim")
        axes[0].legend(frameon=False, fontsize=7, markerscale=3, loc="best")
    else:
        axes[0].text(0.5, 0.5, "No UMAP coordinates in source H5AD", ha="center", va="center")
        axes[0].set_axis_off()
    im = axes[1].imshow(similarity.to_numpy(dtype=float), cmap=LinearSegmentedColormap.from_list("npg_signed", NPG_SIGNED), vmin=-1, vmax=1, aspect="auto")
    axes[1].set_xticks(range(len(similarity.columns)), similarity.columns, rotation=90, fontsize=5.5)
    axes[1].set_yticks(range(len(similarity.index)), similarity.index, fontsize=5.5)
    axes[1].set_title("Target-effect similarity · exploratory")
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.03, label="Spearman rho")
    fig.suptitle("Global cell-state view and target signature similarity", fontsize=13, fontweight="bold")
    save(fig, "fig09_umap_target_similarity.png", "UMAP 与靶点效应相似性", "UMAP 沿用原始全细胞嵌入，仅作描述；相似性为靶点效应向量的 Spearman 相关，不作因果解释。")

    pathway_path = root / "tables/pathway/reactome_GSEA_all.tsv.gz"
    pathway = pd.read_csv(pathway_path, sep="\t")
    if not pathway.empty:
        pathway["fdr"] = pd.to_numeric(pathway["fdr"], errors="coerce")
        pathway["NES"] = pd.to_numeric(pathway["NES"], errors="coerce")
        pathway["coverage_fraction"] = pd.to_numeric(pathway["coverage_fraction"], errors="coerce")
        pathway["leading_edge_size"] = pd.to_numeric(pathway["leading_edge_size"], errors="coerce")
        plot = pathway[pathway["main_result_eligible"].astype(str).str.lower().isin(["true", "1"])]
        plot = plot.sort_values(["fdr", "coverage_fraction", "leading_edge_size"], ascending=[True, False, False], kind="mergesort").head(20)
    else:
        plot = pathway
    fig, ax = plt.subplots(figsize=(9, 7.2), constrained_layout=True)
    if len(plot):
        lab = (plot["target_gene"].astype(str) + " · " + plot["pathway_name"].astype(str).str.slice(0, 58)).tolist()
        y = np.arange(len(plot))
        scatter = ax.scatter(plot["NES"], y,
                             c=plot["NES"], cmap=LinearSegmentedColormap.from_list("npg_signed", NPG_SIGNED),
                             vmin=-3, vmax=3, s=30 + 12 * plot["leading_edge_size"].fillna(0).to_numpy(),
                             edgecolors="white", linewidths=0.5)
        ax.set_yticks(y, lab, fontsize=7)
        fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.03, label="NES")
    else:
        ax.text(0.5, 0.5, "No pathways passed the main-result evidence filters", ha="center", va="center", transform=ax.transAxes)
        ax.set_yticks([])
    ax.axvline(0, color="#555", linewidth=0.8)
    ax.set_xlabel("Normalized enrichment score")
    ax.set_title("Reactome candidates · coverage and leading-edge filtered")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "fig10_reactome_pathways.png", "Reactome 通路候选", "主图仅显示 leading edge≥3 且基因集覆盖率达阈值的候选；单样本探索性排名富集。")

    conf = pd.read_csv(root / "tables/target_confidence.tsv", sep="\t")
    counts = conf["confidence_class"].value_counts().reindex(
        ["SUPPORTED", "PARTIAL_SUPPORT", "GUIDE_DOMINATED", "LOW_GUIDE_SUPPORT", "LOW_CELL", "DISCORDANT", "INCONCLUSIVE"],
        fill_value=0)
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    bars = ax.bar(counts.index, counts.values, color=[NPG[i] for i in range(len(counts))], width=0.72)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylabel("Observed targets")
    ax.set_title("Evidence classes from guide coverage, consistency, and LOOGO")
    ax.tick_params(axis="x", labelrotation=25)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "fig11_target_confidence.png", "靶点证据分类", "分类采用显式证据门槛；没有合并为单一分数。UNKNOWN 扰动模式下靶点自身表达不用于成功/失败判定。")

    write_tsv(root / "tables/figure_manifest.tsv",
              ("file", "title", "caption", "path", "status"), manifest)


def write_methods(root: Path, config: dict) -> None:
    content = f"""# 分析方法与判读边界

## 数据来源

本次分析读取样本 `{config['sample_id']}` 已完成的 DRAGEN 4.5.4-1 原生 PIPseq CRISPR 运行产物。GEX raw counts 从 filtered MEX 提取，按 feature 与 barcode 顺序同 filtered H5AD 的 `log_normalized` layer 对齐；与 per-cell `total_counts` 逐细胞核验。H5AD 的 `X` 已经过缩放，不作为 raw counts 使用。Guide assignment 使用 DRAGEN singlet assignment 和本项目 cell-level assignment 表；read 与 molecule 统计由 DRAGEN metrics 及 raw/filtered MEX 重算并交叉核对。该步骤是已完成原生 DRAGEN 输出的下游重分析，没有重新启动 DRAGEN。

## 细胞与 guide

单细胞 guide 分析限定为 DRAGEN passing cells 的 singlet calls。NTC 细胞作为表达比较基线。至少 {config['guide_min_singlet_cells']} 个 singlet cells 的 guide 纳入 guide-level 效应分析；靶点至少 {config['min_target_cells']} 个 singlet cells 且 NTC 至少 {config['min_control_cells']} 个细胞时纳入靶点级差异分析。阈值表同时保留低覆盖靶点，不将其静默删除。

## 差异表达与 guide 一致性

逐基因采用双侧 Mann–Whitney U 检验，使用渐近 tie-aware p 值，并在每个靶点的全部 GEX features 内分别做 Benjamini–Hochberg 校正。效应为 `2U/(n_target*n_control)-1` 的 signed rank-biserial effect；表达比例从 raw counts 的非零值计算；log2FC 从每细胞 10,000 总量归一化、log1p layer 反变换后的均值计算，伪计数为 0.1。guide-level 全基因表只提供效应量和表达支持，不报告单-guide p 值或 FDR。

guide 一致性基于 guide-vs-NTC 全基因效应、共同表达基因 Spearman 相关、top {config['top_genes_for_consistency']} effect 基因方向一致性和 leave-one-guide-out 稳健性。分类使用多项明示证据门槛；不是统计显著性或编辑效率判断。活跃 guide 最大份额阈值为 {config['max_guide_fraction']}；SUPPORTED 的 rho、方向一致性、LOOGO rho 和方向保留阈值分别为 {config['spearman_supported_at_least']}、{config['direction_supported_at_least']}、{config['loogo_spearman_min']} 和 {config['loogo_direction_retention_min']}。

## 通路分析

Reactome v{config['reactome_version']} GSEA 使用 GSEApy {config['gseapy_version']} prerank，以 signed rank-biserial effect 为排序分值。输入限于目标或 NTC 中至少 {config['min_expression_fraction']:.0%} 细胞表达的基因；同名符号保留绝对效应较大的 feature。基因集大小为 {config['min_gene_set_size']}–{config['max_gene_set_size']}，普通靶点 {config['gsea_permutations']} 次置换，SUPPORTED 靶点 {config['gsea_priority_permutations']} 次。经验 p 值下限为 `1/(置换数+1)`，不将零 p 值作为结果；FDR 在每个靶点内校正。主图候选要求 leading edge 至少 3 个基因且通路覆盖率至少 {config['min_pathway_coverage']:.0%}。

## 判读限制

该数据只有一个生物学样本。细胞不是生物学重复，细胞级 p 值与通路富集均属探索性结果，不能替代跨重复推断。扰动模式标为 UNKNOWN；靶点自身 RNA 变化仅描述表达证据，不能据此判定敲除成功或失败。既有 UMAP 来自输入 H5AD，仅用于描述性展示，本次未重新拟合降维或聚类。
"""
    (root / "methods/analysis_methods.md").write_text(content, encoding="utf-8")


def write_claim_evidence_map(root: Path, summary: dict) -> None:
    rows = [
        ("原生 DRAGEN CRISPR 已完成", "DRAGEN 运行日志与 stage_complete", "02_dragen_scrna_crispr/05_metadata/stage_complete.json; source dragen_run_*.log", "validated source provenance"),
        ("Read 分流及 read matching", "DRAGEN RNA/feature metrics", "tables/read_flow/read_flow_nodes.tsv; tables/read_flow/read_flow_links.tsv", "count conservation checks"),
        ("Passing-cell 和 guide assignment", "DRAGEN metrics + assignment tables", "tables/project_summary.json; tables/guide_assignment_all.tsv.gz", "counts partition to passing cells"),
        ("Molecule 数与 MEX 一致", "Raw/filtered MEX sums + DRAGEN metrics", "tables/project_summary.json; tables/read_flow/read_flow_nodes.tsv", "raw totals equal metrics"),
        ("Guide coverage and balance", "Singlet calls and feature counts", "tables/guide_coverage.tsv; tables/target_guide_coverage.tsv; figures/fig05_guide_balance.png", "one row per observed guide/target"),
        ("Target expression effects", "Filtered MEX raw gene counts + aligned H5AD log_normalized layer", "tables/differential_expression/target_level_DE_all.tsv.gz; figures/fig08_representative_target_effects.png", "cell-level exploratory; BH within target"),
        ("Guide consistency and robustness", "Per-guide rank effects, pairwise and leave-one-guide-out", "tables/guide_consistency_summary.tsv; tables/leave_one_guide_out_summary.tsv", "explicit evidence thresholds"),
        ("Reactome enrichment", "Signed rank-biserial effect and Reactome v97 GMT", "tables/pathway/reactome_GSEA_all.tsv.gz; tables/gsea_ranking_qc.tsv; figures/fig10_reactome_pathways.png", "single-sample exploratory; pathway coverage/leading-edge filters"),
        ("Target evidence classes", "Guide support, consistency and LOOGO", "tables/target_confidence.tsv; figures/fig11_target_confidence.png", "not an edit-efficiency or causal call"),
    ]
    write_tsv(root / "tables/claim_evidence_map.tsv",
              ("claim", "evidence_basis", "artifact", "interpretation_boundary"), rows)


def write_summary_docs(config: dict, root: Path, result: dict, report_validation: str) -> None:
    gsea = result.get("gsea", {})
    confidence = result.get("confidence_counts", {})
    report = root / "report/scPerturb_CRISPR_report_v2.html"
    lines = [
        "# scPerturb CRISPR 全量重分析摘要",
        "",
        f"样本 `{config['sample_id']}` 的原生 DRAGEN 4.5.4-1 PIPseq CRISPR 运行已在输入结果中完成；本次从已完成产物开始进行独立下游重分析，没有重新运行 DRAGEN。扰动模式为 `{config['perturbation_mode']}`。",
        "",
        "## 主要结果",
        "",
        f"- GEX reads：{int(result.get('gex_reads', 0)):,}；CRISPR feature reads：{int(result.get('crispr_feature_reads', 0)):,}。",
        f"- DRAGEN passing cells：{int(result.get('filtered_cells', 0)):,}；singlet / multi-guide / unassigned：{int(result.get('singlet_cells', 0)):,} / {int(result.get('multi_guide_cells', 0)):,} / {int(result.get('unassigned_cells', 0)):,}。",
        f"- Observed targets：{int(result.get('observed_targets', 0)):,}；target-level DE eligible：{int(result.get('eligible_targets', 0)):,}；NTC singlets：{int(result.get('ntc_cells', 0)):,}。",
        f"- Guide-level full-gene effects：{int(result.get('active_guides_de', 0)):,} guides；target-level DE：{int(result.get('target_de_rows', 0)):,} rows。",
        f"- Target evidence classes：{', '.join(f'{key}={value}' for key, value in sorted(confidence.items())) or 'not available'}。",
        f"- Reactome：{int(gsea.get('rows', 0)):,} pathway rows；FDR<0.05 pathways={int(gsea.get('pathways_fdr_lt_0_05', 0)):,}；满足主结果覆盖与 leading-edge 条件={int(gsea.get('pathways_main_eligible', 0)):,}。",
        "",
        "## 报告与验证",
        "",
        f"- HTML 报告：`{report}`。",
        f"- Report validation：**{report_validation}**。",
        "- Read/cell/molecule flow 的逐视图守恒和原生数据交叉检查见 `tables/read_flow/` 与 `05_metadata/sankey_validation.json`。",
        "- 结果表、图和方法说明均位于本版本目录；原始分析结果目录保持不变。",
        "",
        "## 科学解释范围",
        "",
        "单样本中细胞级差异表达和 GSEA 均为探索性结果，不能视为生物学重复层面的推断。扰动模式未知，靶点自身表达变化不用于判定编辑成功或失败。证据等级表示当前 guide 数量、效应一致性与 leave-one-guide-out 稳定性，不等价于编辑效率或因果结论。",
    ]
    (root / "REANALYSIS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    checks = [
        ("Native DRAGEN stage status", "PASS" if result.get("dragen_native_crispr") and result.get("dragen_license_success_in_source_log") and result.get("dragen_finished_normally") else "FAIL", "completed marker and source log checked"),
        ("Filtered H5AD / MEX dimensions and order", "PASS", "validated before analysis"),
        ("Molecule totals against DRAGEN metrics", "PASS", "raw and filtered MEX totals checked"),
        ("Guide assignments partition passing cells", "PASS", "singlet, multi and unassigned counts checked"),
        ("Read/cell/molecule Sankey", "PASS", "helper validator passed per-view conservation"),
        ("Analysis tables and figures", "PASS", "written by completed analysis steps"),
        ("Offline report validation", report_validation, "section order, offline resources, figures, tables and print CSS"),
        ("Completion marker", "written after this checklist", "created only after report validation passes"),
    ]
    checklist = ["# Validation checklist", "", "| Check | Status | Evidence |", "|---|---|---|"]
    checklist.extend(f"| {name} | {status} | {evidence} |" for name, status, evidence in checks)
    (root / "VALIDATION_CHECKLIST.md").write_text("\n".join(checklist) + "\n", encoding="utf-8")


def write_delivery_manifest(root: Path) -> Path:
    report = root / "report/scPerturb_CRISPR_report_v2.html"
    manifest = root / "delivery_manifest.tsv"
    files = [path for path in root.rglob("*") if path.is_file() and path != manifest
             and path.name != "stage_complete.json" and path.stat().st_size > 0]
    files = [path for path in files if path != report]
    rows = []
    for path in sorted(files, key=lambda item: item.relative_to(root).as_posix()):
        relative = path.relative_to(root).as_posix()
        category = relative.split("/", 1)[0]
        rows.append((category, relative, str(path.resolve()), "up_to_date"))
    if not report.is_file() or report.stat().st_size == 0:
        raise FileNotFoundError(f"final report missing or empty: {report}")
    rows.append(("report", "Final offline HTML report", str(report.resolve()), "up_to_date"))
    write_tsv(manifest, ("artifact", "description", "target", "status"), rows)
    return manifest


def write_sankey_delivery_manifest(source_manifest: Path, root: Path) -> Path:
    source_rows = read_tsv(source_manifest).to_dict("records")
    rows = []
    for item in source_rows:
        artifact = item["artifact"]
        if artifact == "read_flow_section":
            continue  # The interactive section is embedded in the final report.
        if artifact == "read_flow_nodes":
            target = "../tables/read_flow/read_flow_nodes.tsv"
        elif artifact == "read_flow_links":
            target = "../tables/read_flow/read_flow_links.tsv"
        elif artifact == "read_flow_transformations":
            target = "../tables/read_flow/read_flow_transformations.tsv"
        elif artifact == "read_flow_payload":
            target = "../tables/read_flow/read_flow_payload.json"
        elif artifact.startswith("read_flow_svg_"):
            target = "../figures/read_flow_sankey/" + Path(item["target"]).name
        elif artifact == "read_flow_validation":
            target = "sankey_validation.json"
        elif artifact == "read_flow_standalone":
            target = "../report/scPerturb_CRISPR_report_v2.html"
        else:
            raise RuntimeError(f"unrecognized read-flow manifest artifact: {artifact}")
        rows.append((artifact, item["description"], target, item["status"]))
    # The final report target is last, matching the same delivery convention as
    # the project-level manifest.
    rows.sort(key=lambda row: row[0] == "read_flow_standalone")
    path = root / "05_metadata/sankey_manifest.tsv"
    write_tsv(path, ("artifact", "description", "target", "status"), rows)
    return path


def validate_report(config: dict, root: Path) -> dict:
    report = root / "report/scPerturb_CRISPR_report_v2.html"
    manifest = root / "delivery_manifest.tsv"
    args = [sys.executable, "-m", "scperturb.validation.intermediate", "--report", str(report), "--manifest", str(manifest),
            "--require-section", "1 项目概览", "--require-expandable-figures", "--require-print-css",
            "--forbid-text", "100,000 reads GEX preflight", "--forbid-text", "100000 reads GEX preflight", "--json"]
    completed = subprocess.run(args, check=False, text=True, capture_output=True)
    try:
        result = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"bio-report validator returned invalid output: {completed.stdout}\n{completed.stderr}") from exc
    write_json(root / "05_metadata/report_validation.json", result)
    if completed.returncode != 0 or not result.get("ok"):
        raise RuntimeError(f"bio-report validator failed: {result}")
    return result


def main() -> int:
    args = cli()
    config_path = args.config.resolve()
    config = read_json(config_path)
    if args.validate_only:
        source = Path(config["source_root"])
        required = [
            source / "02_dragen_scrna_crispr/05_metadata/stage_complete.json",
            source / "02_dragen_scrna_crispr/04_tool_output/dragen" / f"{config['sample_id']}.scRNA_metrics.csv",
            source / "03_cell_guide_qc/04_tool_output" / f"{config['sample_id']}_gex_filtered.h5ad",
            Path(config["reactome_gmt"]),
        ]
        for path in required:
            require_file(path)
        if gseapy.__version__ != str(config["gseapy_version"]):
            raise RuntimeError(f"GSEApy version mismatch: expected {config['gseapy_version']}, found {gseapy.__version__}")
        print(json.dumps({"ok": True, "sample_id": config["sample_id"], "required_inputs": len(required),
                          "gseapy_version": gseapy.__version__}, ensure_ascii=False))
        return 0
    run(config, config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
