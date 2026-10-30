#!/usr/bin/env python3
"""Build cell-level GEX/guide QC tables from DRAGEN CRISPR output."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

from scperturb.common import (
    git_provenance,
    load_analysis,
    load_json,
    output_manifest,
    require_file,
    sample_row,
    set_thread_env,
    stage_root,
    write_completion,
    write_json,
    write_tsv,
)


def read_guide_source(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"gene", "sgRNA_id", "sequence"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"guide source must contain {sorted(required)}: {path}")
    return rows


def read_mex(prefix: Path):
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread

    matrix_path = Path(str(prefix) + ".matrix.mtx.gz")
    barcodes_path = Path(str(prefix) + ".barcodes.tsv.gz")
    features_path = Path(str(prefix) + ".features.tsv.gz")
    for path in (matrix_path, barcodes_path, features_path):
        require_file(path, f"filtered MEX component {path.name}")

    matrix = mmread(str(matrix_path)).tocsr()
    barcodes = []
    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            value = line.rstrip("\r\n").split("\t", 1)[0]
            if value:
                barcodes.append(value)
    feature_rows = []
    with gzip.open(features_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            feature_rows.append(fields)
    if matrix.shape[0] != len(feature_rows) or matrix.shape[1] != len(barcodes):
        raise ValueError(f"filtered MEX dimensions disagree: {matrix.shape}, features={len(feature_rows)}, barcodes={len(barcodes)}")
    feature_ids = [row[0] if row else "" for row in feature_rows]
    feature_names = [row[1] if len(row) > 1 else row[0] for row in feature_rows]
    feature_types = [row[2] if len(row) > 2 else "Gene Expression" for row in feature_rows]
    features = pd.DataFrame({"feature_id": feature_ids, "gene_name": feature_names, "feature_type": feature_types})
    return matrix, barcodes, features


def select_gene_features(features):
    feature_type = features["feature_type"].fillna("").astype(str).str.lower()
    mask = feature_type.str.contains("gene expression")
    if not mask.any():
        mask = ~feature_type.str.contains("crispr|guide|antibody|feature barcode")
    if not mask.any():
        raise ValueError("filtered MEX contains no identifiable gene-expression features")
    return mask.to_numpy()


def barcode_column(columns: list[str]) -> str:
    for column in columns:
        if "barcode" in column.lower():
            return column
    return columns[0]


def parse_assignment_file(path: Path, guide_ids: set[str]) -> dict[str, set[str]]:
    """Read the documented DRAGEN guide column or the diagnostic fallback column."""
    import pandas as pd

    if not path.is_file():
        return {}
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if frame.empty or frame.shape[1] == 0:
        return {}
    barcode_key = "cell_barcode" if "cell_barcode" in frame else "barcode"
    guide_key = next((key for key in ("feature_call", "guide_ids", "guide_id") if key in frame), None)
    if barcode_key not in frame or guide_key is None:
        raise ValueError(f"assignment CSV is missing barcode or guide-call column: {path}")
    assignments: dict[str, set[str]] = defaultdict(set)
    split_pattern = re.compile(r"[;,|\s]+")
    for _, record in frame.iterrows():
        barcode = str(record[barcode_key]).strip()
        if not barcode:
            continue
        for token in split_pattern.split(str(record[guide_key]).strip()):
            if token in guide_ids:
                assignments[barcode].add(token)
    return dict(assignments)


def make_figures(stage: Path, adata, status_counts: Counter[str], target_counts: Counter[str]) -> list[Path]:
    import matplotlib.pyplot as plt
    import numpy as np

    figure_dir = stage / "02_figures"
    figures: list[Path] = []
    qc_path = figure_dir / "cell_qc_distributions.png"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].hist(np.log10(adata.obs["total_counts"].to_numpy() + 1), bins=50, color="#2a6fbb", alpha=0.85)
    axes[0].set_xlabel("log10(total RNA counts + 1)")
    axes[0].set_ylabel("Cells")
    axes[1].hist(np.log10(adata.obs["n_genes_by_counts"].to_numpy() + 1), bins=50, color="#d95f02", alpha=0.85)
    axes[1].set_xlabel("log10(detected genes + 1)")
    axes[1].set_ylabel("Cells")
    fig.savefig(qc_path, dpi=160)
    plt.close(fig)
    figures.append(qc_path)

    status_path = figure_dir / "guide_assignment_status.png"
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    labels = list(status_counts)
    values = [status_counts[label] for label in labels]
    ax.bar(labels, values, color=["#8da0cb", "#1b9e77", "#d95f02"][: len(labels)])
    ax.set_ylabel("Filtered cells")
    ax.set_title("CRISPR guide assignment status")
    fig.savefig(status_path, dpi=160)
    plt.close(fig)
    figures.append(status_path)

    target_path = figure_dir / "top_target_cells.png"
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    items = target_counts.most_common(25)
    if items:
        labels, values = zip(*items)
        ax.barh(list(reversed(labels)), list(reversed(values)), color="#4c78a8")
    ax.set_xlabel("Singlet cells")
    ax.set_title("Top perturbation targets")
    fig.savefig(target_path, dpi=160)
    plt.close(fig)
    figures.append(target_path)

    return figures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage-dir", required=True, type=Path)
    args = parser.parse_args()
    config = load_json(args.config)
    analysis = load_analysis(config)
    row = sample_row(config, args.sample_id)
    stage = args.stage_dir.resolve()
    for name in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (stage / name).mkdir(parents=True, exist_ok=True)
    upstream = stage_root(config, args.sample_id, "02_dragen_scrna_crispr")
    prefix = upstream / "04_tool_output" / "dragen" / f"{args.sample_id}.scRNA.filtered"
    matrix, barcodes, features = read_mex(prefix)
    gene_mask = select_gene_features(features)

    env = set_thread_env(config, args.sample_id)
    # The current process inherited a safe cache policy from the Snakemake rule;
    # keep it explicit for scanpy/matplotlib/numba libraries loaded below.
    import os

    os.environ.update({key: value for key, value in env.items() if key in {"TMPDIR", "TMP", "TEMP", "XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"}})
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    gene_matrix = matrix[gene_mask, :].tocsr().T
    gene_features = features.loc[gene_mask].copy()
    gene_features.index = gene_features["feature_id"].astype(str)
    adata = ad.AnnData(X=gene_matrix, obs=pd.DataFrame(index=pd.Index(barcodes, name="barcode")), var=gene_features)
    adata.var_names_make_unique()
    adata.obs["total_counts"] = np.asarray(gene_matrix.sum(axis=1)).ravel()
    adata.obs["n_genes_by_counts"] = np.asarray((gene_matrix > 0).sum(axis=1)).ravel()
    mt_mask = adata.var["gene_name"].astype(str).str.upper().str.startswith("MT-").to_numpy()
    adata.obs["pct_counts_mt"] = 100 * np.asarray(gene_matrix[:, mt_mask].sum(axis=1)).ravel() / np.maximum(adata.obs["total_counts"].to_numpy(), 1)

    guide_rows = read_guide_source(config["guide_source"])
    guide_map = {record["sgRNA_id"]: record["gene"] for record in guide_rows}
    guide_ids = set(guide_map)
    assignment_dir = upstream / "04_tool_output" / "dragen"
    positive = parse_assignment_file(assignment_dir / f"{args.sample_id}.scRNA.positive_cell_guide_assignments.csv", guide_ids)
    singlet = parse_assignment_file(assignment_dir / f"{args.sample_id}.scRNA.singlet_cell_guide_assignments.csv", guide_ids)
    merged = {barcode: set(guides) for barcode, guides in positive.items()}
    for barcode, guides in singlet.items():
        merged.setdefault(barcode, set()).update(guides)
    if not merged:
        raise ValueError("DRAGEN guide assignment files contained no guide IDs matching the supplied feature reference")
    matrix_barcode_set = {str(barcode) for barcode in adata.obs_names}
    assignment_barcode_set = {str(barcode) for barcode in merged}
    linkage_overlap = matrix_barcode_set & assignment_barcode_set
    linkage_config = analysis.get("linkage", {})
    expected_random_overlap = len(matrix_barcode_set) * len(assignment_barcode_set) / float(96 ** 4)
    overlap_fraction = len(linkage_overlap) / max(1, len(assignment_barcode_set))
    enrichment = len(linkage_overlap) / max(expected_random_overlap, 1e-12)
    min_overlap_cells = int(linkage_config.get("min_overlap_cells", 10))
    min_overlap_fraction = float(linkage_config.get("min_overlap_fraction", 0.01))
    min_enrichment = float(linkage_config.get("min_enrichment", 10.0))
    quantitative_linked = (
        len(linkage_overlap) >= min_overlap_cells
        and overlap_fraction >= min_overlap_fraction
        and enrichment >= min_enrichment
    )
    mode = config.get("quantification_mode", "dragen")
    linkage_status = "NOT_VALIDATED" if mode == "fallback_diagnostic" else ("LINKED" if quantitative_linked else "LOW_LINKAGE")
    linkage_warning = ""
    if linkage_status == "NOT_VALIDATED":
        linkage_warning = (
            "Diagnostic fallback does not implement the documented CRISPR-to-GEX barcode transformation; "
            "observed overlap cannot establish whether the libraries are paired."
        )
    elif linkage_status != "LINKED":
        linkage_warning = (
            "Native DRAGEN barcodes have insufficient overlap between guide calls and filtered GEX cells; "
            "verify the library identity and DRAGEN feature metrics."
        )

    # Diagnostic fallback has no validated CRISPR-to-GEX barcode mapping.
    expression_assignments = merged if linkage_status == "LINKED" else {}

    assignment_rows = []
    target_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    for barcode in adata.obs_names:
        guides = sorted(expression_assignments.get(str(barcode), set()))
        genes = sorted({guide_map[guide] for guide in guides})
        if linkage_status == "NOT_VALIDATED":
            status = "not_validated"
        elif len(guides) == 0:
            status = "none"
        elif len(guides) == 1:
            status = "singlet"
        else:
            status = "multi"
        status_counts[status] += 1
        if status == "singlet":
            target_counts[genes[0]] += 1
        assignment_rows.append((str(barcode), status, len(guides), ";".join(guides), ";".join(genes), len(genes)))
    assignments_path = stage / "01_tables" / "cell_guide_assignments.tsv"
    write_tsv(assignments_path, ("barcode", "guide_status", "guide_count", "guide_ids", "target_genes", "target_gene_count"), assignment_rows)
    assignments = pd.DataFrame(assignment_rows, columns=["barcode", "guide_status", "guide_count", "guide_ids", "target_genes", "target_gene_count"]).set_index("barcode")
    adata.obs = adata.obs.join(assignments, how="left")
    adata.obs["guide_status"] = adata.obs["guide_status"].fillna("none").astype(str)
    adata.obs["target_genes"] = adata.obs["target_genes"].fillna("").astype(str)
    adata.obs["guide_ids"] = adata.obs["guide_ids"].fillna("").astype(str)

    try:
        sc.pp.normalize_total(adata, target_sum=float(analysis["expression"]["normalization_target_sum"]))
        sc.pp.log1p(adata)
        adata.layers["log_normalized"] = adata.X.copy()
        n_top = min(int(analysis["expression"]["n_top_genes"]), max(1, adata.n_vars - 1))
        detected_gene_mask = np.asarray(adata.X.getnnz(axis=0)).ravel() > 0
        hvg_adata = adata[:, detected_gene_mask].copy()
        hvg_top = min(n_top, max(1, hvg_adata.n_vars - 1))
        sc.pp.highly_variable_genes(hvg_adata, n_top_genes=hvg_top, flavor="cell_ranger")
        adata.var["highly_variable"] = False
        adata.var.loc[hvg_adata.var_names, "highly_variable"] = hvg_adata.var["highly_variable"].to_numpy()
        if int(adata.var["highly_variable"].sum()) < 2:
            adata.var["highly_variable"] = True
        sc.pp.scale(adata, max_value=10, zero_center=False)
        sc.tl.pca(adata, n_comps=min(int(analysis["expression"]["n_pcs"]), max(1, adata.n_obs - 1), max(1, adata.n_vars - 1)), svd_solver="arpack", random_state=int(analysis["expression"]["random_seed"]))
        sc.pp.neighbors(adata, n_neighbors=min(int(analysis["expression"]["n_neighbors"]), max(2, adata.n_obs - 1)), n_pcs=adata.obsm["X_pca"].shape[1])
        sc.tl.umap(adata, random_state=int(analysis["expression"]["random_seed"]))
        try:
            sc.tl.leiden(adata, resolution=float(analysis["expression"]["leiden_resolution"]), random_state=int(analysis["expression"]["random_seed"]), flavor="igraph", directed=False)
        except Exception:
            adata.obs["leiden"] = "0"
    except Exception as exc:
        raise RuntimeError(f"PCA/neighbors/UMAP failed: {exc}") from exc

    h5ad_path = stage / "04_tool_output" / f"{args.sample_id}_gex_filtered.h5ad"
    adata.write_h5ad(h5ad_path, compression="gzip")
    expression_qc = stage / "01_tables" / "expression_qc.tsv"
    write_tsv(expression_qc, ("metric", "value"), [
        ("filtered_cells", adata.n_obs),
        ("gene_features", adata.n_vars),
        ("median_total_counts", float(adata.obs["total_counts"].median())),
        ("median_detected_genes", float(adata.obs["n_genes_by_counts"].median())),
        ("median_pct_counts_mt", float(adata.obs["pct_counts_mt"].median())),
        ("singlet_cells", status_counts.get("singlet", 0)),
        ("multi_guide_cells", status_counts.get("multi", 0)),
        ("unassigned_cells", status_counts.get("none", 0)),
        ("assignment_cells_observed", len(assignment_barcode_set)),
        ("assignment_barcode_overlap", len(linkage_overlap)),
        ("assignment_barcode_overlap_fraction", overlap_fraction),
        ("expected_random_overlap", expected_random_overlap),
        ("observed_to_random_overlap_enrichment", enrichment),
        ("linkage_status", linkage_status),
    ])
    guide_summary_rows = []
    for record in guide_rows:
        guide = record["sgRNA_id"]
        positive_cells = sum(guide in guides for guides in positive.values())
        singlet_cells = sum(guide in guides and len(guides) == 1 for guides in merged.values())
        guide_summary_rows.append((guide, record["gene"], record["sequence"], positive_cells, singlet_cells))
    guide_summary = stage / "01_tables" / "guide_summary.tsv"
    write_tsv(guide_summary, ("sgRNA_id", "gene", "sequence", "positive_cells", "singlet_cells"), guide_summary_rows)
    qc_json = stage / "01_tables" / "cell_guide_qc.json"
    write_json(qc_json, {
        "sample_id": args.sample_id,
        "filtered_cells": int(adata.n_obs),
        "gene_features": int(adata.n_vars),
        "guide_source_records": len(guide_rows),
        "guide_source_genes": len({record["gene"] for record in guide_rows}),
        "assignment_status": dict(status_counts),
        "positive_assignment_cells": len(positive),
        "singlet_assignment_cells": len(singlet),
        "assignment_cells_observed": len(assignment_barcode_set),
        "assignment_barcode_overlap": len(linkage_overlap),
        "assignment_barcode_overlap_fraction": overlap_fraction,
        "expected_random_overlap": expected_random_overlap,
        "observed_to_random_overlap_enrichment": enrichment,
        "linkage_status": linkage_status,
        "linkage_thresholds": {
            "min_overlap_cells": min_overlap_cells,
            "min_overlap_fraction": min_overlap_fraction,
            "min_enrichment": min_enrichment,
        },
        "linkage_warning": linkage_warning,
        "top_targets": target_counts.most_common(40),
        "scanpy_version": sc.__version__,
    })
    figures = make_figures(stage, adata, status_counts, target_counts)
    provenance = git_provenance(config["project_root"], stage)
    outputs = [assignments_path, expression_qc, guide_summary, qc_json, h5ad_path, provenance, *figures]
    manifest = output_manifest(stage, outputs)
    outputs.append(manifest)
    write_completion(stage, config, args.sample_id, outputs, interface="03-cell-guide-qc-v3")
    print(f"PASS\t{args.sample_id}\tfiltered_cells={adata.n_obs}\tsinglet={status_counts.get('singlet', 0)}\tmulti={status_counts.get('multi', 0)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
