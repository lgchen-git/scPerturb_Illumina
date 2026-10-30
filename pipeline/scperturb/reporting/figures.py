#!/usr/bin/env python3
"""Create report-only UMAP and signed-effect figures from the existing V3 data."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CONTROL = "Non-Targeting Control"
from scperturb.reporting.style import NPG_CATEGORICAL, apply_matplotlib_style
NPG = list(NPG_CATEGORICAL)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="V3.2 delivery root")
    parser.add_argument("--source-root", required=True, type=Path, help="Native DRAGEN analysis root")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--control-gene", default="Non-Targeting Control")
    return parser.parse_args()


def natural_key(value: str):
    return (0, int(value)) if str(value).isdigit() else (1, str(value))


def save_figure(fig, path: Path):
    fig.savefig(path, dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def refresh_leiden_figures(delivery_root: Path) -> None:
    """Refresh descriptive Leiden plots from already-delivered tables only."""
    root = delivery_root.resolve()
    tables = root / "03_Tables"
    figures = root / "02_Figures/figures"
    composition = pd.read_csv(tables / "cell_state_composition.tsv", sep="\t", dtype=str, keep_default_na=False)
    summary = pd.read_csv(tables / "target_level_summary.tsv", sep="\t", dtype=str, keep_default_na=False)
    metadata = pd.read_csv(tables / "umap_cell_metadata.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    apply_matplotlib_style()

    summary["target_singlet_cells"] = pd.to_numeric(summary["target_singlet_cells"], errors="raise")
    display_groups = [CONTROL] + summary.sort_values(
        ["target_singlet_cells", "target_gene"], ascending=[False, True], kind="mergesort"
    )["target_gene"].astype(str).tolist()[:7]
    displayed = composition[composition["group"].isin(display_groups)].copy()
    clusters = sorted(set(composition["leiden_cluster"].astype(str)), key=natural_key)
    fig, ax = plt.subplots(figsize=(9.6, 4.8), constrained_layout=True)
    bottoms = np.zeros(len(display_groups))
    for index, cluster in enumerate(clusters):
        values = []
        for group in display_groups:
            row = displayed[(displayed["group"] == group) & (displayed["leiden_cluster"].astype(str) == cluster)]
            values.append(float(row.iloc[0]["fraction_within_group"]) if not row.empty else 0.0)
        ax.bar(display_groups, np.asarray(values) * 100, bottom=bottoms, color=NPG[index % len(NPG)],
               label=f"Leiden cluster {cluster}", width=.72)
        bottoms += np.asarray(values) * 100
    ax.set_ylabel("Singlet cells per group (%)")
    ax.tick_params(axis="x", labelrotation=30)
    ax.legend(frameon=False, fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, figures / "fig08_cell_state_composition.png")

    coords = metadata[["UMAP1", "UMAP2"]].apply(pd.to_numeric, errors="raise").to_numpy()
    cluster_labels = metadata["leiden"].astype(str).to_numpy()
    if coords.shape != (len(metadata), 2) or not np.isfinite(coords).all():
        raise ValueError("Delivered UMAP coordinates are malformed")
    levels = sorted(set(cluster_labels), key=natural_key)
    colors = dict(zip(levels, NPG))
    x_span = float(coords[:, 0].max() - coords[:, 0].min())
    y_span = float(coords[:, 1].max() - coords[:, 1].min())
    x_pad, y_pad = max(x_span * .035, .1), max(y_span * .035, .1)
    fig, ax = plt.subplots(figsize=(8.6, 6.5), constrained_layout=True)
    for level in levels:
        mask = cluster_labels == level
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2.1, c=colors[level], alpha=.84,
                   linewidths=0, rasterized=True, label=f"Leiden cluster {level}")
        center = np.median(coords[mask], axis=0)
        ax.text(center[0], center[1], level, ha="center", va="center", fontsize=8, fontweight="bold",
                color="#263238", bbox={"boxstyle": "round,pad=.22", "facecolor": "white",
                                      "edgecolor": "none", "alpha": .82})
    ax.set_xlim(float(coords[:, 0].min() - x_pad), float(coords[:, 0].max() + x_pad))
    ax.set_ylim(float(coords[:, 1].min() - y_pad), float(coords[:, 1].max() + y_pad))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    ax.set_title(f"Existing UMAP by Leiden cluster (n={len(metadata):,})")
    ax.legend(frameon=False, fontsize=8, title="Leiden cluster", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, figures / "fig09_global_cell_state_umap.png")

    manifest_path = tables / "figure_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    updates = {
        "fig08_cell_state_composition.png": (
            "Leiden-cluster composition of NTC and high-coverage targets",
            "Uses existing, unannotated Leiden cluster labels in the input H5AD. Bars are descriptive within-group fractions; no cell-level significance or mechanism is inferred.",
        ),
        "fig09_global_cell_state_umap.png": (
            "Existing UMAP with Leiden cluster labels",
            "Uses the existing embedding from the input H5AD; UMAP was not recomputed. Leiden clusters are unannotated grouping labels; UMAP proximity is descriptive and does not establish cell state or perturbation mechanism.",
        ),
    }
    for filename, (title, caption) in updates.items():
        mask = manifest["file"].eq(filename)
        manifest.loc[mask, "title"] = title
        manifest.loc[mask, "caption"] = caption
    manifest.to_csv(manifest_path, sep="\t", index=False)


def main() -> int:
    args = parse_args()
    global CONTROL
    root, source = args.root.resolve(), args.source_root.resolve()
    CONTROL = args.control_gene
    tables, figures = root / "tables", root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    apply_matplotlib_style()

    assignments = pd.read_csv(tables / "guide_assignment_all.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    target_summary = pd.read_csv(tables / "target_level_summary.tsv", sep="\t", dtype=str, keep_default_na=False)
    representatives = pd.read_csv(tables / "representative_target_effects.tsv", sep="\t", dtype=str, keep_default_na=False)
    h5ad_path = source / "03_cell_guide_qc/04_tool_output" / f"{args.sample_id}_gex_filtered.h5ad"
    adata = ad.read_h5ad(h5ad_path, backed="r")
    barcodes = list(map(str, adata.obs_names))
    if barcodes != assignments["barcode"].astype(str).tolist():
        adata.file.close()
        raise RuntimeError("H5AD barcode order differs from the V3 guide-assignment table")
    for column in ("guide_status", "target_genes"):
        if column not in adata.obs:
            adata.file.close()
            raise RuntimeError(f"Required H5AD obs column is absent: {column}")
        if not np.array_equal(adata.obs[column].astype(str).to_numpy(), assignments[column].astype(str).to_numpy()):
            adata.file.close()
            raise RuntimeError(f"H5AD {column} values/order differ from the V3 assignment table")
    if "leiden" not in adata.obs or "X_umap" not in adata.obsm:
        adata.file.close()
        raise RuntimeError("Input H5AD must contain the existing Leiden labels and X_umap")

    coords = np.asarray(adata.obsm["X_umap"])
    if coords.shape != (len(assignments), 2) or not np.isfinite(coords).all():
        adata.file.close()
        raise RuntimeError(f"Unexpected X_umap shape or non-finite coordinates: {coords.shape}")
    obs = adata.obs.copy()
    clusters = obs["leiden"].astype(str).to_numpy()
    statuses = assignments["guide_status"].astype(str).to_numpy()
    targets = assignments["target_genes"].astype(str).to_numpy()
    if set(np.unique(statuses)) - {"none", "singlet", "multi"}:
        adata.file.close()
        raise RuntimeError(f"Unrecognized guide assignment classes: {sorted(set(np.unique(statuses)))}")

    cell_meta = pd.DataFrame({
        "barcode": barcodes, "UMAP1": coords[:, 0], "UMAP2": coords[:, 1],
        "leiden": clusters, "guide_status": statuses, "target_gene": targets,
    })
    cell_meta.to_csv(tables / "umap_cell_metadata.tsv.gz", sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0})

    cluster_levels = sorted(set(clusters), key=natural_key)
    cluster_colors = dict(zip(cluster_levels, NPG))
    x_span = float(coords[:, 0].max() - coords[:, 0].min())
    y_span = float(coords[:, 1].max() - coords[:, 1].min())
    x_pad, y_pad = max(x_span * .035, .1), max(y_span * .035, .1)
    xlim = (float(coords[:, 0].min() - x_pad), float(coords[:, 0].max() + x_pad))
    ylim = (float(coords[:, 1].min() - y_pad), float(coords[:, 1].max() + y_pad))

    def finish_umap(ax):
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(False)

    fig, ax = plt.subplots(figsize=(8.6, 6.5), constrained_layout=True)
    for level in cluster_levels:
        mask = clusters == level
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2.1, c=cluster_colors[level], alpha=.84,
                   linewidths=0, rasterized=True, label=f"Leiden cluster {level}")
    for level in cluster_levels:
        mask = clusters == level
        center = np.median(coords[mask], axis=0)
        ax.text(center[0], center[1], str(level), ha="center", va="center", fontsize=8, fontweight="bold",
                color="#263238", bbox={"boxstyle": "round,pad=.22", "facecolor": "white", "edgecolor": "none", "alpha": .82})
    finish_umap(ax)
    ax.set_title(f"Global UMAP by Leiden cluster (n={len(cell_meta):,})")
    ax.legend(frameon=False, fontsize=8, title="Leiden cluster", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, figures / "fig09_global_cell_state_umap.png")

    assignment_specs = [("none", "Unassigned", "#B8BDC6"), ("singlet", "Singlet", "#00A087"),
                        ("multi", "Multi-guide", "#E64B35")]
    fig, ax = plt.subplots(figsize=(8.6, 6.5), constrained_layout=True)
    assignment_count_rows = []
    passing = int(len(assignments))
    for key, label, color in assignment_specs:
        mask = statuses == key
        count = int(mask.sum())
        assignment_count_rows.append({"view": "guide_assignment", "group": label, "cell_count": count,
                                      "denominator_cells": passing, "fraction_of_view": count / passing})
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2.0 if key == "none" else 3.2,
                   c=color, alpha=.52 if key == "none" else .86, linewidths=0, rasterized=True,
                   label=f"{label} (n={count:,})")
    finish_umap(ax); ax.set_title("Guide assignment across the existing embedding")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, figures / "fig10_guide_assignment_umap.png")

    singlet = statuses == "singlet"
    ntc = singlet & (targets == CONTROL)
    perturbation = singlet & (targets != "") & (targets != CONTROL)
    background = ~(ntc | perturbation)
    fig, ax = plt.subplots(figsize=(8.6, 6.5), constrained_layout=True)
    ax.scatter(coords[background, 0], coords[background, 1], s=1.6, c="#E5E7EB", alpha=.32,
               linewidths=0, rasterized=True, label=f"Other passing cells (n={int(background.sum()):,})")
    ax.scatter(coords[perturbation, 0], coords[perturbation, 1], s=3.4, c="#E64B35", alpha=.68,
               linewidths=0, rasterized=True, label=f"Perturbation singlets (n={int(perturbation.sum()):,})")
    ax.scatter(coords[ntc, 0], coords[ntc, 1], s=9, c="#4DBBD5", alpha=.95,
               linewidths=0, rasterized=True, label=f"NTC singlets (n={int(ntc.sum()):,})")
    finish_umap(ax); ax.set_title("NTC and perturbation singlets in the existing embedding")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, figures / "fig11_ntc_vs_perturbation_umap.png")
    assignment_count_rows.extend([
        {"view": "ntc_vs_perturbation", "group": "Other passing cells", "cell_count": int(background.sum()),
         "denominator_cells": passing, "fraction_of_view": float(background.sum() / passing)},
        {"view": "ntc_vs_perturbation", "group": "Perturbation singlets", "cell_count": int(perturbation.sum()),
         "denominator_cells": passing, "fraction_of_view": float(perturbation.sum() / passing)},
        {"view": "ntc_vs_perturbation", "group": "NTC singlets", "cell_count": int(ntc.sum()),
         "denominator_cells": passing, "fraction_of_view": float(ntc.sum() / passing)},
    ])

    targets_by_coverage = target_summary.copy()
    targets_by_coverage["target_singlet_cells"] = pd.to_numeric(targets_by_coverage["target_singlet_cells"], errors="raise")
    targets_by_coverage = targets_by_coverage.sort_values(["target_singlet_cells", "target_gene"],
                                                         ascending=[False, True], kind="mergesort")
    target_counts = pd.Series(targets[singlet & (targets != "") & (targets != CONTROL)]).value_counts()
    top_targets = targets_by_coverage.head(5)["target_gene"].astype(str).tolist()
    target_rows = []
    for rank, target in enumerate(top_targets, start=1):
        summary_row = targets_by_coverage.loc[targets_by_coverage["target_gene"].astype(str).eq(target)].iloc[0]
        source_n = int(summary_row["target_singlet_cells"])
        observed_n = int(target_counts.get(target, 0))
        if source_n != observed_n:
            adata.file.close()
            raise RuntimeError(f"Target cell-count mismatch for {target}: summary={source_n}; H5AD assignment={observed_n}")
        target_rows.append({"rank_by_coverage": rank, "target_gene": target, "target_singlet_cells": source_n,
                            "active_guides": int(summary_row["active_guides"]),
                            "largest_guide_share": float(summary_row["largest_guide_share"]),
                            "evidence_class": str(summary_row["evidence_class"]),
                            "coordinate_source": "existing X_umap"})
        mask = singlet & (targets == target)
        assignment_count_rows.append({"view": "representative_target_panels", "group": target, "cell_count": observed_n,
                                      "denominator_cells": passing, "fraction_of_view": observed_n / passing})

    fig = plt.figure(figsize=(12.6, 8.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 6)
    panel_slots = [(0, 0, 2), (0, 2, 4), (0, 4, 6), (1, 1, 3), (1, 3, 5)]
    axes = []
    for target, (row, col_start, col_stop) in zip(top_targets, panel_slots):
        ax = fig.add_subplot(grid[row, col_start:col_stop], sharex=axes[0] if axes else None,
                             sharey=axes[0] if axes else None)
        axes.append(ax)
        mask = singlet & (targets == target)
        ax.scatter(coords[:, 0], coords[:, 1], s=1.1, c="#E5E7EB", alpha=.24, linewidths=0,
                   rasterized=True, label="Other passing cells")
        ax.scatter(coords[mask, 0], coords[mask, 1], s=17, c="#00A087", alpha=.95, linewidths=0,
                   rasterized=True, label="Target singlets")
        n = int(mask.sum())
        ax.set_title(f"{target} (n={n:,})", fontsize=10)
        finish_umap(ax)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color="#E5E7EB", markersize=6, label="Other passing cells"),
               plt.Line2D([], [], marker="o", linestyle="", color="#00A087", markersize=6, label="Target singlets")]
    fig.legend(handles=handles, frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(.5, -.01))
    fig.suptitle("Representative target singlets on the same existing coordinates", y=1.01, fontsize=12)
    save_figure(fig, figures / "fig12_representative_target_umap.png")
    adata.file.close()

    pd.DataFrame(assignment_count_rows).to_csv(tables / "umap_group_counts.tsv", sep="\t", index=False)
    pd.DataFrame(target_rows).to_csv(tables / "umap_representative_target_counts.tsv", sep="\t", index=False)

    de_path = tables / "differential_expression/target_level_DE_all.tsv.gz"
    de_targets = representatives["target_gene"].astype(str).tolist()
    required_de = {"target_gene", "gene_id", "gene_name", "rank_biserial", "fdr", "pct_target", "pct_control", "delta_pct"}
    de_parts: dict[str, list[pd.DataFrame]] = {target: [] for target in de_targets}
    for chunk in pd.read_csv(de_path, sep="\t", dtype=str, keep_default_na=False, chunksize=200_000):
        missing = required_de - set(chunk.columns)
        if missing:
            raise RuntimeError(f"Primary DE source is missing required columns: {sorted(missing)}")
        part = chunk[chunk["target_gene"].astype(str).isin(de_parts)]
        for target, frame in part.groupby("target_gene", sort=False):
            de_parts[target].append(frame.copy())

    summary_by_target = targets_by_coverage.set_index("target_gene")
    de_records = []
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0), constrained_layout=True)
    for ax, target in zip(axes.flat, de_targets):
        pieces = de_parts[target]
        if not pieces:
            raise RuntimeError(f"No rows found in the primary DE table for representative target {target}")
        frame = pd.concat(pieces, ignore_index=True)
        effect = pd.to_numeric(frame["rank_biserial"], errors="coerce")
        fdr = pd.to_numeric(frame["fdr"], errors="coerce")
        pct_target = pd.to_numeric(frame["pct_target"], errors="coerce")
        pct_control = pd.to_numeric(frame["pct_control"], errors="coerce")
        delta = pd.to_numeric(frame["delta_pct"], errors="coerce")
        frame = frame.assign(_effect=effect, _fdr=fdr, _support=pct_target + pct_control, _delta=delta)
        supported = frame[(frame["_effect"].notna()) & (frame["_support"] >= .10) & (frame["_effect"].abs() > 0)].copy()
        if supported.empty:
            supported = frame[(frame["_effect"].notna()) & (frame["_support"] >= .10)].copy()
        selected = supported.sort_values(["_effect", "_support", "_delta", "_fdr", "gene_name"],
                                         key=lambda col: col.abs() if col.name == "_effect" or col.name == "_delta" else col,
                                         ascending=[False, False, False, True, True], kind="mergesort").head(6)
        if selected.empty:
            raise RuntimeError(f"No expression-supported genes for effect plot target {target}")
        selected = selected.sort_values("_effect", ascending=True, kind="mergesort")
        colors = ["#E64B35" if x > 0 else "#3C5488" for x in selected["_effect"]]
        ypos = np.arange(len(selected))
        ax.barh(ypos, selected["_effect"], color=colors, height=.68)
        ax.set_yticks(ypos, selected["gene_name"].astype(str), fontsize=8)
        ax.axvline(0, color="#555555", linewidth=.8)
        ax.set_xlim(-1, 1)
        counts_row = summary_by_target.loc[target]
        fdr_count = int(float(counts_row["de_genes_fdr_lt_0_05"]))
        ax.set_title(f"{target} (n={int(counts_row['target_singlet_cells'])}, guides={int(counts_row['active_guides'])}; FDR<0.05 genes={fdr_count})",
                     fontsize=9.5)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_xlabel("Signed rank-biserial effect")
        for record in selected.to_dict("records"):
            de_records.append({
                "target_gene": target, "target_singlet_cells": int(counts_row["target_singlet_cells"]),
                "active_guides": int(counts_row["active_guides"]), "evidence_class": str(counts_row["evidence_class"]),
                "gene_id": str(record["gene_id"]), "gene_name": str(record["gene_name"]),
                "rank_biserial": float(record["_effect"]), "fdr": float(record["_fdr"]),
                "pct_target": float(record["pct_target"]), "pct_control": float(record["pct_control"]),
                "delta_pct": float(record["_delta"]),
                "selection_rule": "combined expression fraction >=0.10; nonzero effect preferred; ranked by |rank_biserial|, expression support, |delta_pct|, then FDR",
            })
    for ax in list(axes.flat)[len(de_targets):]:
        ax.axis("off")
    save_figure(fig, figures / "fig13_representative_de_effects.png")
    pd.DataFrame(de_records).to_csv(tables / "representative_DE_plot_genes.tsv", sep="\t", index=False)

    old_umap = figures / "fig09_existing_umap_exploratory.png"
    if old_umap.exists():
        old_umap.unlink()
    manifest_path = tables / "figure_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    replaced_files = {"fig09_existing_umap_exploratory.png", "fig09_global_cell_state_umap.png",
                      "fig10_guide_assignment_umap.png", "fig11_ntc_vs_perturbation_umap.png",
                      "fig12_representative_target_umap.png", "fig13_representative_de_effects.png"}
    manifest = manifest[~manifest["file"].isin(replaced_files)].copy()
    if "report_revision" not in manifest:
        manifest["report_revision"] = "v3.2.1"
    else:
        manifest["report_revision"] = "v3.2.1"
    caption_umap = ("Uses the existing embedding from the input H5AD; UMAP was not recomputed. "
                    "Leiden clusters are unannotated grouping labels. UMAP proximity is descriptive and should not be interpreted as evidence of cell state or perturbation mechanism.")
    new_rows = [
        {"file": "fig09_global_cell_state_umap.png", "title": "Global UMAP by Leiden cluster", "section": "exploratory",
         "caption": caption_umap, "source_analysis_version": "v3", "report_revision": "v3.2.1"},
        {"file": "fig10_guide_assignment_umap.png", "title": "Guide assignment UMAP", "section": "exploratory",
         "caption": caption_umap, "source_analysis_version": "v3", "report_revision": "v3.2.1"},
        {"file": "fig11_ntc_vs_perturbation_umap.png", "title": "NTC vs perturbation UMAP", "section": "exploratory",
         "caption": caption_umap, "source_analysis_version": "v3", "report_revision": "v3.2.1"},
        {"file": "fig12_representative_target_umap.png", "title": "Representative target UMAP small multiples", "section": "exploratory",
         "caption": caption_umap, "source_analysis_version": "v3", "report_revision": "v3.2.1"},
        {"file": "fig13_representative_de_effects.png", "title": "Representative signed expression effects", "section": "effects",
         "caption": "Existing primary DE results; display genes were selected by expression support and signed effect size. Cell-level statistics are descriptive within one biological sample; this plot is not mechanism validation and does not alter evidence class.",
         "source_analysis_version": "v3", "report_revision": "v3.2.1"},
    ]
    manifest = pd.concat([manifest, pd.DataFrame(new_rows)], ignore_index=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)

    group_counts = pd.DataFrame([{"cluster": level, "cell_count": int(np.sum(clusters == level)),
                                  "denominator_cells": passing, "fraction_of_view": float(np.mean(clusters == level))}
                                 for level in cluster_levels])
    group_counts.insert(0, "view", "global_cell_state")
    group_counts.to_csv(tables / "umap_group_counts.tsv", sep="\t", index=False, mode="a", header=False)
    metadata = {
        "report_revision": "v3.2.1", "analysis_version": "v3", "embedding_key": "X_umap",
        "embedding_source_h5ad": str(h5ad_path), "embedding_recomputed": False,
        "same_coordinates_reused_in_all_umap_figures": True, "cells": int(len(coords)),
        "coordinate_columns": ["UMAP1", "UMAP2"], "coordinate_min": coords.min(axis=0).astype(float).tolist(),
        "coordinate_max": coords.max(axis=0).astype(float).tolist(), "leiden_column": "leiden",
        "leiden_levels": cluster_levels, "representative_umap_targets": target_rows,
        "representative_de_targets": de_targets,
        "de_display_rule": "Existing primary DE table; combined target/control expression fraction >=0.10; nonzero signed rank-biserial effects preferred; rank by absolute effect, expression support, absolute delta_pct, then FDR; no evidence class is modified.",
        "new_figures": [row["file"] for row in new_rows],
    }
    (root / "05_metadata/report_v3_2_1_figures.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "OK", "cells": len(coords), "clusters": len(cluster_levels),
                      "umap_targets": top_targets, "de_targets": de_targets,
                      "figure_count": len(manifest), "new_figures": len(new_rows)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
