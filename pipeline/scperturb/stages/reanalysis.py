#!/usr/bin/env python3
"""Run the version-pinned primary scPerturb reanalysis on native DRAGEN outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from scperturb.common import load_analysis, load_json, sample_row, write_completion, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage-dir", type=Path, required=True)
    args = parser.parse_args()
    run = load_json(args.config)
    analysis = load_analysis(run)
    sample = sample_row(run, args.sample_id)
    source = Path(run["results_root"]) / args.sample_id
    stage_dir = args.stage_dir.resolve()
    temp_root = Path(run["tmp_root"]) / args.sample_id / "primary_reanalysis"
    temp_root.mkdir(parents=True, exist_ok=True)
    analysis_values = analysis.get("reanalysis", {})
    pathway = analysis["pathway"]
    project_id = sample["project_id"]
    source_markers = []
    for name in ("02_dragen_scrna_crispr", "03_cell_guide_qc"):
        marker_path = source / name / "05_metadata/stage_complete.json"
        if not marker_path.is_file():
            raise FileNotFoundError(f"upstream stage marker missing: {marker_path}")
        source_markers.append(json.loads(marker_path.read_text(encoding="utf-8")))
    seed = {
        "sample": sample, "analysis": analysis,
        "native_stage": source_markers[0], "cell_qc_stage": source_markers[1],
    }
    fingerprint = hashlib.sha256(json.dumps(seed, sort_keys=True, separators=(",", ":")).encode()).hexdigest()[:16]
    output_root = temp_root / fingerprint
    config = {
        "schema_version": 1, "sample_id": args.sample_id, "project_id": project_id,
        "source_root": str(source), "output_root": str(output_root),
        "pipeline_assets": str(Path(run["project_root"]) / "pipeline/assets"),
        "runtime_root": str(Path(run["tmp_root"]) / args.sample_id / "primary_runtime" / fingerprint),
        "reactome_gmt": str(pathway["gmt_file"]), "control_gene": analysis["guide_reference"]["control_gene"],
        "perturbation_mode": "UNKNOWN",
        "guide_min_singlet_cells": int(analysis_values.get("guide_min_singlet_cells", 5)),
        "min_target_cells": int(analysis["qc"]["min_target_cells"]),
        "min_control_cells": int(analysis["qc"]["min_control_cells"]),
        "min_expression_fraction": float(analysis_values.get("min_expression_fraction", 0.05)),
        "top_genes_for_consistency": int(analysis_values.get("top_genes_for_consistency", 100)),
        "max_guide_fraction": float(analysis_values.get("max_guide_fraction", 0.5)),
        "direction_discordant_below": float(analysis_values.get("direction_discordant_below", 0.5)),
        "direction_supported_at_least": float(analysis_values.get("direction_supported_at_least", 0.8)),
        "spearman_supported_at_least": float(analysis_values.get("spearman_supported_at_least", 0.3)),
        "loogo_spearman_min": float(analysis_values.get("loogo_spearman_min", 0.8)),
        "loogo_direction_retention_min": float(analysis_values.get("loogo_direction_retention_min", 0.8)),
        "gsea_permutations": int(pathway["permutation_num"]),
        "gsea_priority_permutations": int(analysis_values.get("gsea_priority_permutations", 5000)),
        "gsea_priority_confidence": str(analysis_values.get("gsea_priority_confidence", "SUPPORTED")),
        "gseapy_version": str(pathway["gseapy_version"]),
        "random_seed": int(pathway["random_seed"]),
        "analysis_threads": int(run["resources"]["analysis_threads"]),
        "gene_block_size": int(analysis_values.get("gene_block_size", 512)),
        "reactome_database": str(pathway["database"]), "reactome_version": str(pathway["version"]),
        "min_gene_set_size": int(pathway["min_gene_set_size"]),
        "max_gene_set_size": int(pathway["max_gene_set_size"]),
        "min_pathway_coverage": float(analysis_values.get("min_pathway_coverage", 0.1)),
        "ntc_cells": int(analysis_values.get("ntc_cells", analysis["qc"]["min_control_cells"])),
        "build_v2_report": True,
    }
    stage_dir.mkdir(parents=True, exist_ok=True)
    metadata = stage_dir / "05_metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    config_path = temp_root / f"reanalysis.{fingerprint}.json"
    write_json(config_path, config)
    subprocess.run([sys.executable, "-m", "scperturb.analysis.primary", "--config", str(config_path)], check=True,
                   cwd=run["project_root"])
    if not (output_root / "05_metadata/stage_complete.json").is_file():
        raise RuntimeError("primary reanalysis did not write its completion marker")
    write_json(metadata / "reanalysis_root.json", {
        "sample_id": args.sample_id, "fingerprint": fingerprint,
        "source_root": str(source), "v2_root": str(output_root),
    })
    marker = write_completion(stage_dir, run, args.sample_id,
                              [output_root / "05_metadata/stage_complete.json", metadata / "reanalysis_root.json"],
                              interface="04-primary-reanalysis-v3")
    print(json.dumps({"status": "complete", "sample_id": args.sample_id,
                      "fingerprint": fingerprint, "v2_root": str(output_root),
                      "marker": str(marker)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
