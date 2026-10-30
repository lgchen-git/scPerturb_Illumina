#!/usr/bin/env python3
"""Run DRAGEN PIPseq CRISPR mode with the GEX and sgRNA libraries together."""

from __future__ import annotations

import argparse
import csv
import gzip
import stat
import subprocess
import sys
from pathlib import Path

from scperturb.dragen.runtime import make_probe_fastq_list, run_dragen
from scperturb.common import (
    git_provenance,
    load_json,
    output_manifest,
    project_path,
    require_file,
    sample_row,
    set_thread_env,
    stage_root,
    write_completion,
    write_tsv,
)


def validate_mex(prefix: Path) -> tuple[int, int, int]:
    matrix = Path(str(prefix) + ".matrix.mtx.gz")
    barcodes = Path(str(prefix) + ".barcodes.tsv.gz")
    features = Path(str(prefix) + ".features.tsv.gz")
    for path, label in ((matrix, "matrix"), (barcodes, "barcodes"), (features, "features")):
        require_file(path, label)
    dimensions = None
    with gzip.open(matrix, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("%"):
                dimensions = tuple(int(value) for value in line.split())
                break
    if not dimensions or len(dimensions) != 3:
        raise ValueError(f"invalid Matrix Market dimensions: {matrix}")
    with gzip.open(barcodes, "rt", encoding="utf-8") as handle:
        barcode_count = sum(1 for _ in handle)
    with gzip.open(features, "rt", encoding="utf-8") as handle:
        feature_count = sum(1 for _ in handle)
    genes, cells, nnz = dimensions
    if genes != feature_count or cells != barcode_count:
        raise ValueError(f"MEX dimensions disagree: {dimensions} vs features={feature_count}, barcodes={barcode_count}")
    return genes, cells, nnz


def require_hdf5(path: Path, label: str) -> None:
    require_file(path, label)
    with path.open("rb") as handle:
        if handle.read(8) != b"\x89HDF\r\n\x1a\n":
            raise ValueError(f"invalid HDF5 signature: {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage-dir", required=True, type=Path)
    parser.add_argument("--license-check", choices=("strict", "product-validation"), default="product-validation")
    parser.add_argument("--license-probe-only", action="store_true")
    args = parser.parse_args()
    config = load_json(args.config)
    mode = config.get("quantification_mode", "dragen")
    if mode != "dragen":
        raise ValueError("Production workflow requires native DRAGEN quantification; fallback modes are not supported")
    row = sample_row(config, args.sample_id)
    stage = args.stage_dir.resolve()
    for name in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (stage / name).mkdir(parents=True, exist_ok=True)
    output_dir = stage / "04_tool_output" / "dragen"
    output_dir.mkdir(parents=True, exist_ok=True)
    feature_reference = stage_root(config, args.sample_id, "01_input_preflight") / "04_tool_output" / "dragen_feature_reference.csv"
    require_file(feature_reference, "converted CRISPR feature reference")
    dragen = require_file(config["dragen"]["binary"], "DRAGEN binary")
    star_binary = require_file(config["reference"]["star_binary"], "bundled STAR mapper")
    credential = require_file(config["dragen"]["credential"], "DRAGEN credential")
    if stat.S_IMODE(credential.stat().st_mode) != 0o600:
        raise PermissionError(f"DRAGEN credential must have mode 600: {credential}")
    ref_hash = require_file(Path(config["reference"]["dragen_hash_dir"]) / "hash_table.cfg", "DRAGEN hash table")
    require_file(Path(config["reference"]["dragen_hash_dir"]) / "reference.bin", "DRAGEN reference binary")
    gtf = require_file(config["reference"]["gtf"], "reference GTF")
    star_reference = require_file(Path(config["reference"]["star_reference"]) / "Genome", "STAR Genome index")
    fastq_list = stage / "05_metadata" / "fastq_list.csv"
    with fastq_list.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("RGID", "RGSM", "RGLB", "Lane", "Read1File", "Read2File"))
        for label, read1_key, read2_key in (
            ("GEX_1", "gex_r1", "gex_r2"),
            ("CRISPR_1", "crispr_r1", "crispr_r2"),
        ):
            writer.writerow((label, args.sample_id, row["rglb"], row["lane"],
                             str(project_path(config, row[read1_key])),
                             str(project_path(config, row[read2_key]))))
    intermediate = Path(config["tmp_root"]) / args.sample_id / "dragen_intermediate"
    intermediate.mkdir(parents=True, exist_ok=True)
    threads = int(config["resources"]["dragen_threads"])
    common_command = [
        str(dragen), "-f",
        "--lic-credentials", str(credential),
        "--ref-dir", str(Path(config["reference"]["dragen_hash_dir"])),
        "--annotation-file", str(gtf),
        "--output-directory", str(output_dir),
        "--output-file-prefix", args.sample_id,
        "--intermediate-results-dir", str(intermediate),
        "--scrna-enable-pipseq-mode", "true",
        "--scrna-enable-pipseq-crispr-mode", "true",
        "--scrna-feature-barcode-groups", "CRISPR_1",
        "--scrna-feature-barcode-reference", str(feature_reference),
        "--enable-sort", "true",
        "--rna-library-type", str(config["runtime"]["rna_library_type"]),
        "--scrna-filter-gene-biotypes", "true",
        "--single-cell-threshold", str(config["runtime"]["cell_threshold"]),
        "--enable-map-align-output", "false",
        "--single-cell-enable-h5ad-output", "true",
        "--single-cell-compress-intermediate-fastq", str(config["runtime"]["compress_intermediate_fastq"]).lower(),
        "--logging-to-output-dir", "true",
    ]
    # CRISPR mode must see both GEX and feature FASTQs. A GEX-only STAR BAM
    # cannot carry the feature reads or DRAGEN's CRISPR barcode transformation.
    command = common_command + [
        "--fastq-list", str(fastq_list),
        "--fastq-list-sample-id", args.sample_id,
        "--umi-source", "read1",
        "--single-cell-enable-star-mapper", "true",
        "--single-cell-star-mapper-path", str(star_binary),
        "--single-cell-star-mapper-reference-dir", str(Path(config["reference"]["star_reference"])),
        "--single-cell-star-mapper-options", f"--alignIntronMax 25000 --clip5pNbases 1 --runThreadN {threads}",
    ]
    command_record = write_tsv(stage / "05_metadata" / "dragen_command.tsv", ("position", "argument"), [(index, value if value != str(credential) else "<credential-path>") for index, value in enumerate(command)])
    license_result = subprocess.run(
        [sys.executable, "-m", "scperturb.dragen.license",
         "--config", str(args.config), "--license-check", args.license_check],
        check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    license_record = stage / "05_metadata" / "license_preflight.txt"
    license_record.write_text(license_result.stdout, encoding="utf-8")
    print(license_result.stdout, end="")
    if license_result.returncode:
        raise subprocess.CalledProcessError(license_result.returncode, license_result.args)
    version_result = run_dragen(config, [str(dragen), "--version"], stage=stage, intermediate=intermediate, capture=True)
    if version_result.returncode:
        raise subprocess.CalledProcessError(version_result.returncode, [str(dragen), "--version"], output=version_result.stdout)
    version = version_result.stdout or ""
    (stage / "05_metadata" / "dragen_version.txt").write_text(version, encoding="utf-8")
    env = set_thread_env(config, args.sample_id)
    env["RESOURCE_THREADS"] = str(threads)
    probe_root = Path(config["tmp_root"]) / args.sample_id / "license_probe"
    probe_output_dir = probe_root / "output"
    probe_intermediate = probe_root / "intermediate"
    probe_output_dir.mkdir(parents=True, exist_ok=True)
    probe_intermediate.mkdir(parents=True, exist_ok=True)
    probe_fastq_list = make_probe_fastq_list(
        config, row, probe_root, int(config["runtime"].get("preflight_reads", 100000)),
    )
    probe_command = list(command)
    replacements = {
        "--fastq-list": str(probe_fastq_list),
        "--output-directory": str(probe_output_dir),
        "--output-file-prefix": f"{args.sample_id}_license_probe",
        "--intermediate-results-dir": str(probe_intermediate),
    }
    for index, argument in enumerate(probe_command[:-1]):
        if argument in replacements:
            probe_command[index + 1] = replacements[argument]
    probe = run_dragen(config, probe_command, stage=stage, intermediate=probe_intermediate, capture=True)
    probe_output = probe.stdout or ""
    probe_record = stage / "05_metadata" / "dragen_license_probe.txt"
    probe_record.write_text(probe_output, encoding="utf-8")
    print(probe_output, end="")
    if probe.returncode or "error code 98" in probe_output.lower() or "error code 99" in probe_output.lower():
        raise RuntimeError(f"DRAGEN 100k-read license pilot failed; see {probe_record}")
    if args.license_probe_only:
        pilot_outputs = [probe_record, license_record, stage / "05_metadata" / "dragen_version.txt"]
        write_completion(stage, config, args.sample_id, pilot_outputs, interface="02-dragen-license-pilot-v1")
        print("PASS\tDRAGEN 100k-read GEX+CRISPR pilot completed; full FASTQ analysis was not started")
        return 0
    print("RUN\tDRAGEN 4.5.4-1 PIPseq CRISPR")
    run_result = run_dragen(config, command, stage=stage, intermediate=intermediate)
    if run_result.returncode:
        raise subprocess.CalledProcessError(run_result.returncode, command)

    filtered_prefix = output_dir / f"{args.sample_id}.scRNA.filtered"
    filtered_dimensions = validate_mex(filtered_prefix)
    if filtered_dimensions[1] <= 0 or filtered_dimensions[2] <= 0:
        raise ValueError("filtered matrix has no passing cells or nonzero entries")
    prefix = output_dir / f"{args.sample_id}.scRNA"
    raw_dimensions = validate_mex(prefix)
    for path, label in (
        (Path(str(prefix) + ".h5ad"), "raw h5ad"),
        (Path(str(filtered_prefix) + ".h5ad"), "filtered h5ad"),
        (Path(str(prefix) + ".moleculeInfo.h5"), "moleculeInfo h5"),
        (Path(str(prefix) + "_metrics.csv"), "scRNA metrics"),
        (Path(str(prefix) + ".barcodeSummary.tsv"), "barcode summary"),
        (Path(str(prefix) + ".positive_cell_guide_assignments.csv"), "positive guide assignments"),
        (Path(str(prefix) + ".singlet_cell_guide_assignments.csv"), "singlet guide assignments"),
        (Path(str(prefix) + ".guide_metrics.csv"), "guide metrics"),
    ):
        require_file(path, label)
    metrics = write_tsv(stage / "01_tables" / "dragen_matrix_qc.tsv", ("matrix", "features", "cells", "nonzero_entries"), [("raw", *raw_dimensions), ("filtered", *filtered_dimensions)])

    provenance = git_provenance(config["project_root"], stage)
    outputs = [fastq_list, command_record, stage / "05_metadata" / "dragen_version.txt", license_record, metrics, provenance]
    probe_record = stage / "05_metadata" / "dragen_license_probe.txt"
    if probe_record.is_file():
        outputs.append(probe_record)
    outputs.extend(path for path in output_dir.iterdir() if path.is_file())
    manifest = output_manifest(stage, outputs)
    outputs.append(manifest)
    write_completion(stage, config, args.sample_id, outputs, interface="02-dragen-pipseq-crispr-v4")
    print(f"PASS\t{args.sample_id}\tmode=dragen\traw_cells={raw_dimensions[1]}\tfiltered_cells={filtered_dimensions[1]}\tfiltered_features={filtered_dimensions[0]}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
