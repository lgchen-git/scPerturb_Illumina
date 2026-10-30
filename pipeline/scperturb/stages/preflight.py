#!/usr/bin/env python3
"""Validate GEX/CRISPR FASTQs and convert the supplied guide list."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path

from scperturb.common import (
    file_digest,
    git_provenance,
    load_analysis,
    load_json,
    output_manifest,
    require_file,
    sample_row,
    stage_root,
    write_completion,
    write_json,
    write_tsv,
)


EXPECTED_HOOK = "GTTTTAGAGCTAGAAATAGCA"
GUIDE_FIELDS = ("id", "name", "read", "target", "position", "sequence", "feature_type", "gene_id", "gene_name")
GEX_BARCODE_POSITIONS = ((0, 8), (11, 17), (20, 26), (31, 39))
CRISPR_BARCODE_POSITIONS = ((0, 8), (11, 17), (20, 26), (29, 37))


def barcode_blocks(path: Path) -> list[list[str]]:
    blocks: list[list[str]] = []
    current: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("#-Block"):
            if current:
                blocks.append(current)
            current = []
        elif line and not line.startswith("#"):
            current.append(line)
    if current:
        blocks.append(current)
    if len(blocks) != 4 or any(len(block) != 96 for block in blocks):
        raise ValueError(f"unexpected PIPseq barcode whitelist: {path}")
    return blocks


def near_sets(blocks: list[list[str]]) -> list[set[str]]:
    result = []
    for block in blocks:
        values = set(block)
        for sequence in block:
            for index, base in enumerate(sequence):
                for alternate in "ACGT":
                    if alternate != base:
                        values.add(sequence[:index] + alternate + sequence[index + 1 :])
        result.append(values)
    return result


def inspect_pair(read1: Path, read2: Path, limit: int, *, barcode_near: list[set[str]] | None = None, barcode_positions=GEX_BARCODE_POSITIONS, guide_by_seq: dict[str, str] | None = None) -> dict:
    records = 0
    r1_lengths: Counter[int] = Counter()
    r2_lengths: Counter[int] = Counter()
    phases: Counter[object] = Counter()
    guide_positions: Counter[int] = Counter()
    hooks: Counter[str] = Counter()
    guide_hits = 0
    unique_guides: set[str] = set()
    with gzip.open(read1, "rt", encoding="ascii", errors="strict") as left, gzip.open(read2, "rt", encoding="ascii", errors="strict") as right:
        while records < limit:
            h1 = left.readline()
            h2 = right.readline()
            if not h1 and not h2:
                break
            if not h1 or not h2:
                raise ValueError("paired FASTQs ended at different records")
            s1 = left.readline().strip().upper()
            s2 = right.readline().strip().upper()
            plus1 = left.readline()
            plus2 = right.readline()
            q1 = left.readline().strip()
            q2 = right.readline().strip()
            if h1.split()[0] != h2.split()[0]:
                raise ValueError(f"R1/R2 header mismatch: {h1.strip()} != {h2.strip()}")
            if not plus1.startswith("+") or not plus2.startswith("+") or len(s1) != len(q1) or len(s2) != len(q2):
                raise ValueError(f"malformed FASTQ record at sampled record {records + 1}")
            r1_lengths[len(s1)] += 1
            r2_lengths[len(s2)] += 1
            if barcode_near is not None:
                matches = []
                for phase in range(4):
                    if all(s1[start + phase:end + phase] in barcode_near[i] for i, (start, end) in enumerate(barcode_positions)):
                        matches.append(phase)
                phases[matches[0] if len(matches) == 1 else ("multi" if matches else "unmatched")] += 1
            if guide_by_seq is not None:
                found = []
                for position in range(10):
                    guide_id = guide_by_seq.get(s2[position:position + 20])
                    if guide_id is not None:
                        found.append((position, guide_id))
                if found:
                    guide_hits += 1
                    position, guide_id = found[0]
                    guide_positions[position] += 1
                    unique_guides.add(guide_id)
                    hooks[s2[position + 20:position + 20 + len(EXPECTED_HOOK)]] += 1
            records += 1
    if records == 0:
        raise ValueError(f"FASTQ pair contains no records: {read1}")
    matched = sum(phases[i] for i in range(4)) + phases["multi"]
    return {
        "records_tested": records,
        "r1_lengths": dict(r1_lengths),
        "r2_lengths": dict(r2_lengths),
        "phases": {str(key): value for key, value in phases.items()},
        "barcode_match_fraction": matched / records if barcode_near is not None else None,
        "guide_hits": guide_hits,
        "guide_hit_fraction": guide_hits / records if guide_by_seq is not None else None,
        "unique_guides": len(unique_guides),
        "guide_positions": dict(guide_positions),
        "hooks": dict(hooks),
    }


def load_guides(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"gene", "sgRNA_id", "sequence"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"guide source must contain columns gene,sgRNA_id,sequence: {path}")
    ids = [row["sgRNA_id"].strip() for row in rows]
    sequences = [row["sequence"].strip().upper() for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("guide source contains duplicate sgRNA_id")
    if len(sequences) != len(set(sequences)):
        raise ValueError("guide source contains duplicate guide sequences")
    if any(not re.fullmatch(r"[ACGT]{20}", sequence) for sequence in sequences):
        raise ValueError("guide source contains a non-20bp A/C/G/T sequence")
    normalized = []
    for row, sequence in zip(rows, sequences, strict=True):
        gene = row["gene"].strip()
        guide_id = row["sgRNA_id"].strip()
        if not gene or not guide_id:
            raise ValueError("guide source contains a blank gene or sgRNA_id")
        normalized.append({"gene": gene, "sgRNA_id": guide_id, "sequence": sequence})
    return normalized


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
    stage.mkdir(parents=True, exist_ok=True)
    for name in ("01_tables", "02_figures", "03_logs", "04_tool_output", "05_metadata"):
        (stage / name).mkdir(parents=True, exist_ok=True)

    guide_source = require_file(config["guide_source"], "guide source")
    guides = load_guides(guide_source)
    guide_by_seq = {guide["sequence"]: guide["sgRNA_id"] for guide in guides}
    feature_reference = stage / "04_tool_output" / "dragen_feature_reference.csv"
    with feature_reference.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GUIDE_FIELDS, lineterminator="\n")
        writer.writeheader()
        for guide in guides:
            writer.writerow({
                "id": guide["sgRNA_id"],
                "name": guide["sgRNA_id"],
                "read": analysis["guide_reference"]["read"],
                "target": analysis["guide_reference"]["target"],
                "position": analysis["guide_reference"]["position"],
                "sequence": guide["sequence"],
                "feature_type": analysis["guide_reference"]["feature_type"],
                "gene_id": guide["gene"],
                "gene_name": guide["gene"],
            })

    barcode_file = Path(config["dragen"]["install_root"]) / "dragen-softwaremode-4.5.4_1.el9.x86_64" / "resources" / "single_cell" / "pipseq_barcode_sequence.txt"
    require_file(barcode_file, "DRAGEN PIPseq barcode whitelist")
    barcode_near = near_sets(barcode_blocks(barcode_file))
    gex = inspect_pair(Path(row["gex_r1"]), Path(row["gex_r2"]), int(config["runtime"]["preflight_reads"]), barcode_near=barcode_near, barcode_positions=GEX_BARCODE_POSITIONS)
    crispr = inspect_pair(Path(row["crispr_r1"]), Path(row["crispr_r2"]), int(config["runtime"]["preflight_reads"]), barcode_near=barcode_near, barcode_positions=CRISPR_BARCODE_POSITIONS, guide_by_seq=guide_by_seq)
    if gex["barcode_match_fraction"] < 0.8:
        raise ValueError(f"GEX PIPseq barcode match fraction is too low: {gex['barcode_match_fraction']:.4f}")
    if crispr["guide_hits"] == 0:
        raise ValueError("no supplied guide sequence was found in sampled CRISPR R2 reads")
    if crispr["barcode_match_fraction"] < 0.8:
        raise ValueError(f"CRISPR PIPseq barcode match fraction is too low: {crispr['barcode_match_fraction']:.4f}")
    if not crispr["hooks"] or max(crispr["hooks"], key=crispr["hooks"].get) != EXPECTED_HOOK:
        raise ValueError(f"CRISPR hook does not match configured hook {EXPECTED_HOOK}")

    input_paths = [
        ("GEX_R1", Path(row["gex_r1"])), ("GEX_R2", Path(row["gex_r2"])),
        ("CRISPR_R1", Path(row["crispr_r1"])), ("CRISPR_R2", Path(row["crispr_r2"])),
    ]
    input_rows = [(args.sample_id, role, str(path.resolve()), path.stat().st_size, file_digest(path, "md5"), file_digest(path)) for role, path in input_paths]
    input_manifest = write_tsv(stage / "05_metadata" / "input_manifest.tsv", ("sample_id", "role", "path", "size_bytes", "md5", "sha256"), input_rows)
    guide_manifest = write_tsv(stage / "05_metadata" / "guide_source.tsv", ("source", "rows", "unique_genes", "sha256"), [(str(guide_source), len(guides), len({g['gene'] for g in guides}), file_digest(guide_source))])
    qc_rows = [
        ("guide_rows", len(guides)),
        ("guide_unique_ids", len({g["sgRNA_id"] for g in guides})),
        ("guide_unique_sequences", len({g["sequence"] for g in guides})),
        ("guide_unique_genes", len({g["gene"] for g in guides})),
        ("gex_records_tested", gex["records_tested"]),
        ("gex_r1_lengths", gex["r1_lengths"]),
        ("gex_r2_lengths", gex["r2_lengths"]),
        ("gex_barcode_match_fraction", f"{gex['barcode_match_fraction']:.6f}"),
        ("crispr_records_tested", crispr["records_tested"]),
        ("crispr_r1_lengths", crispr["r1_lengths"]),
        ("crispr_r2_lengths", crispr["r2_lengths"]),
        ("crispr_barcode_match_fraction", f"{crispr['barcode_match_fraction']:.6f}"),
        ("crispr_exact_guide_hit_fraction", f"{crispr['guide_hit_fraction']:.6f}"),
        ("crispr_unique_guides_seen", crispr["unique_guides"]),
        ("crispr_guide_positions", crispr["guide_positions"]),
        ("crispr_top_hook", max(crispr["hooks"], key=crispr["hooks"].get)),
        ("crispr_top_hook_count", max(crispr["hooks"].values())),
    ]
    qc = write_tsv(stage / "01_tables" / "input_qc.tsv", ("metric", "value"), qc_rows)
    qc_json = write_json(stage / "01_tables" / "input_qc.json", {
        "sample_id": args.sample_id,
        "guide_rows": len(guides),
        "gex_reads_checked": gex["records_tested"],
        "crispr_reads_checked": crispr["records_tested"],
        "gex_barcode_match_fraction": gex["barcode_match_fraction"],
        "crispr_barcode_match_fraction": crispr["barcode_match_fraction"],
        "crispr_guide_hit_fraction": crispr["guide_hit_fraction"],
        "crispr_unique_guides_seen": crispr["unique_guides"],
        "crispr_top_hook": max(crispr["hooks"], key=crispr["hooks"].get),
    })
    reference = write_tsv(stage / "05_metadata" / "reference_contract.tsv", ("field", "path"), [(key, value) for key, value in config["reference"].items()])
    provenance = git_provenance(config["project_root"], stage)
    manifest = output_manifest(stage, [feature_reference, input_manifest, guide_manifest, qc, qc_json, reference, provenance])
    outputs = [feature_reference, input_manifest, guide_manifest, qc, qc_json, reference, provenance, manifest]
    write_completion(stage, config, args.sample_id, outputs, interface="01-input-preflight-crispr-v2")
    print(f"PASS\t{args.sample_id}\tguides={len(guides)}\tgex_barcode_match={gex['barcode_match_fraction']:.4f}\tcrispr_guides_seen={crispr['unique_guides']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
