#!/usr/bin/env python3
"""Host-side path validation and final artifact checks for the Docker runner."""
from __future__ import annotations
import argparse, csv, json, stat, subprocess, sys
from pathlib import Path
from scperturb.common import delivery_name
from scperturb.icm import handoff_root, prepare_handoff, validate_returned, validate_stage2

REQUIRED = ("config", "inputs", "log", "outputs", "pipeline", "reference", "tmp")

def inside_project(path: str, project: Path) -> Path:
    p = Path(path)
    if p.is_absolute() and str(p).startswith("/project/"):
        return project / str(p).removeprefix("/project/")
    if p.is_absolute():
        return p
    return project / p

def validate(args):
    root = args.project.resolve()
    missing = [name for name in REQUIRED if not (root / name).is_dir()]
    if missing: raise SystemExit("missing required project directories: " + ", ".join(missing))
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if cfg.get("project_root") != "/project": raise SystemExit("run config project_root must be /project inside Docker")
    for key, label in (("sample_sheet", "sample sheet"), ("guide_source", "guide reference"), ("analysis_config", "analysis config")):
        path = inside_project(cfg[key], root)
        if not path.is_file(): raise SystemExit(f"missing {label}: {path}")
    with inside_project(cfg["sample_sheet"], root).open(encoding="utf-8", newline="") as h:
        samples = list(csv.DictReader(h))
    if not samples or len({r["sample_id"] for r in samples}) != len(samples): raise SystemExit("sample sheet is empty or has duplicate sample IDs")
    for row in samples:
        for key in ("gex_r1", "gex_r2", "crispr_r1", "crispr_r2"):
            path = inside_project(row[key], root)
            if not path.is_file() or path.stat().st_size == 0: raise SystemExit(f"missing {key} FASTQ for {row['sample_id']}: {path}")
    reference = args.reference.resolve()
    for rel in ("reference/refdata-gex-GRCh38-2024-A-DRAGEN-v4.5/validated/dragen_hash/hash_table.cfg",
                "reference/refdata-gex-GRCh38-2024-A-DRAGEN-v4.5/validated/dragen_hash/reference.bin",
                "reference/refdata-gex-GRCh38-2024-A-DRAGEN-v4.5/validated/genes.gtf", "star_reference/Genome"):
        if not (reference / rel).is_file(): raise SystemExit(f"missing required reference file: {reference / rel}")
    pathway = inside_project(json.loads(inside_project(cfg["analysis_config"], root).read_text(encoding="utf-8"))["pathway"]["gmt_file"], root)
    if not pathway.is_file(): raise SystemExit(f"missing Reactome GMT: {pathway}")
    credential = args.credential.resolve()
    if not credential.is_file() or stat.S_IMODE(credential.stat().st_mode) != 0o600:
        raise SystemExit("DRAGEN credential must exist and have mode 600")
    if cfg.get("dragen", {}).get("credential") != "/run/secrets/dragen-credential":
        raise SystemExit("run config must use /run/secrets/dragen-credential")
    print(f"PRECHECK_PASS samples={len(samples)} layout=7-roots reference=present credential_mode=600")

def configured_samples(args):
    root = args.project.resolve()
    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    if cfg.get("project_root") != "/project":
        raise SystemExit("run config project_root must be /project")
    with inside_project(cfg["sample_sheet"], root).open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle))
    if not samples or len({row["sample_id"] for row in samples}) != len(samples):
        raise SystemExit("sample sheet is empty or has duplicate sample IDs")
    return root, cfg, samples


def reports_ready(args):
    root, cfg, samples = configured_samples(args)
    for row in samples:
        sample = row["sample_id"]
        source = inside_project(cfg["results_root"], root) / sample
        for stage in ("02_dragen_scrna_crispr", "03_cell_guide_qc", "04_primary_reanalysis", "05_report"):
            marker = source / stage / "05_metadata/stage_complete.json"
            if not marker.is_file():
                raise SystemExit(f"report-only requires completed stage: {marker}")
            doc = json.loads(marker.read_text(encoding="utf-8"))
            if doc.get("status") != "complete" or doc.get("sample_id") != sample:
                raise SystemExit(f"invalid completion marker: {marker}")
        delivery = source / "05_report/04_tool_output" / delivery_name(cfg, sample)
        required = [source / "02_dragen_scrna_crispr/04_tool_output/dragen",
                    source / "03_cell_guide_qc/01_tables/cell_guide_qc.json",
                    delivery / "03_Tables/project_summary.json",
                    delivery / f"01_Report/{sample}_scPerturb_CRISPR_report.html",
                    delivery / f"05_DRAGEN_Report/{sample}_DRAGEN_report.html"]
        for path in required:
            if not path.exists():
                raise SystemExit(f"report-only missing existing result: {path}")
        try:
            validate_returned(handoff_root(source / "02_dragen_scrna_crispr"))
        except (OSError, ValueError) as exc:
            raise SystemExit(f"report-only requires completed ICM Perturb-seq exports for {sample}: {exc}") from exc
    print(f"REPORT_READY samples={len(samples)} reuse=completed_analysis")
    return root, cfg, samples


def prepare_icm(args):
    root, cfg, samples = configured_samples(args)
    stage_roots = [(row["sample_id"], inside_project(cfg["results_root"], root) /
                    row["sample_id"] / "02_dragen_scrna_crispr") for row in samples]
    # Check every configured sample before writing any bundle.
    for sample, stage in stage_roots:
        validate_stage2(stage, sample)
    bundles = [prepare_handoff(stage, sample, dry_run=args.dry_run) for sample, stage in stage_roots]
    status = "ICM_HANDOFF_DRY_RUN" if args.dry_run else "ICM_HANDOFF_READY"
    print(f"{status} samples={len(bundles)} paths=" + ",".join(map(str, bundles)))


def refresh(args):
    root, cfg, samples = reports_ready(args)
    # Sequential samples keep memory bounded; no upstream stage is scheduled.
    for row in samples:
        stage = inside_project(cfg["results_root"], root) / row["sample_id"] / "05_report"
        subprocess.run([sys.executable, "-m", "scperturb", "report", "--config", str(args.config),
                        "--sample-id", row["sample_id"], "--stage-dir", str(stage), "--refresh-report"],
                       check=True, cwd=root)
    verify(args)


def review(args):
    root, cfg, samples = reports_ready(args)
    for row in samples:
        subprocess.run([sys.executable, "-m", "scperturb.validation.review", "--config", str(args.config),
                        "--sample-id", row["sample_id"]], check=True, cwd=root)
    verify(args)


def verify(args):
    from scperturb.dragen.alignment import SPECS
    root, cfg, samples = configured_samples(args)
    review_states = []
    for row in samples:
        sample = row["sample_id"]
        source = inside_project(cfg["results_root"], root) / sample
        for stage in ("01_input_preflight", "02_dragen_scrna_crispr", "03_cell_guide_qc", "04_primary_reanalysis", "05_report"):
            marker = source / stage / "05_metadata/stage_complete.json"
            if not marker.is_file(): raise SystemExit(f"missing completion marker: {marker}")
            payload = json.loads(marker.read_text(encoding="utf-8"))
            if payload.get("status") != "complete" or payload.get("sample_id") != sample:
                raise SystemExit(f"invalid completion marker: {marker}")
        base = source / "05_report/04_tool_output" / delivery_name(cfg, sample)
        for relative in (f"01_Report/{sample}_scPerturb_CRISPR_report.html",
                         f"05_DRAGEN_Report/{sample}_DRAGEN_report.html",
                         "03_Tables/dragen_official_metrics.tsv", "03_Tables/dragen_barcode_rank.json",
                         "02_Figures/figures/fig14_dragen_barcode_rank.png",
                         "04_Methods_Reproducibility/delivery_manifest.tsv"):
            path = base / relative
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"missing final artifact: {path}")
        validation = json.loads((base / "04_Methods_Reproducibility/validation/report_v3_2_1_validation.json").read_text())
        if validation.get("result") != "PASS" or not validation.get("checks") or any(c.get("status") != "PASS" for c in validation["checks"]):
            raise SystemExit(f"delivery validation not PASS: {sample}")
        with (base / "04_Methods_Reproducibility/validation/dragen_report_parity.tsv").open(newline="",encoding="utf-8") as handle:
            parity = list(csv.DictReader(handle, delimiter="\t"))
        if len(parity) != len(SPECS) or any(r.get("status") != "PASS" or r.get("sample_id") != sample for r in parity):
            raise SystemExit(f"official metric parity not PASS: {sample}")
        review_root = base / "04_Methods_Reproducibility/review"
        for suffix in ("md", "tsv", "json"):
            path = review_root / f"{sample}_DELIVERY_REVIEW.{suffix}"
            if not path.is_file() or not path.stat().st_size:
                raise SystemExit(f"missing delivery review: {path}; use --review-only")
        reviewed = json.loads((review_root / f"{sample}_DELIVERY_REVIEW.json").read_text())
        if reviewed.get("sample_id") != sample or reviewed.get("overall_result") not in {"PASS", "NOT_PASS", "PENDING"}:
            raise SystemExit(f"invalid delivery review: {sample}")
        review_states.append(f"{sample}:{reviewed['overall_result']}")
    print(f"COMPLETE_PASS samples={len(samples)} stages=5 official_metrics={len(SPECS)} validation=PASS delivery_review={','.join(review_states)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "reports-ready", "refresh", "review", "verify", "prepare-icm"):
        command = sub.add_parser(name)
        command.add_argument("--project", type=Path, required=True)
        command.add_argument("--config", type=Path, required=True)
        if name == "prepare-icm":
            command.add_argument("--dry-run", action="store_true")
        if name == "validate":
            command.add_argument("--reference", type=Path, required=True)
            command.add_argument("--credential", type=Path, required=True)
    args = parser.parse_args()
    {"validate": validate, "reports-ready": reports_ready, "refresh": refresh, "review": review,
     "verify": verify, "prepare-icm": prepare_icm}[args.command](args)


if __name__ == "__main__":
    main()
