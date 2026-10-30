"""Report-only isolation and preservation of the existing delivery on failure."""
import argparse
import csv
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scperturb import runner
from scperturb.common import delivery_name
from scperturb.icm import HANDOFF_DIRNAME, MANIFEST_NAME, RETURN_ASSETS, install_returned_assets
from scperturb.stages import report


class ReportLifecycleTests(unittest.TestCase):
    def fixture(self, root, samples=("sample_a", "sample_b")):
        (root / "config").mkdir()
        sheet = root / "config/samples.csv"
        sheet.write_text("sample_id\n" + "\n".join(samples) + "\n")
        cfg = {"project_root": "/project", "sample_sheet": "/project/config/samples.csv",
               "results_root": "/project/outputs", "tmp_root": "/project/tmp"}
        config = root / "config/run.json"
        config.write_text(json.dumps(cfg))
        for sample in samples:
            source = root / "outputs" / sample
            for stage in ("02_dragen_scrna_crispr", "03_cell_guide_qc", "04_primary_reanalysis", "05_report"):
                marker = source / stage / "05_metadata/stage_complete.json"
                marker.parent.mkdir(parents=True)
                marker.write_text(json.dumps({"sample_id": sample, "status": "complete"}))
            (source / "02_dragen_scrna_crispr/04_tool_output/dragen").mkdir(parents=True)
            returned = source / "02_dragen_scrna_crispr/04_tool_output" / HANDOFF_DIRNAME / "returned"
            returned.mkdir(parents=True)
            rows = []
            for role, asset in RETURN_ASSETS.items():
                content = b"\x89PNG\r\n\x1a\nfixture" if asset["filename"].endswith(".png") else b"%PDF-1.4\nfixture"
                (returned / asset["filename"]).write_bytes(content)
                rows.append((role, asset["filename"], asset["title"], asset["caption"]))
            with (returned / MANIFEST_NAME).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(("role", "filename", "title", "caption"))
                writer.writerows(rows)
            files = ["03_cell_guide_qc/01_tables/cell_guide_qc.json",
                     "05_report/04_tool_output/delivery/03_Tables/project_summary.json",
                     f"05_report/04_tool_output/delivery/01_Report/{sample}_scPerturb_CRISPR_report.html",
                     f"05_report/04_tool_output/delivery/05_DRAGEN_Report/{sample}_DRAGEN_report.html"]
            for relative in files:
                path = source / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("existing result")
        return argparse.Namespace(project=root, config=config)

    def test_report_ready_needs_no_fastq_reference_or_credential(self):
        with tempfile.TemporaryDirectory() as d:
            args = self.fixture(Path(d))
            _, _, samples = runner.reports_ready(args)
            self.assertEqual(len(samples), 2)
            self.assertFalse((Path(d) / "inputs").exists())
            self.assertFalse((Path(d) / "reference").exists())

    def test_configured_delivery_name_is_used_by_report_readiness(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = self.fixture(root)
            cfg = json.loads(args.config.read_text())
            cfg["delivery_name"] = "OUT_20261014_{sample_id}"
            args.config.write_text(json.dumps(cfg))
            self.assertEqual(delivery_name({}, "sample_a"), "delivery")
            for sample in ("sample_a", "sample_b"):
                parent = root / "outputs" / sample / "05_report/04_tool_output"
                (parent / "delivery").rename(parent / delivery_name(cfg, sample))
            _, _, samples = runner.reports_ready(args)
            self.assertEqual(len(samples), 2)

    def test_refresh_dispatches_only_report_for_each_sample(self):
        with tempfile.TemporaryDirectory() as d:
            args = self.fixture(Path(d))
            with patch("scperturb.runner.subprocess.run") as run, patch("scperturb.runner.verify"):
                runner.refresh(args)
            self.assertEqual(run.call_count, 2)
            for call in run.call_args_list:
                command = call.args[0]
                self.assertEqual(command[1:4], ["-m", "scperturb", "report"])
                self.assertIn("--refresh-report", command)

    def test_missing_upstream_stops_before_any_sample_refresh(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = self.fixture(root)
            (root / "outputs/sample_b/04_primary_reanalysis/05_metadata/stage_complete.json").unlink()
            with patch("scperturb.runner.subprocess.run") as run:
                with self.assertRaisesRegex(SystemExit, "requires completed stage"):
                    runner.refresh(args)
                run.assert_not_called()

    def test_missing_icm_exports_stops_before_any_sample_refresh(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = self.fixture(root)
            missing = (root / "outputs/sample_b/02_dragen_scrna_crispr/04_tool_output" /
                       HANDOFF_DIRNAME / "returned/analysis_summary.pdf")
            missing.unlink()
            with patch("scperturb.runner.subprocess.run") as run:
                with self.assertRaisesRegex(SystemExit, "requires completed ICM Perturb-seq exports"):
                    runner.refresh(args)
                run.assert_not_called()

    def test_icm_exports_are_staged_as_local_report_assets(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = self.fixture(root, samples=("sample_a",))
            handoff = (root / "outputs/sample_a/02_dragen_scrna_crispr/04_tool_output" /
                       HANDOFF_DIRNAME)
            delivery = root / "delivery"
            (delivery / "03_Tables").mkdir(parents=True)
            (delivery / "03_Tables/figure_manifest.tsv").write_text(
                "file\ttitle\tsection\tcaption\tsource_analysis_version\treport_revision\n",
                encoding="utf-8",
            )
            installed = install_returned_assets(handoff, delivery)
            self.assertEqual(set(installed), set(RETURN_ASSETS))
            self.assertTrue((delivery / "02_Figures/figures/icm_expression_overview.png").is_file())
            self.assertTrue((delivery / "02_Figures/figures/icm_guide_capture_overview.png").is_file())
            self.assertTrue((delivery / "04_Methods_Reproducibility/icm/analysis_summary.pdf").is_file())
            self.assertTrue((delivery / "04_Methods_Reproducibility/icm" / MANIFEST_NAME).is_file())

    def test_replaced_local_qc_and_umap_figures_leave_tables_untouched(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "prepared"
            figures, tables = root / "figures", root / "tables"
            figures.mkdir(parents=True)
            tables.mkdir()
            replaced = "fig02_gex_crispr_qc.png"
            retained = "fig13_representative_de_effects.png"
            (figures / replaced).write_bytes(b"figure")
            (figures / retained).write_bytes(b"figure")
            manifest = tables / "figure_manifest.tsv"
            manifest.write_text("file\ttitle\n" + replaced + "\tQC\n" + retained + "\tEffect\n", encoding="utf-8")
            report.remove_replaced_figures(root)
            self.assertFalse((figures / replaced).exists())
            self.assertTrue((figures / retained).is_file())
            self.assertEqual(manifest.read_text(encoding="utf-8"), "file\ttitle\n" + retained + "\tEffect\n")

    def test_validation_failure_preserves_old_delivery(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            stage, work = root / "stage", root / "work"
            final = stage / "04_tool_output/delivery"
            package = work / "delivery"
            final.mkdir(parents=True)
            (final / "keep.txt").write_text("original")
            (package / "04_Methods_Reproducibility").mkdir(parents=True)
            with patch("scperturb.stages.report.prepare_review", return_value={}), patch("scperturb.stages.report.finish_review", return_value={}), patch("scperturb.stages.report.validate", return_value={"result": "FAIL", "checks": [{"name": "official parity", "status": "FAIL"}]}):
                with self.assertRaisesRegex(RuntimeError, "existing delivery retained"):
                    report.publish_report({}, "sample", stage, package, work, status="refreshed")
            self.assertEqual((final / "keep.txt").read_text(), "original")
            self.assertTrue(package.is_dir())
            self.assertFalse((work / "previous_delivery").exists())

    def test_move_failure_restores_previous_delivery(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            stage, work = root / "stage", root / "work"
            final, package = stage / "04_tool_output/OUT_test_sample", work / "delivery"
            final.mkdir(parents=True)
            (final / "keep.txt").write_text("original")
            (package / "04_Methods_Reproducibility").mkdir(parents=True)
            move = shutil.move
            def fail_package(source, destination):
                if Path(source) == package:
                    raise OSError("simulated move failure")
                return move(source, destination)
            with patch("scperturb.stages.report.validate", return_value={"result": "PASS", "checks": []}), patch("scperturb.stages.report.prepare_review", return_value={}), patch("scperturb.stages.report.finish_review", return_value={"overall_result":"PASS"}), patch("scperturb.stages.report.shutil.move", side_effect=fail_package):
                with self.assertRaisesRegex(OSError, "simulated move failure"):
                    report.publish_report({"delivery_name": "OUT_test_sample"}, "sample", stage, package, work, status="refreshed")
            self.assertEqual((final / "keep.txt").read_text(), "original")


if __name__ == "__main__":
    unittest.main()
