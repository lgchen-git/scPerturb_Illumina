from __future__ import annotations

import json
import ast
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline"))
from scperturb.common import project_path

class WorkflowContractTests(unittest.TestCase):
    def test_templates_are_container_relative_and_native(self):
        run = json.loads((ROOT / "config/run.template.json").read_text())
        self.assertEqual(run["project_root"], "/project")
        self.assertEqual(run["dragen"]["credential"], "/run/secrets/dragen-credential")
        self.assertEqual(run["runtime"]["engine"], "docker")
        self.assertEqual(run["runtime"]["preflight_reads"], 100000)
        analysis = json.loads((ROOT / "config/perturbation.template.json").read_text())
        self.assertEqual(analysis["pathway"]["gseapy_version"], "1.2.1")
        self.assertEqual(analysis["guide_reference"]["control_gene"], "Non-Targeting Control")
        self.assertIn("illumina-pipseq-scrna:4.5.4-1-reports", (ROOT / "pipeline/Dockerfile").read_text())

    def test_report_stage_requires_official_html_and_metric_parity(self):
        snake = (ROOT / "pipeline/Snakefile").read_text()
        runner = (ROOT / "pipeline/scperturb/stages/report.py").read_text()
        validator = (ROOT / "pipeline/scperturb/validation/delivery.py").read_text()
        self.assertIn("official_report_path", snake)
        self.assertIn("parity_path", snake)
        self.assertIn("generate_official_report", runner)
        self.assertIn("dragen_report_parity.tsv", runner)
        self.assertIn("dragen_report_shared_metrics", validator)

    def test_fastq_paths_in_dragen_list_are_project_absolute(self):
        config = {"project_root": "/project"}
        self.assertEqual(project_path(config, "inputs/sample_R1.fq.gz"),
                         Path("/project/inputs/sample_R1.fq.gz"))
        self.assertEqual(project_path(config, "/project/inputs/sample_R2.fq.gz"),
                         Path("/project/inputs/sample_R2.fq.gz"))
        source = (ROOT / "pipeline/scperturb/stages/native_crispr.py").read_text()
        self.assertIn("str(project_path(config, row[read1_key]))", source)

    def test_metric_int_calls_include_sample_argument(self):
        source = (ROOT / "pipeline/scperturb/analysis/primary.py").read_text()
        tree = ast.parse(source)
        calls = [node for node in ast.walk(tree)
                 if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                 and node.func.id == "metric_int"]
        self.assertTrue(calls)
        self.assertTrue(all(len(node.args) == 4 for node in calls),
                        "metric_int calls must pass the sample selector")

    def test_report_builder_reads_read_flow_tables_from_subdirectory(self):
        source = (ROOT / "pipeline/scperturb/reporting/html.py").read_text()
        self.assertIn('read_flow_tables = tables / "read_flow"', source)
        self.assertIn('read_flow_tables / "read_flow_nodes.tsv"', source)
        self.assertIn('read_flow_tables / "read_flow_transformations.tsv"', source)

    def test_report_stage_flattens_preparation_metadata_into_validation(self):
        path = ROOT / "pipeline/scperturb/stages/report.py"
        spec = importlib.util.spec_from_file_location("report_stage", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "prepared/05_metadata"
            destination = root / "validation"
            source.mkdir(parents=True)
            destination.mkdir()
            (source / "v3_preparation.json").write_text("{}", encoding="utf-8")
            module.move_contents(source, destination)
            self.assertTrue((destination / "v3_preparation.json").is_file())
            self.assertFalse((destination / "05_metadata").exists())

    def test_active_pipeline_is_native_and_ai_independent(self):
        snake = (ROOT / "pipeline/Snakefile").read_text()
        dragen = (ROOT / "pipeline/scperturb/stages/native_crispr.py").read_text()
        self.assertIn("rule dragen_scrna_crispr", snake)
        self.assertIn("rule primary_reanalysis", snake)
        self.assertIn("--lic-credentials", dragen)
        self.assertIn("mode != \"dragen\"", dragen)
        for forbidden in ("openai", "anthropic", "langchain", "llm"):
            self.assertNotIn(forbidden, (snake + dragen).lower())

    def test_npg_style_palette_is_shared(self):
        source = (ROOT / "pipeline/scperturb/reporting/prepare.py").read_text()
        style = (ROOT / "pipeline/scperturb/reporting/style.py").read_text()
        self.assertIn("from scperturb.reporting.style import", source)
        self.assertIn("NPG_CATEGORICAL", style)
        self.assertIn("NPG_EVIDENCE", style)

    def test_analysis_code_has_no_sample_specific_identifiers(self):
        active = ["analysis/primary.py", "reporting/prepare.py",
                  "reporting/html.py", "reporting/figures.py"]
        for name in active:
            self.assertNotIn("CS2607604Per", (ROOT / "pipeline/scperturb" / name).read_text())

    def test_python_sources_compile_without_writing_bytecode(self):
        files = [ROOT / "pipeline/Snakefile"] + list((ROOT / "pipeline/scperturb").rglob("*.py"))
        for path in files:
            if path.name == "Snakefile":
                continue
            compile(path.read_text(encoding="utf-8"), str(path), "exec")

if __name__ == "__main__":
    unittest.main()
