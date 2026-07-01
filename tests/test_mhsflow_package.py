import os
from pathlib import Path
import tempfile
import unittest

from mfsflow.cli import build_parser
from mfsflow.runtime import PipelineRuntime
from mfsflow.stages import COUNTING, FILTERING, MAPPING, STAGE_ORDER, SUMMARISING
from mfsflow.report import _select_report_template
from mfsflow.timer import format_duration


class MfsflowPackageTests(unittest.TestCase):
    def test_stage_order_is_stable(self):
        self.assertEqual(STAGE_ORDER, (FILTERING, MAPPING, COUNTING, SUMMARISING))

    def test_cli_parser_accepts_current_run_shape(self):
        args = build_parser().parse_args([
            "--fastqs", "raw",
            "--genomeDir", "ref",
            "--sample", "sample1",
            "--threads", "8",
            "--manual", "9",
        ])

        self.assertEqual(args.fastqs, "raw")
        self.assertEqual(args.genomeDir, "ref")
        self.assertEqual(args.sample, "sample1")
        self.assertEqual(args.threads, 8)
        self.assertEqual(args.manual, "9")

    def test_runtime_paths_are_derived_from_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config = {
                "project": "P1",
                "out_dir": os.path.join(tmpdir, "XPRESS_PROCESSING"),
                "num_threads": 4,
                "which_Stage": FILTERING,
                "toolkit_directory": tmpdir,
            }
            runtime = PipelineRuntime.from_config(config, "/tmp/run_config.yaml")

            self.assertEqual(runtime.project, "P1")
            self.assertEqual(runtime.num_threads, 4)
            self.assertTrue(runtime.log_path.endswith("pipeline.log"))
            self.assertTrue(runtime.timing_path.endswith("pipeline_timing.tsv"))

    def test_runtime_accepts_snake_case_stage_key(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config = {
                "project": "P1",
                "out_dir": os.path.join(tmpdir, "XPRESS_PROCESSING"),
                "num_threads": 4,
                "which_stage": MAPPING,
                "toolkit_directory": tmpdir,
            }
            runtime = PipelineRuntime.from_config(config, "/tmp/run_config.yaml")

            self.assertEqual(runtime.which_stage, MAPPING)

    def test_duration_formatting(self):
        self.assertEqual(format_duration(12.345), "12.35s")
        self.assertEqual(format_duration(65), "1m05.00s")

    def test_report_template_selection_uses_split_templates(self):
        template_dir = Path("/tmp/report")
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)

            self.assertEqual(
                _select_report_template("auto", outdir, template_dir)[0].name,
                "template_auto.html",
            )
            self.assertEqual(
                _select_report_template("manual", outdir, template_dir)[0].name,
                "template_manual.html",
            )
            self.assertEqual(
                _select_report_template("custom", outdir, template_dir)[0].name,
                "template_manual.html",
            )
            self.assertEqual(
                _select_report_template(
                    "discover",
                    outdir,
                    template_dir,
                    {"sample": {"discovered_sample_type": "manual"}},
                )[0].name,
                "template_manual.html",
            )


if __name__ == "__main__":
    unittest.main()
