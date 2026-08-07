import json
import os
import sys
import tempfile
import unittest
from unittest import mock

# Keep this orchestration test runnable in the lightweight local test
# environment; production and CI install PyYAML normally.
try:
    import yaml  # noqa: F401
except ImportError:
    sys.modules["yaml"] = mock.Mock()

from mfsflow.stages import COUNTING, FILTERING
from mfsflow.pipeline import runner


class PipelineRunnerTests(unittest.TestCase):
    def write_config(self, directory, stage=FILTERING):
        config_path = os.path.join(directory, "run_config.yaml")
        config = {
            "project": "sample",
            "out_dir": os.path.join(directory, "XPRESS_PROCESSING"),
            "num_threads": 2,
            "which_Stage": stage,
            "toolkit_directory": directory,
            "make_stats": True,
        }
        # JSON is valid YAML 1.2 and avoids coupling this orchestration test to
        # the optional test import shims used by barcode-focused tests.
        with open(config_path, "w", encoding="utf-8") as handle:
            json.dump(config, handle)
        return config_path, config

    def test_full_pipeline_records_each_stage_before_cleanup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path, config = self.write_config(tmpdir)
            events = []

            def record_success(_runtime, stage):
                events.append(f"record:{stage}")
                return stage

            def mapping_stage(*_args, **_kwargs):
                events.append("run:Mapping")
                return (["umi.bam"], ["internal.bam"])

            def cleanup_mapping(*_args):
                events.append("cleanup:Mapping")

            def counting_stage(*_args):
                events.append("run:Counting")

            def cleanup_counting(*_args):
                events.append("cleanup:Counting")

            def statistics_stage(*_args):
                events.append("run:Summarising")

            def cleanup_statistics(*_args):
                events.append("cleanup:Statistics")

            with mock.patch.object(runner.yaml, "safe_load", return_value=config), \
                 mock.patch.object(runner, "run_preflight"), \
                 mock.patch.object(runner, "run_filtering_stage", side_effect=lambda *_args: (events.append("run:Filtering") or (["raw.bam"], ["raw.bam"]))), \
                 mock.patch.object(runner, "run_mapping_stage", side_effect=mapping_stage), \
                 mock.patch.object(runner, "run_counting_stage", side_effect=counting_stage), \
                 mock.patch.object(runner, "run_statistics_stage", side_effect=statistics_stage), \
                 mock.patch.object(runner, "record_stage_success", side_effect=record_success), \
                 mock.patch.object(runner, "cleanup_mapping_inputs", side_effect=cleanup_mapping), \
                 mock.patch.object(runner, "cleanup_counting_inputs", side_effect=cleanup_counting), \
                 mock.patch.object(runner, "cleanup_statistics_inputs", side_effect=cleanup_statistics):
                runner.run_pipeline_stages(config_path)

            self.assertEqual(
                events,
                [
                    "run:Filtering",
                    "record:Filtering",
                    "run:Mapping",
                    "record:Mapping",
                    "cleanup:Mapping",
                    "run:Counting",
                    "record:Counting",
                    "run:Summarising",
                    "record:Summarising",
                    "cleanup:Counting",
                    "cleanup:Statistics",
                ],
            )

    def test_counting_failure_does_not_run_success_cleanup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path, config = self.write_config(tmpdir, stage=COUNTING)

            with mock.patch.object(runner.yaml, "safe_load", return_value=config), \
                 mock.patch.object(runner, "run_preflight"), \
                 mock.patch.object(runner, "run_counting_stage", side_effect=RuntimeError("DGE failed")), \
                 mock.patch.object(runner, "cleanup_counting_inputs") as cleanup_counting, \
                 mock.patch.object(runner, "cleanup_statistics_inputs") as cleanup_statistics:
                with self.assertRaisesRegex(RuntimeError, "DGE failed"):
                    runner.run_pipeline_stages(config_path)

            cleanup_counting.assert_not_called()
            cleanup_statistics.assert_not_called()


if __name__ == "__main__":
    unittest.main()
