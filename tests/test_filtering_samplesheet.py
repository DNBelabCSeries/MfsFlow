"""Tests for samplesheet-specific filtering flow."""

import contextlib
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

sys.modules.setdefault("yaml", mock.Mock())

from mfsflow.stages import filtering


class _Timer:
    def section(self, *_args, **_kwargs):
        return contextlib.nullcontext()


class _Proc:
    returncode = 0

    def wait(self):
        return None


class FilteringSamplesheetTests(unittest.TestCase):
    def test_samplesheet_discover_runs_discovery_after_bcstats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fastq = os.path.join(tmpdir, "r1.fq")
            with open(fastq, "w") as handle:
                handle.write("@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")

            analysis_dir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            tmp_merge = os.path.join(analysis_dir, "intermediate", "tmp_merge")
            os.makedirs(tmp_merge)
            os.makedirs(os.path.join(analysis_dir, "barcodes"))
            os.makedirs(os.path.join(analysis_dir, "config"))

            runtime = SimpleNamespace(
                config={
                    "barcode_source": "samplesheet_barcode",
                    "sample": {"sample_type": "discover"},
                    "sequence_files": {
                        "file1": {"name": fastq},
                        "file2": {"name": ""},
                    },
                    "counting_opts": {"max_reads": 0},
                },
                project="proj",
                out_dir=analysis_dir,
                num_threads=1,
                python_exec="python3",
                tools=SimpleNamespace(samtools="samtools", pigz="pigz", seqkit="seqkit"),
                exec_env={},
                analysis_dir=analysis_dir,
                tmp_merge_path=tmp_merge,
                yaml_file=os.path.join(analysis_dir, "config", "run_config.yaml"),
                resolve_script=lambda name: name,
                log_path=os.path.join(analysis_dir, "pipeline.log"),
            )

            run_stage_cmd = mock.Mock()
            with mock.patch("mfsflow.stages.filtering.pipeline_modules.split_fastq", return_value=[]), \
                 mock.patch("mfsflow.stages.filtering.pipeline_modules.merge_bam_stats"), \
                 mock.patch("mfsflow.stages.filtering.run_barcode_discovery") as discover, \
                 open(os.devnull, "w") as run_log:
                umi_chunks, int_chunks = filtering.run_filtering_stage(runtime, _Timer(), run_stage_cmd, run_log)

            discover.assert_called_once_with(runtime.config, "proj", analysis_dir)
            run_stage_cmd.assert_not_called()
            self.assertEqual(umi_chunks, [])
            self.assertEqual(int_chunks, [])

    def test_many_samplesheet_groups_skip_fastq_split(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            analysis_dir = os.path.join(tmpdir, "XPRESS_PROCESSING")
            tmp_merge = os.path.join(analysis_dir, "intermediate", "tmp_merge")
            os.makedirs(tmp_merge)
            os.makedirs(os.path.join(analysis_dir, "barcodes"))
            os.makedirs(os.path.join(analysis_dir, "config"))

            groups = []
            r1_files = []
            r2_files = []
            for idx in range(6):
                r1 = os.path.join(tmpdir, f"g{idx}_R1.fq")
                r2 = os.path.join(tmpdir, f"g{idx}_R2.fq")
                for path in (r1, r2):
                    with open(path, "w") as handle:
                        handle.write("@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
                r1_files.append(r1)
                r2_files.append(r2)
                groups.append({"read1": r1, "read2": r2, "barcode": f"BC{idx:018d}"})

            runtime = SimpleNamespace(
                config={
                    "barcode_source": "samplesheet_barcode",
                    "sample": {"sample_type": "manual"},
                    "sequence_files": {
                        "file1": {"name": ",".join(r1_files)},
                        "file2": {"name": ",".join(r2_files)},
                    },
                    "fastq_groups": groups,
                    "counting_opts": {"max_reads": 0},
                    "performance_opts": {"stream_bc_correction": True},
                },
                project="proj",
                out_dir=analysis_dir,
                num_threads=6,
                python_exec="python3",
                tools=SimpleNamespace(samtools="samtools", pigz="pigz", seqkit="seqkit"),
                exec_env={},
                analysis_dir=analysis_dir,
                tmp_merge_path=tmp_merge,
                yaml_file=os.path.join(analysis_dir, "config", "run_config.yaml"),
                resolve_script=lambda name: name,
                log_path=os.path.join(analysis_dir, "pipeline.log"),
            )

            launched = []

            def fake_popen(cmd, **_kwargs):
                launched.append(cmd)
                suffix = cmd[4]
                raw_bam = os.path.join(tmp_merge, f"proj{suffix}.raw.tagged.bam")
                open(raw_bam, "w").close()
                return _Proc()

            with mock.patch("mfsflow.stages.filtering.pipeline_modules.split_fastq") as split_fastq, \
                 mock.patch("mfsflow.stages.filtering.pipeline_modules.merge_bam_stats"), \
                 mock.patch("mfsflow.stages.filtering.subprocess.Popen", side_effect=fake_popen), \
                 open(os.devnull, "w") as run_log:
                umi_chunks, int_chunks = filtering.run_filtering_stage(runtime, _Timer(), mock.Mock(), run_log)

            split_fastq.assert_not_called()
            self.assertEqual(2, len(launched))
            self.assertTrue(all("--direct-fastq" in cmd for cmd in launched))
            self.assertEqual(umi_chunks, int_chunks)
            self.assertEqual(2, len(umi_chunks))


if __name__ == "__main__":
    unittest.main()
