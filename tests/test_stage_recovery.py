import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from mfsflow.stages.counting import cleanup_counting_inputs, run_counting_stage
from mfsflow.stages.mapping import cleanup_mapping_inputs
from mfsflow.stages.statistics import cleanup_statistics_inputs, run_statistics_stage


class CountingRecoveryTests(unittest.TestCase):
    def make_runtime(self, outdir):
        return SimpleNamespace(
            project="sample",
            analysis_dir=outdir,
            yaml_file=os.path.join(outdir, "run_config.yaml"),
            python_exec="python3",
            tools=SimpleNamespace(samtools="samtools"),
            resolve_script=lambda name: name,
            config={"make_stats": True},
        )

    def create_counting_inputs(self, outdir):
        names = [
            "sample.filtered.tagged.umi.Aligned.out.bam",
            "sample.filtered.tagged.internal.Aligned.out.bam",
            "sample.filtered.tagged.umi.Aligned.toTranscriptome.out.bam",
            "sample.filtered.tagged.internal.Aligned.toTranscriptome.out.bam",
        ]
        paths = []
        for name in names:
            path = os.path.join(outdir, name)
            open(path, "w").close()
            paths.append(path)
        return paths

    def test_dge_failure_keeps_mapping_outputs_for_resume(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = self.create_counting_inputs(tmpdir)

            def fail_dge(_cmd, stage_name):
                if stage_name == "dge_analysis.py":
                    raise RuntimeError("DGE failed")

            with self.assertRaisesRegex(RuntimeError, "DGE failed"):
                run_counting_stage(self.make_runtime(tmpdir), fail_dge)

            self.assertTrue(all(os.path.exists(path) for path in inputs))

    def test_counting_stage_keeps_mapping_outputs_until_manifest_is_written(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = self.create_counting_inputs(tmpdir)
            run_counting_stage(self.make_runtime(tmpdir), lambda _cmd, _stage_name: None)
            self.assertTrue(all(os.path.exists(path) for path in inputs))
            cleanup_counting_inputs(self.make_runtime(tmpdir))
            self.assertTrue(all(not os.path.exists(path) for path in inputs))

    def test_counting_rerun_removes_stale_optional_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = self.create_counting_inputs(tmpdir)
            stale_expression = os.path.join(tmpdir, "expression", "sample.intron.umi")
            os.makedirs(stale_expression)
            stale_h5ad = os.path.join(tmpdir, "expression", "sample.h5ad")
            open(stale_h5ad, "w").close()

            run_counting_stage(self.make_runtime(tmpdir), lambda _cmd, _stage_name: None)

            self.assertFalse(os.path.exists(stale_expression))
            self.assertFalse(os.path.exists(stale_h5ad))
            self.assertTrue(all(os.path.exists(path) for path in inputs))

    def test_disabled_statistics_rerun_removes_old_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.config["make_stats"] = False
            stats_dir = os.path.join(tmpdir, "stats")
            os.makedirs(stats_dir)
            stale = os.path.join(stats_dir, "sample.stats.tsv")
            open(stale, "w").close()
            run_cmd = mock.Mock()

            run_statistics_stage(runtime, run_cmd)

            self.assertFalse(os.path.exists(stale))
            run_cmd.assert_not_called()

    def test_mapping_cleanup_removes_only_completed_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            umi = os.path.join(tmpdir, "umi.bam")
            internal = os.path.join(tmpdir, "internal.bam")
            unrelated = os.path.join(tmpdir, "keep.bam")
            for path in (umi, internal, unrelated):
                open(path, "w").close()

            cleanup_mapping_inputs([umi], [internal])

            self.assertFalse(os.path.exists(umi))
            self.assertFalse(os.path.exists(internal))
            self.assertTrue(os.path.exists(unrelated))

    def test_statistics_cleanup_removes_gene_tagged_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            path = os.path.join(tmpdir, "sample.filtered.Aligned.GeneTagged.bam")
            open(path, "w").close()

            cleanup_statistics_inputs(runtime)

            self.assertFalse(os.path.exists(path))


if __name__ == "__main__":
    unittest.main()
