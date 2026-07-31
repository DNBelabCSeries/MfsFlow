import os
import tempfile
import unittest
from types import SimpleNamespace

from mfsflow.stages.counting import cleanup_counting_inputs, run_counting_stage
from mfsflow.stages.mapping import cleanup_mapping_inputs
from mfsflow.stages.statistics import cleanup_statistics_inputs


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
