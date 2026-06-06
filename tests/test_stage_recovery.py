import os
import tempfile
import unittest
from types import SimpleNamespace

from mfsflow.stages.counting import run_counting_stage


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

    def test_success_removes_mapping_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = self.create_counting_inputs(tmpdir)
            run_counting_stage(self.make_runtime(tmpdir), lambda _cmd, _stage_name: None)
            self.assertTrue(all(not os.path.exists(path) for path in inputs))


if __name__ == "__main__":
    unittest.main()
