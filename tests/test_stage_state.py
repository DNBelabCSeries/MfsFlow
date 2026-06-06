import json
import os
import tempfile
import unittest
from types import SimpleNamespace

from mfsflow.path_layout import stage_state_dir
from mfsflow.stage_state import invalidate_stage_success, record_stage_success
from mfsflow.stages import MAPPING


class StageStateTests(unittest.TestCase):
    def make_runtime(self, outdir):
        return SimpleNamespace(
            project="sample",
            out_dir=outdir,
            tmp_merge_path=os.path.join(outdir, "intermediate", "tmp_merge"),
            config={"make_stats": True},
        )

    def test_mapping_success_writes_manifest_and_marker(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            bam = os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")

            manifest = record_stage_success(runtime, MAPPING)
            marker = os.path.join(stage_state_dir(tmpdir), "Mapping.success")
            self.assertTrue(os.path.exists(marker))
            payload = json.loads(manifest.read_text())
            self.assertEqual(payload["status"], "success")
            self.assertEqual(payload["artifacts"][0]["size_bytes"], 3)

            invalidate_stage_success(runtime, MAPPING)
            self.assertFalse(manifest.exists())
            self.assertFalse(os.path.exists(marker))

    def test_mapping_success_rejects_empty_required_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            open(os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam"), "wb").close()
            with self.assertRaisesRegex(RuntimeError, "missing or empty"):
                record_stage_success(runtime, MAPPING)


if __name__ == "__main__":
    unittest.main()
