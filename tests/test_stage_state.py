import json
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from mfsflow.path_layout import stage_state_dir
from mfsflow.stage_state import _quickcheck_bams, invalidate_stage_success, record_stage_success
from mfsflow.stages import FILTERING, MAPPING


class StageStateTests(unittest.TestCase):
    def make_runtime(self, outdir):
        return SimpleNamespace(
            project="sample",
            out_dir=outdir,
            tmp_merge_path=os.path.join(outdir, "intermediate", "tmp_merge"),
            config={"make_stats": True},
        )

    def test_filtering_quickcheck_allows_unmapped_header(self):
        runtime = SimpleNamespace(tools=SimpleNamespace(samtools="/tools/samtools"))
        completed = SimpleNamespace(returncode=0, stdout="")
        with mock.patch("mfsflow.stage_state.os.path.isfile", return_value=True), \
             mock.patch("mfsflow.stage_state.os.access", return_value=True), \
             mock.patch("mfsflow.stage_state.subprocess.run", return_value=completed) as run:
            _quickcheck_bams(runtime, ["/tmp/raw.tagged.bam"], unmapped=True)

        self.assertEqual(
            run.call_args.args[0],
            ["/tools/samtools", "quickcheck", "-v", "-u", "/tmp/raw.tagged.bam"],
        )

    @staticmethod
    def _write_no_sq_bam(path):
        """Write a valid BAM whose header has no @SQ lines (mimics fqfilter output)."""
        import pysam
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}}
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            a = pysam.AlignedSegment()
            a.query_name = "read1"
            a.flag = 4
            a.query_sequence = "ACGT"
            a.query_qualities = [30, 30, 30, 30]
            bam.write(a)

    def test_quickcheck_pysam_fallback_when_samtools_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = os.path.join(tmpdir, "raw.tagged.bam")
            self._write_no_sq_bam(bam)
            runtime = SimpleNamespace(tools=SimpleNamespace(samtools=None))
            _quickcheck_bams(runtime, [bam], unmapped=True)

    def test_quickcheck_pysam_fallback_when_samtools_old_rejects_u_flag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = os.path.join(tmpdir, "raw.tagged.bam")
            self._write_no_sq_bam(bam)
            runtime = SimpleNamespace(tools=SimpleNamespace(samtools="/tools/samtools"))
            completed = SimpleNamespace(returncode=2, stdout="unknown option -u")
            with mock.patch("mfsflow.stage_state.os.path.isfile", return_value=True), \
                 mock.patch("mfsflow.stage_state.os.access", return_value=True), \
                 mock.patch("mfsflow.stage_state.subprocess.run", return_value=completed):
                _quickcheck_bams(runtime, [bam], unmapped=True)

    def test_quickcheck_pysam_fallback_rejects_corrupt_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bam = os.path.join(tmpdir, "raw.tagged.bam")
            self._write_no_sq_bam(bam)
            with open(bam, "rb") as handle:
                data = handle.read()
            with open(bam, "wb") as handle:
                handle.write(data[: len(data) // 2])
            runtime = SimpleNamespace(tools=SimpleNamespace(samtools=None))
            with self.assertRaisesRegex(RuntimeError, "BAM integrity check failed"):
                _quickcheck_bams(runtime, [bam], unmapped=True)

    def test_filtering_manifest_uses_unmapped_quickcheck_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            os.makedirs(runtime.tmp_merge_path)
            bam = os.path.join(runtime.tmp_merge_path, "sample.part_001.raw.tagged.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")
            with mock.patch("mfsflow.stage_state._quickcheck_bams") as quickcheck:
                record_stage_success(runtime, FILTERING)

            self.assertTrue(quickcheck.call_args.kwargs["unmapped"])

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
