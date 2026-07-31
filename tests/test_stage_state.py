import json
import os
import tempfile
import unittest
import importlib.util
from types import SimpleNamespace
from unittest import mock

from mfsflow.path_layout import stage_state_dir
from mfsflow.stage_state import (
    _quickcheck_bams,
    _validate_gzip,
    invalidate_stage_success,
    record_stage_success,
    validate_stage_manifest,
)
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
    @unittest.skipUnless(importlib.util.find_spec("pysam"), "pysam is not installed")
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

    def test_manifest_validation_rejects_changed_artifact(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            bam = os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")

            with mock.patch("mfsflow.stage_state._quickcheck_bams"):
                record_stage_success(runtime, MAPPING)
                with open(bam, "ab") as handle:
                    handle.write(b"changed")
                with self.assertRaisesRegex(RuntimeError, "size changed"):
                    validate_stage_manifest(runtime, MAPPING)

    def test_manifest_validation_missing_artifact_advises_rerun_from_filtering(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            bam = os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")

            with mock.patch("mfsflow.stage_state._quickcheck_bams"):
                record_stage_success(runtime, MAPPING)
                os.remove(bam)
                with self.assertRaisesRegex(RuntimeError, "missing or empty"):
                    validate_stage_manifest(runtime, MAPPING)

    def test_manifest_validation_allows_resume_stage_change(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.config["which_Stage"] = "Filtering"
            bam = os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")

            with mock.patch("mfsflow.stage_state._quickcheck_bams"):
                record_stage_success(runtime, MAPPING)
                runtime.config["which_Stage"] = "Mapping"
                validate_stage_manifest(runtime, MAPPING)

    def test_manifest_validation_accepts_stats_disabled_noop(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.config["make_stats"] = False
            with mock.patch("mfsflow.stage_state.stage_artifacts", return_value=[]):
                record_stage_success(runtime, "Summarising")
                payload = validate_stage_manifest(runtime, "Summarising")
            self.assertEqual(payload["artifacts"], [])

    def test_validate_gzip_accepts_valid_small_file(self):
        import gzip

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "valid.mtx.gz")
            with gzip.open(path, "wt") as handle:
                handle.write("1 2 3\n")
            _validate_gzip(path, "matrix")  # should not raise

    def test_validate_gzip_rejects_empty_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "empty.mtx.gz")
            open(path, "wb").close()
            with self.assertRaisesRegex(RuntimeError, "empty"):
                _validate_gzip(path, "matrix")

    def test_validate_gzip_rejects_non_gzip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "bad.mtx.gz")
            with open(path, "wb") as handle:
                handle.write(b"not gzip")
            with self.assertRaisesRegex(RuntimeError, "corrupt or unreadable"):
                _validate_gzip(path, "matrix")

    def test_validate_gzip_large_file_light_check_only_samples(self):
        import gzip

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "big.mtx.gz")
            with gzip.open(path, "wb") as handle:
                handle.write(b"x" * 4096)
            # Force the large-file branch by shrinking the full-check threshold,
            # then assert the file is still accepted.
            with mock.patch("mfsflow.stage_state._GZIP_FULL_CHECK_MAX_BYTES", 64):
                _validate_gzip(path, "matrix", full_check=False)  # should not raise

    def test_validate_gzip_full_check_rejects_truncated_large_file(self):
        import gzip
        import os

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "truncated.mtx.gz")
            with gzip.open(path, "wb") as handle:
                handle.write(os.urandom(4096))
            with open(path, "rb") as handle:
                data = handle.read()
            with open(path, "wb") as handle:
                handle.write(data[:-8])

            with mock.patch("mfsflow.stage_state._GZIP_FULL_CHECK_MAX_BYTES", 64):
                with self.assertRaisesRegex(RuntimeError, "corrupt or unreadable"):
                    _validate_gzip(path, "matrix", full_check=True)


if __name__ == "__main__":
    unittest.main()
