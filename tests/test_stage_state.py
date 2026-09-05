import gzip
import json
import os
import tempfile
import unittest
import importlib.util
from types import SimpleNamespace
from unittest import mock

from mfsflow.path_layout import expression_dir, stage_state_dir
from mfsflow.stage_state import (
    _quickcheck_bams,
    _validate_gzip,
    invalidate_stage_success,
    record_stage_success,
    validate_resume_inputs,
    validate_stage_manifest,
)
from mfsflow.stages import COUNTING, FILTERING, MAPPING


class StageStateTests(unittest.TestCase):
    def make_runtime(self, outdir):
        return SimpleNamespace(
            project="sample",
            out_dir=outdir,
            tmp_merge_path=os.path.join(outdir, "intermediate", "tmp_merge"),
            config={
                "make_stats": True,
                "make_h5ad": False,
                "counting_opts": {"introns": False},
            },
            which_stage=COUNTING,
        )

    @staticmethod
    def write_mex_bundle(
        outdir,
        rows=2,
        columns=2,
        feature_rows=None,
        barcode_rows=None,
        entry="1 1 1",
        declared_entries=1,
    ):
        bundles = [
            os.path.join(expression_dir(outdir), "sample.exon.umi"),
            os.path.join(expression_dir(outdir), "sample.exon.read"),
        ]
        for bundle in bundles:
            os.makedirs(bundle)
        feature_rows = feature_rows if feature_rows is not None else rows
        barcode_rows = barcode_rows if barcode_rows is not None else columns
        for bundle in bundles:
            with gzip.open(os.path.join(bundle, "matrix.mtx.gz"), "wt") as handle:
                handle.write("%%MatrixMarket matrix coordinate integer general\n")
                handle.write("%\n")
                handle.write(f"{rows} {columns} {declared_entries}\n{entry}\n")
            with gzip.open(os.path.join(bundle, "features.tsv.gz"), "wt") as handle:
                for index in range(feature_rows):
                    handle.write(f"gene{index}\tGene {index}\tGene Expression\n")
            with gzip.open(os.path.join(bundle, "barcodes.tsv.gz"), "wt") as handle:
                for index in range(barcode_rows):
                    handle.write(f"BC{index}\n")
        return bundles[0]

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

    def test_mapped_quickcheck_uses_pysam_when_samtools_missing(self):
        runtime = SimpleNamespace(tools=SimpleNamespace(samtools=None))
        with mock.patch("mfsflow.stage_state._pysam_validate_bams") as fallback:
            _quickcheck_bams(runtime, ["/tmp/aligned.bam"])

        fallback.assert_called_once_with(["/tmp/aligned.bam"])

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

            with mock.patch("mfsflow.stage_state._quickcheck_bams"):
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

    def test_counting_success_validates_complete_mex_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            self.write_mex_bundle(tmpdir)

            manifest = record_stage_success(runtime, COUNTING)
            self.assertTrue(manifest.exists())
            validate_stage_manifest(runtime, COUNTING)

    def test_counting_rejects_mex_dimension_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            self.write_mex_bundle(tmpdir, rows=2, columns=2, feature_rows=1)

            with self.assertRaisesRegex(RuntimeError, "dimensions do not match"):
                record_stage_success(runtime, COUNTING)

    def test_counting_rejects_invalid_h5ad(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            self.write_mex_bundle(tmpdir)
            h5ad = os.path.join(expression_dir(tmpdir), "sample.h5ad")
            with open(h5ad, "wb") as handle:
                handle.write(b"not an hdf5 file")

            with self.assertRaisesRegex(RuntimeError, "not a valid HDF5 file"):
                record_stage_success(runtime, COUNTING)

    def test_counting_requires_intron_bundles_when_enabled(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.config["counting_opts"]["introns"] = True
            self.write_mex_bundle(tmpdir)

            with self.assertRaisesRegex(RuntimeError, "missing required expression matrix"):
                record_stage_success(runtime, COUNTING)

    def test_resume_validates_matrix_indices_and_nnz(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            self.write_mex_bundle(tmpdir, entry="3 1 1")

            record_stage_success(runtime, COUNTING)
            with self.assertRaisesRegex(RuntimeError, "index out of range"):
                validate_stage_manifest(runtime, COUNTING)

    def test_summarising_resume_validates_all_mex_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.which_stage = "Summarising"
            bundle = self.write_mex_bundle(tmpdir)
            with gzip.open(os.path.join(bundle, "features.tsv.gz"), "at") as handle:
                handle.write("partial\n")

            with self.assertRaisesRegex(RuntimeError, "dimensions do not match"):
                validate_resume_inputs(runtime)

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

    def test_manifest_validation_allows_resume_thread_and_tmp_changes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = self.make_runtime(tmpdir)
            runtime.config.update({"num_threads": 8, "performance_opts": {"tmp_root": "/tmp/a"}})
            bam = os.path.join(tmpdir, "sample.filtered.tagged.umi.Aligned.out.bam")
            with open(bam, "wb") as handle:
                handle.write(b"BAM")

            with mock.patch("mfsflow.stage_state._quickcheck_bams"):
                record_stage_success(runtime, MAPPING)
                runtime.config["num_threads"] = 20
                runtime.config["performance_opts"]["tmp_root"] = "/tmp/b"
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
