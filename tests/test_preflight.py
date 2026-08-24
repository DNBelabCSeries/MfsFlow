import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from mfsflow.preflight import (
    check_disk_space,
    check_external_tools,
    check_python_dependencies,
    check_reference_integrity,
)
from mfsflow.stage_state import validate_resume_inputs
from mfsflow.stages import COUNTING, FILTERING, MAPPING, SUMMARISING


class PreflightTests(unittest.TestCase):
    def test_missing_python_dependency_reports_install_command(self):
        with mock.patch("mfsflow.preflight.importlib.util.find_spec", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "pip install -r requirements.txt"):
                check_python_dependencies({"missing_module": "missing-package"})

    def test_counting_without_h5ad_does_not_require_scipy(self):
        def find_spec(name):
            return None if name == "scipy" else object()

        config = {"which_Stage": COUNTING, "make_h5ad": False}
        with mock.patch("mfsflow.preflight.importlib.util.find_spec", side_effect=find_spec):
            check_python_dependencies(config=config)

    def test_counting_with_h5ad_requires_scipy(self):
        def find_spec(name):
            return None if name == "scipy" else object()

        config = {"which_Stage": COUNTING, "make_h5ad": True}
        with mock.patch("mfsflow.preflight.importlib.util.find_spec", side_effect=find_spec):
            with self.assertRaisesRegex(RuntimeError, "scipy"):
                check_python_dependencies(config=config)

    def test_external_tool_check_accepts_executable_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tools = {}
            for name in ("samtools", "pigz", "STAR", "featureCounts"):
                path = os.path.join(tmpdir, name)
                with open(path, "w") as handle:
                    handle.write("#!/bin/sh\n")
                os.chmod(path, 0o755)
                tools[name] = path
            check_external_tools({
                "which_Stage": FILTERING,
                "samtools_exec": tools["samtools"],
                "pigz_exec": tools["pigz"],
                "STAR_exec": tools["STAR"],
                "featureCounts_exec": tools["featureCounts"],
            })

    def test_external_tool_check_rejects_exec_format_errors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tool = os.path.join(tmpdir, "samtools")
            with open(tool, "w") as handle:
                handle.write("not a native executable")
            os.chmod(tool, 0o755)
            config = {
                "which_Stage": FILTERING,
                "samtools_exec": tool,
                "pigz_exec": tool,
                "STAR_exec": tool,
                "featureCounts_exec": tool,
            }
            with mock.patch("mfsflow.preflight.subprocess.run", side_effect=OSError("Exec format error")):
                with self.assertRaisesRegex(RuntimeError, "cannot execute"):
                    check_external_tools(config)

    def test_later_stage_preflight_preserves_recorded_tool_versions(self):
        config = {
            "which_Stage": SUMMARISING,
            "tool_versions": {"STAR": "STAR_2.7", "samtools": "samtools 1.20"},
        }

        versions = check_external_tools(config)

        self.assertEqual(versions, {"STAR": "STAR_2.7", "samtools": "samtools 1.20"})
        self.assertEqual(config["tool_versions"], versions)

    def test_incomplete_star_index_fails_before_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gtf = os.path.join(tmpdir, "genes.gtf")
            with open(gtf, "w") as handle:
                handle.write('chr1\tT\texon\t1\t10\t.\t+\t.\tgene_id "g";\n')
            with self.assertRaisesRegex(RuntimeError, "STAR index is incomplete"):
                check_reference_integrity({
                    "which_Stage": MAPPING,
                    "reference": {"STAR_index": tmpdir, "GTF_file": gtf},
                })

    def test_disk_check_can_use_configured_zero_minimum(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = SimpleNamespace(out_dir=tmpdir, tmp_merge_path=tmpdir, which_stage=SUMMARISING)
            check_disk_space({"performance_opts": {"min_free_gb": 0}}, runtime)

    def test_disk_check_rejects_insufficient_space(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = SimpleNamespace(out_dir=tmpdir, tmp_merge_path=tmpdir, which_stage=SUMMARISING)
            usage = SimpleNamespace(total=10, used=10, free=0)
            with mock.patch("mfsflow.preflight.shutil.disk_usage", return_value=usage):
                with self.assertRaisesRegex(RuntimeError, "Insufficient disk space"):
                    check_disk_space({"performance_opts": {"min_free_gb": 1}}, runtime)

    def test_counting_resume_requires_nonempty_mapping_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = SimpleNamespace(
                which_stage=COUNTING,
                project="sample",
                out_dir=tmpdir,
                tmp_merge_path=os.path.join(tmpdir, "tmp"),
            )
            with self.assertRaisesRegex(RuntimeError, "Cannot resume from Counting"):
                validate_resume_inputs(runtime)

    def test_summarising_resume_rejects_corrupt_matrix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix_dir = os.path.join(tmpdir, "expression", "sample.exon.umi")
            os.makedirs(matrix_dir)
            matrix = os.path.join(matrix_dir, "matrix.mtx.gz")
            with open(matrix, "wb") as handle:
                handle.write(b"not gzip")
            runtime = SimpleNamespace(
                which_stage=SUMMARISING,
                project="sample",
                out_dir=tmpdir,
                tmp_merge_path=os.path.join(tmpdir, "tmp"),
            )
            with self.assertRaisesRegex(RuntimeError, "corrupt or unreadable"):
                validate_resume_inputs(runtime)


if __name__ == "__main__":
    unittest.main()
