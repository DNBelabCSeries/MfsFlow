import gzip
import os
import subprocess
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from mfsflow.preflight import (
    check_disk_space,
    check_external_tools,
    check_python_dependencies,
    check_reference_integrity,
    run_preflight,
)
from mfsflow.stage_state import validate_resume_inputs
from mfsflow.stages import COUNTING, FILTERING, MAPPING, SUMMARISING


class PreflightTests(unittest.TestCase):
    def test_pigz_compatibility_fallback_supports_current_commands(self):
        compatibility = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "mfsflow",
            "scripts",
            "pigz_compat.py",
        )
        payload = b"MfsFlow fallback compression test\n"
        compressed = subprocess.check_output(
            [compatibility, "-p", "2", "-c"],
            input=payload,
        )
        self.assertEqual(gzip.decompress(compressed), payload)

        with tempfile.NamedTemporaryFile(suffix=".gz") as handle:
            handle.write(compressed)
            handle.flush()
            decompressed = subprocess.check_output(
                [compatibility, "-p", "2", "-dc", handle.name]
            )
        self.assertEqual(decompressed, payload)

        compressed = subprocess.check_output(
            [compatibility, "--processes=2", "--stdout"],
            input=payload,
        )
        self.assertEqual(gzip.decompress(compressed), payload)

    def test_unusable_pigz_switches_to_compatibility_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = {
                name: os.path.join(tmpdir, name)
                for name in ("samtools", "pigz", "STAR", "featureCounts")
            }
            for name, path in paths.items():
                with open(path, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    if name == "pigz":
                        handle.write("if [ \"$1\" = \"--version\" ]; then exit 0; fi\n")
                        handle.write("echo 'zlib version less than 1.2.3' >&2\nexit 22\n")
                os.chmod(path, 0o755)

            config = {
                "which_Stage": FILTERING,
                "samtools_exec": paths["samtools"],
                "pigz_exec": paths["pigz"],
                "STAR_exec": paths["STAR"],
                "featureCounts_exec": paths["featureCounts"],
            }
            with mock.patch("mfsflow.preflight.shutil.which", return_value=None):
                check_external_tools(config)

            self.assertTrue(config["pigz_exec"].endswith("pigz_compat.py"))

    def test_missing_pigz_uses_compatibility_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = {}
            for name in ("samtools", "STAR", "featureCounts"):
                path = os.path.join(tmpdir, name)
                with open(path, "w") as handle:
                    handle.write("#!/bin/sh\n")
                os.chmod(path, 0o755)
                paths[name] = path
            config = {
                "which_Stage": FILTERING,
                "samtools_exec": paths["samtools"],
                "pigz_exec": os.path.join(tmpdir, "missing-pigz"),
                "STAR_exec": paths["STAR"],
                "featureCounts_exec": paths["featureCounts"],
            }
            with mock.patch("mfsflow.preflight.shutil.which", return_value=None):
                check_external_tools(config)
            self.assertTrue(config["pigz_exec"].endswith("pigz_compat.py"))

    def test_unusable_seqkit_switches_to_gnu_split_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = {}
            for name in ("samtools", "pigz", "seqkit", "STAR", "featureCounts"):
                path = os.path.join(tmpdir, name)
                with open(path, "w") as handle:
                    handle.write("#!/bin/sh\n")
                    if name == "seqkit":
                        handle.write("exit 127\n")
                os.chmod(path, 0o755)
                paths[name] = path
            config = {
                "which_Stage": FILTERING,
                "samtools_exec": paths["samtools"],
                "pigz_exec": paths["pigz"],
                "seqkit_exec": paths["seqkit"],
                "STAR_exec": paths["STAR"],
                "featureCounts_exec": paths["featureCounts"],
            }
            check_external_tools(config)
            self.assertEqual(config["seqkit_exec"], "")

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

    def test_preflight_syncs_resolved_tool_paths_to_runtime(self):
        runtime = SimpleNamespace(
            tools=SimpleNamespace(samtools="old-samtools", pigz="old-pigz", seqkit="old-seqkit")
        )
        config = {
            "samtools_exec": "new-samtools",
            "pigz_exec": "new-pigz",
            "seqkit_exec": "new-seqkit",
        }
        with mock.patch("mfsflow.preflight.check_python_dependencies"), \
             mock.patch("mfsflow.preflight.check_external_tools", side_effect=lambda value: value.update({
                 "samtools_exec": "new-samtools",
                 "pigz_exec": "new-pigz",
                 "seqkit_exec": "new-seqkit",
             })), \
             mock.patch("mfsflow.preflight.check_reference_integrity"), \
             mock.patch("mfsflow.preflight.check_disk_space"), \
             mock.patch("mfsflow.preflight.validate_resume_inputs"):
            run_preflight(config, runtime)

        self.assertEqual(runtime.tools.samtools, "new-samtools")
        self.assertEqual(runtime.tools.pigz, "new-pigz")
        self.assertEqual(runtime.tools.seqkit, "new-seqkit")

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
