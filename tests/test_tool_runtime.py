import os
import stat
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

from mfsflow.runtime import PipelineRuntime
from mfsflow.tool_runtime import ensure_executable_file, resolve_config_tool_paths


class ToolRuntimeTests(unittest.TestCase):
    def test_missing_execute_bit_is_repaired_in_place(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tool = Path(tmpdir) / "samtools"
            tool.write_bytes(b"tool")
            tool.chmod(stat.S_IRUSR | stat.S_IWUSR)

            resolved = ensure_executable_file(tool)

            self.assertEqual(resolved, str(tool))
            self.assertTrue(os.access(str(tool), os.X_OK))

    def test_read_only_tool_uses_cache_copy(self):
        with tempfile.TemporaryDirectory() as source_dir, tempfile.TemporaryDirectory() as cache_dir:
            tool = Path(source_dir) / "STAR"
            tool.write_bytes(b"tool")
            tool.chmod(stat.S_IRUSR | stat.S_IWUSR)

            with mock.patch("mfsflow.tool_runtime.Path.chmod", side_effect=PermissionError("read-only")):
                resolved = ensure_executable_file(tool, cache_dir=cache_dir)

            cached = Path(resolved)
            self.assertNotEqual(cached, tool)
            self.assertEqual(cached.read_bytes(), b"tool")
            self.assertTrue(os.access(str(cached), os.X_OK))

    def test_cache_is_reused_when_source_metadata_matches(self):
        with tempfile.TemporaryDirectory() as source_dir, tempfile.TemporaryDirectory() as cache_dir:
            tool = Path(source_dir) / "pigz"
            tool.write_bytes(b"tool")
            tool.chmod(stat.S_IRUSR | stat.S_IWUSR)

            with mock.patch("mfsflow.tool_runtime.Path.chmod", side_effect=PermissionError("read-only")):
                first = ensure_executable_file(tool, cache_dir=cache_dir)
            with mock.patch("mfsflow.tool_runtime.Path.chmod", side_effect=PermissionError("read-only")):
                second = ensure_executable_file(tool, cache_dir=cache_dir)

            self.assertEqual(first, second)

    def test_config_paths_resolve_bundled_tools(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            software = Path(tmpdir) / "software"
            software.mkdir()
            tool = software / "samtools"
            tool.write_bytes(b"tool")
            tool.chmod(0o755)

            config = {"samtools_exec": "samtools"}
            with mock.patch("mfsflow.tool_runtime.sys.platform", "linux"), \
                 mock.patch("mfsflow.tool_runtime.platform.machine", return_value="x86_64"):
                resolved = resolve_config_tool_paths(config, tmpdir)

            self.assertEqual(resolved["samtools_exec"], str(tool))

    def test_resumed_run_persists_repaired_tool_path_for_children(self):
        try:
            import yaml  # noqa: F401
        except ImportError:
            self.skipTest("PyYAML is not installed in the lightweight local environment")
        if not isinstance(yaml, types.ModuleType):
            self.skipTest("PyYAML import is mocked in the lightweight local test suite")

        with tempfile.TemporaryDirectory() as tmpdir:
            software = Path(tmpdir) / "software"
            software.mkdir()
            tool = software / "samtools"
            tool.write_bytes(b"tool")
            tool.chmod(stat.S_IRUSR | stat.S_IWUSR)
            yaml_file = Path(tmpdir) / "run_config.yaml"
            yaml_file.write_text(
                "project: P1\n"
                f"out_dir: {tmpdir}/XPRESS_PROCESSING\n"
                "num_threads: 1\n"
                "which_Stage: Filtering\n"
                f"toolkit_directory: {tmpdir}\n"
                f"samtools_exec: {tool}\n",
                encoding="utf-8",
            )
            config = {
                "project": "P1",
                "out_dir": os.path.join(tmpdir, "XPRESS_PROCESSING"),
                "num_threads": 1,
                "which_Stage": "Filtering",
                "toolkit_directory": tmpdir,
                "samtools_exec": str(tool),
                "performance_opts": {"tool_cache": os.path.join(tmpdir, "tool-cache")},
            }

            with mock.patch("mfsflow.tool_runtime.sys.platform", "linux"), \
                 mock.patch("mfsflow.tool_runtime.platform.machine", return_value="x86_64"), \
                 mock.patch("mfsflow.tool_runtime.Path.chmod", side_effect=PermissionError("read-only")):
                PipelineRuntime.from_config(config, str(yaml_file))

            config_text = yaml_file.read_text(encoding="utf-8")
            self.assertIn("tool-cache", config_text)
            self.assertNotIn(f"samtools_exec: {tool}", config_text)

    def test_non_linux_platform_does_not_select_packaged_elf_tools(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            software = Path(tmpdir) / "software"
            software.mkdir()
            tool = software / "samtools"
            tool.write_bytes(b"ELF")
            tool.chmod(0o755)

            with mock.patch("mfsflow.tool_runtime.sys.platform", "darwin"), \
                 mock.patch("mfsflow.tool_runtime.platform.machine", return_value="arm64"):
                resolved = resolve_config_tool_paths({"samtools_exec": "samtools"}, tmpdir)

            self.assertEqual(resolved["samtools_exec"], "samtools")

    def test_unrepairable_explicit_path_falls_back_to_path_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tool = Path(tmpdir) / "broken-samtools"
            tool.write_bytes(b"tool")

            with mock.patch(
                "mfsflow.tool_runtime.ensure_executable_file",
                side_effect=RuntimeError("read-only"),
            ), mock.patch("mfsflow.tool_runtime.shutil.which", return_value="/usr/bin/samtools"):
                resolved = resolve_config_tool_paths({"samtools_exec": str(tool)}, tmpdir)

            self.assertEqual(resolved["samtools_exec"], "samtools")


if __name__ == "__main__":
    unittest.main()
