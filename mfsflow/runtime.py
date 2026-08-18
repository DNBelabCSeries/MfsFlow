"""Pipeline runtime environment: tool resolution and runtime dataclasses."""

import hashlib
import os
import sys
from dataclasses import dataclass

from mfsflow.path_layout import logs_dir, tmp_merge_dir
from mfsflow.tool_runtime import persist_tool_paths, resolve_config_tool_paths

__all__ = [
    "PipelineRuntime",
    "PipelineTools",
]


@dataclass
class PipelineTools:
    """Data class holding executable paths for external tools."""
    samtools: str
    pigz: str
    seqkit: str


@dataclass
class PipelineRuntime:
    """Runtime configuration and environment for the pipeline.
    
    This data class encapsulates all runtime parameters needed for pipeline
    execution, including configuration, paths, tools, and environment variables.
    """
    config: dict
    yaml_file: str
    project: str
    out_dir: str
    num_threads: int
    which_stage: str
    python_exec: str
    toolkit_dir: str
    tools: PipelineTools
    exec_env: dict
    analysis_dir: str
    log_path: str
    timing_path: str
    tmp_merge_path: str

    @classmethod
    def from_config(cls, config, yaml_file):
        """Create a PipelineRuntime instance from configuration dictionary.
        
        Args:
            config (dict): Pipeline configuration dictionary.
            yaml_file (str): Path to the YAML configuration file.
            
        Returns:
            PipelineRuntime: Configured runtime instance.
        """
        toolkit_dir = config.get("toolkit_directory")
        if not toolkit_dir:
            # Default to mfsflow package location
            toolkit_dir = os.path.dirname(os.path.abspath(__file__))
        elif not os.path.isabs(toolkit_dir):
            # Relative path: resolve relative to YAML file's directory
            toolkit_dir = os.path.join(os.path.dirname(os.path.abspath(yaml_file)), toolkit_dir)
        cache_dir = (config.get("performance_opts", {}) or {}).get("tool_cache")
        resolved_tools = resolve_config_tool_paths(config, toolkit_dir, cache_dir=cache_dir)
        persist_tool_paths(yaml_file, config, resolved_tools)
        config.update(resolved_tools)
        exec_env = os.environ.copy()
        software_dir = os.path.join(toolkit_dir, "software")
        if sys.platform.startswith("linux") and os.path.isdir(software_dir):
            exec_env["PATH"] = software_dir + os.pathsep + exec_env.get("PATH", "")

        out_dir = config["out_dir"]
        tmp_merge_path = cls._resolve_tmp_merge_path(config, out_dir)
        which_stage = config.get("which_stage", config.get("which_Stage"))
        if not which_stage:
            raise KeyError("Missing pipeline stage key: which_Stage")
        try:
            num_threads = int(config["num_threads"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("num_threads must be a positive integer") from exc
        if num_threads < 1:
            raise ValueError("num_threads must be a positive integer")
        os.makedirs(tmp_merge_path, exist_ok=True)
        return cls(
            config=config,
            yaml_file=yaml_file,
            project=config["project"],
            out_dir=out_dir,
            num_threads=num_threads,
            which_stage=which_stage,
            python_exec=sys.executable or "python3",
            toolkit_dir=toolkit_dir,
            tools=PipelineTools(
                samtools=config.get("samtools_exec", "samtools"),
                pigz=config.get("pigz_exec", "pigz"),
                seqkit=config.get("seqkit_exec", "seqkit"),
            ),
            exec_env=exec_env,
            analysis_dir=out_dir,
            log_path=os.path.join(logs_dir(out_dir), "pipeline.log"),
            timing_path=os.path.join(logs_dir(out_dir), "pipeline_timing.tsv"),
            tmp_merge_path=tmp_merge_path,
        )

    @staticmethod
    def _resolve_tmp_merge_path(config, out_dir):
        """Resolve the temporary merge directory path based on configuration.
        
        Args:
            config (dict): Pipeline configuration dictionary.
            out_dir (str): Output directory path.
            
        Returns:
            str: Resolved temporary merge directory path.
        """
        tmp_root = (config.get("performance_opts", {}) or {}).get("tmp_root")
        if not tmp_root:
            return tmp_merge_dir(out_dir)
        project = str(config.get("project", "sample"))
        safe_project = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in project)
        out_hash = hashlib.sha1(os.path.abspath(out_dir).encode("utf-8")).hexdigest()[:10]
        return os.path.join(os.path.abspath(tmp_root), f"mfsflow_{safe_project}_{out_hash}", "tmp_merge")

    def resolve_script(self, script_name):
        """Resolve the path to a pipeline script.
        
        Searches for the script in multiple locations: toolkit directory,
        scripts subdirectory, and mfsflow package directory.
        
        Args:
            script_name (str): Name of the script to resolve.
            
        Returns:
            str: Full path to the script.
            
        Raises:
            FileNotFoundError: If the script cannot be found in any location.
        """
        # mfsflow package location (where runtime.py resides)
        mfsflow_pkg = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(self.toolkit_dir, script_name),
            os.path.join(self.toolkit_dir, "scripts", script_name),
            os.path.join(mfsflow_pkg, "scripts", script_name),
            os.path.join(mfsflow_pkg, script_name),
        ]
        for candidate in candidates:
            if os.path.exists(candidate):
                return candidate
        raise FileNotFoundError(f"Script not found: {script_name}. Tried: {', '.join(candidates)}")
