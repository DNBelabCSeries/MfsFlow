"""
Serialization of pipeline configuration to a YAML run config file with tool path resolution.

This module handles writing the final pipeline configuration to a YAML file,
resolving tool executable paths and applying custom YAML serialization for
boolean and None values.
"""

import copy
import os
import tempfile

import yaml

from mfsflow.tool_runtime import resolve_bundled_tool


def write_run_config(config):
    """Write the pipeline configuration to a YAML run config file.
    
    Resolves tool executable paths and serializes the configuration
    with custom YAML formatting for booleans, None values, and strings.
    
    Args:
        config (dict): Pipeline configuration dictionary.
        
    Returns:
        str: Path to the generated YAML configuration file.
    """
    run_config = copy.deepcopy(config)
    toolkit_dir = config.get("toolkit_directory") or "."
    run_config["toolkit_directory"] = toolkit_dir
    cache_dir = (config.get("performance_opts", {}) or {}).get("tool_cache")

    def resolve_tool(name):
        return resolve_bundled_tool(name, toolkit_dir, fallback=name, cache_dir=cache_dir)

    run_config["samtools_exec"] = resolve_tool("samtools")
    run_config["pigz_exec"] = resolve_tool("pigz")
    run_config["seqkit_exec"] = resolve_tool("seqkit")
    run_config["STAR_exec"] = resolve_tool("STAR")
    run_config["featureCounts_exec"] = resolve_tool("featureCounts")

    class ForceStr:
        def __init__(self, value):
            self.value = value

    def force_str_representer(dumper, data):
        return dumper.represent_scalar("tag:yaml.org,2002:str", data.value, style='"')

    class RunConfigDumper(yaml.SafeDumper):
        pass

    def bool_representer(dumper, value):
        return dumper.represent_scalar("tag:yaml.org,2002:bool", "yes" if value else "no")

    def none_representer(dumper, _value):
        return dumper.represent_scalar("tag:yaml.org,2002:null", "~")

    yaml.add_representer(ForceStr, force_str_representer, Dumper=RunConfigDumper)
    yaml.add_representer(bool, bool_representer, Dumper=RunConfigDumper)
    yaml.add_representer(type(None), none_representer, Dumper=RunConfigDumper)

    final_yaml_path = os.path.join(config["out_dir"], "config", "run_config.yaml")
    os.makedirs(os.path.dirname(final_yaml_path), exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".run_config.",
            suffix=".tmp",
            dir=os.path.dirname(final_yaml_path),
            delete=False,
        ) as handle:
            temporary = handle.name
            yaml.dump(run_config, handle, Dumper=RunConfigDumper, default_flow_style=False, sort_keys=False)
        os.replace(temporary, final_yaml_path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
    return final_yaml_path


def load_run_config(path):
    """Load and validate an existing run configuration."""
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Run configuration is empty or invalid: {path}")
    return config


def persist_run_config(config, path):
    """Atomically persist runtime state without re-resolving tools or inputs."""
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".run_config.",
            suffix=".tmp",
            dir=directory,
            delete=False,
        ) as handle:
            temporary = handle.name
            yaml.safe_dump(config, handle, default_flow_style=False, sort_keys=False)
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
    return path
