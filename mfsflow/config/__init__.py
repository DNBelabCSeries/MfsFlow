"""Configuration helpers for MfsFlow."""

from importlib import import_module

__all__ = [
    "build_base_config",
    "configure_reference",
    "DEFAULT_CONFIG",
    "discover_fastq_pairs",
    "load_samplesheet",
    "require_supported_python",
    "resolve_samplesheet_barcodes",
    "resolve_samplesheet_fastq_groups",
    "validate_input_files",
    "write_run_config",
]


def __getattr__(name):
    if name in {
        "build_base_config",
        "configure_reference",
        "discover_fastq_pairs",
        "load_samplesheet",
        "resolve_samplesheet_barcodes",
        "resolve_samplesheet_fastq_groups",
    }:
        builder = import_module("mfsflow.config.builder")
        return getattr(builder, name)
    if name == "write_run_config":
        serialization = import_module("mfsflow.config.serialization")
        write_run_config = serialization.write_run_config
        return write_run_config
    if name == "DEFAULT_CONFIG":
        defaults = import_module("mfsflow.config.defaults")
        DEFAULT_CONFIG = defaults.DEFAULT_CONFIG
        return DEFAULT_CONFIG
    if name in {"require_supported_python", "validate_input_files"}:
        validation = import_module("mfsflow.config.validation")
        return getattr(validation, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
