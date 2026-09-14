"""Runtime preflight checks for dependencies, tools, disk, and references."""

import gzip
import importlib.util
import logging
import os
import shutil
import subprocess
import zlib

from mfsflow.stage_state import validate_resume_inputs
from mfsflow.stages import COUNTING, FILTERING, MAPPING, STAGE_ORDER, SUMMARISING


PYTHON_DEPENDENCIES = {
    "yaml": "PyYAML",
    "pysam": "pysam",
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "anndata": "anndata",
    "h5py": "h5py",
    "PIL": "Pillow",
}

logger = logging.getLogger(__name__)


TOOL_VERSION_ARGS = {
    "samtools": ("--version",),
    "pigz": ("--version",),
    "seqkit": ("version",),
    "STAR": ("--version",),
    "featureCounts": ("-v",),
}


def _as_bool(value, default=False):
    """Interpret YAML booleans and their common string representations."""
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(value)


def check_python_dependencies(extra=None, config=None):
    """Raise an actionable error for packages needed by the requested run."""
    dependencies = {"yaml": PYTHON_DEPENDENCIES["yaml"]}
    if config is None:
        dependencies.update(PYTHON_DEPENDENCIES)
    else:
        stage = config.get("which_stage", config.get("which_Stage", FILTERING))
        start_index = STAGE_ORDER.index(stage)
        remaining = set(STAGE_ORDER[start_index:])
        if remaining.intersection((FILTERING, MAPPING, COUNTING)):
            dependencies["pysam"] = PYTHON_DEPENDENCIES["pysam"]
        if COUNTING in remaining:
            dependencies["numpy"] = PYTHON_DEPENDENCIES["numpy"]
        elif SUMMARISING in remaining and _as_bool(config.get("make_stats", True), default=True):
            dependencies["numpy"] = PYTHON_DEPENDENCIES["numpy"]
        if FILTERING in remaining and config.get("barcode_source") != "samplesheet_barcode":
            dependencies.update({name: PYTHON_DEPENDENCIES[name] for name in ("numpy", "pandas")})
        if _as_bool(config.get("make_h5ad", True), default=True) and COUNTING in remaining:
            dependencies.update({
                name: PYTHON_DEPENDENCIES[name]
                for name in ("pandas", "scipy", "anndata", "h5py")
            })
    dependencies.update(extra or {})
    missing = [package for module, package in dependencies.items() if importlib.util.find_spec(module) is None]
    if missing:
        raise RuntimeError(
            "Missing required Python packages: "
            + ", ".join(sorted(missing))
            + ". Install project dependencies with: python3 -m pip install -r requirements.txt"
        )


def _tool_available(command):
    command = str(command or "")
    if not command:
        return False
    if os.path.isabs(command) or os.path.sep in command:
        return os.path.isfile(command) and os.access(command, os.X_OK)
    return shutil.which(command) is not None


def _probe_tool(name, command):
    """Execute a lightweight tool probe to catch runtime/linker errors."""
    args = [str(command)] + list(TOOL_VERSION_ARGS[name])
    try:
        result = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"External tool {name} cannot execute ({command}): {exc}") from exc
    if result.returncode != 0:
        details = (result.stdout or "").strip()
        raise RuntimeError(
            f"External tool {name} version check failed with code {result.returncode} ({command}): {details}"
        )
    lines = [line.strip() for line in (result.stdout or "").splitlines() if line.strip()]
    if name == "pigz":
        _probe_pigz_runtime(command)
    return lines[0][:500] if lines else "version command succeeded"


def _probe_pigz_runtime(command):
    """Verify that pigz can actually compress data, not just print its version."""
    try:
        result = subprocess.run(
            [str(command), "-c"],
            input=b"MfsFlow pigz runtime probe\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"External tool pigz cannot run a compression probe ({command}): {exc}") from exc
    if result.returncode != 0:
        output = result.stdout or b""
        details = (
            output.decode("utf-8", errors="replace")
            if isinstance(output, bytes)
            else str(output)
        ).strip()
        raise RuntimeError(
            f"External tool pigz compression probe failed with code {result.returncode} "
            f"({command}): {details}"
        )


def _pigz_fallback(command):
    """Find a working system pigz or the packaged pure-Python fallback."""
    system_pigz = shutil.which("pigz")
    if system_pigz and os.path.realpath(system_pigz) != os.path.realpath(str(command)):
        try:
            _probe_tool("pigz", system_pigz)
            return system_pigz
        except RuntimeError:
            pass

    compatibility = os.path.join(os.path.dirname(__file__), "scripts", "pigz_compat.py")
    if os.path.isfile(compatibility):
        try:
            from mfsflow.tool_runtime import ensure_executable_file

            return ensure_executable_file(compatibility, tool_name="pigz_compat")
        except (OSError, RuntimeError):
            pass
    return None


def check_external_tools(config):
    stage = config.get("which_stage", config.get("which_Stage", FILTERING))
    start_index = STAGE_ORDER.index(stage)
    remaining = set(STAGE_ORDER[start_index:])
    required = {}
    if FILTERING in remaining:
        required.update({"samtools": config.get("samtools_exec", "samtools"), "pigz": config.get("pigz_exec", "pigz")})
    if MAPPING in remaining:
        required.update({"STAR": config.get("STAR_exec", "STAR"), "samtools": config.get("samtools_exec", "samtools")})
    if COUNTING in remaining:
        required.update({
            "featureCounts": config.get("featureCounts_exec", "featureCounts"),
            "samtools": config.get("samtools_exec", "samtools"),
        })
    missing = []
    for name, path in list(required.items()):
        if _tool_available(path):
            continue
        if name == "pigz":
            fallback = _pigz_fallback(path)
            if fallback is not None:
                config["pigz_exec"] = fallback
                required[name] = fallback
                logger.warning("Native pigz is unavailable; using fallback executable: %s", fallback)
                continue
        missing.append(f"{name} ({path})")
    if missing:
        raise RuntimeError("Missing or non-executable external tools: " + ", ".join(missing))
    versions = dict(config.get("tool_versions") or {})
    for name, path in list(required.items()):
        try:
            versions[name] = _probe_tool(name, path)
        except RuntimeError as exc:
            if name != "pigz":
                raise
            fallback = _pigz_fallback(path)
            if fallback is None:
                raise RuntimeError(
                    f"Bundled pigz is not usable ({exc}). Install a compatible pigz or set "
                    "performance_opts.tool_cache to a writable location."
                ) from exc
            config["pigz_exec"] = fallback
            required[name] = fallback
            versions[name] = _probe_tool(name, fallback)
            logger.warning("Native pigz is unusable; using fallback executable: %s", fallback)
    if FILTERING in remaining:
        seqkit = config.get("seqkit_exec", "seqkit")
        if not _tool_available(seqkit):
            logger.warning("SeqKit is unavailable (%s); FASTQ splitting will use the GNU split fallback.", seqkit)
            config["seqkit_exec"] = ""
            versions["seqkit"] = "unavailable (GNU split fallback)"
        else:
            try:
                versions["seqkit"] = _probe_tool("seqkit", seqkit)
            except RuntimeError as exc:
                # SeqKit is an optimization, not a correctness requirement;
                # split_fastq already has a GNU split implementation.
                logger.warning(
                    "SeqKit is not usable (%s); FASTQ splitting will use the GNU split fallback: %s",
                    seqkit,
                    exc,
                )
                config["seqkit_exec"] = ""
                versions["seqkit"] = "unusable (GNU split fallback)"
    config["tool_versions"] = versions
    return versions


def check_reference_integrity(config):
    """Validate GTF readability and the core STAR index files before Mapping."""
    stage = config.get("which_stage", config.get("which_Stage", FILTERING))
    start_index = STAGE_ORDER.index(stage)
    remaining = STAGE_ORDER[start_index:]
    if not any(stage in remaining for stage in (FILTERING, MAPPING, COUNTING)):
        return

    reference = config.get("reference", {})
    if MAPPING in remaining:
        star_index = reference.get("STAR_index", "")
        required_index_files = (
            "Genome",
            "SA",
            "SAindex",
            "chrLength.txt",
            "chrName.txt",
            "chrStart.txt",
            "genomeParameters.txt",
        )
        missing = [
            os.path.join(star_index, name)
            for name in required_index_files
            if not os.path.isfile(os.path.join(star_index, name)) or os.path.getsize(os.path.join(star_index, name)) == 0
        ]
        if missing:
            raise RuntimeError("STAR index is incomplete; missing or empty files: " + ", ".join(missing))

    gtf = reference.get("GTF_file", "")
    if not os.path.isfile(gtf) or os.path.getsize(gtf) == 0:
        raise RuntimeError(f"GTF file is missing or empty: {gtf}")
    opener = gzip.open if str(gtf).endswith(".gz") else open
    try:
        with opener(gtf, "rt") as handle:
            if not any(line.strip() and not line.startswith("#") for line in handle):
                raise RuntimeError(f"GTF file contains no annotation records: {gtf}")
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise RuntimeError(f"GTF file is unreadable: {gtf}: {exc}") from exc


def _input_size(config):
    total = 0
    for group in ("file1", "file2"):
        names = ((config.get("sequence_files") or {}).get(group) or {}).get("name", "")
        for path in str(names).split(","):
            path = path.strip()
            if path and os.path.isfile(path):
                total += os.path.getsize(path)
    return total


def check_disk_space(config, runtime):
    """Check free space on output and temporary filesystems."""
    opts = config.get("performance_opts", {}) or {}
    minimum = float(opts.get("min_free_gb", 5) or 0) * (1024 ** 3)
    multiplier = float(opts.get("disk_space_multiplier", 4.0) or 0)
    estimated = _input_size(config) * multiplier if runtime.which_stage == FILTERING else 0
    required = max(minimum, estimated)
    checked_devices = set()
    failures = []
    for path in (runtime.out_dir, runtime.tmp_merge_path):
        os.makedirs(path, exist_ok=True)
        device = os.stat(path).st_dev
        if device in checked_devices:
            continue
        checked_devices.add(device)
        free = shutil.disk_usage(path).free
        if free < required:
            failures.append(f"{path}: {free / (1024 ** 3):.1f} GiB free, {required / (1024 ** 3):.1f} GiB required")
    if failures:
        raise RuntimeError("Insufficient disk space: " + "; ".join(failures))


def run_preflight(config, runtime):
    """Run all checks required immediately before pipeline execution."""
    check_python_dependencies(config=config)
    check_external_tools(config)
    # Tool resolution can change a native executable to a validated fallback
    # during preflight. Keep the already-created runtime in sync before stages
    # read runtime.tools.
    for config_key, runtime_name in (
        ("samtools_exec", "samtools"),
        ("pigz_exec", "pigz"),
        ("seqkit_exec", "seqkit"),
    ):
        if hasattr(runtime, "tools") and config_key in config:
            setattr(runtime.tools, runtime_name, config[config_key])
    check_reference_integrity(config)
    check_disk_space(config, runtime)
    validate_resume_inputs(runtime)
