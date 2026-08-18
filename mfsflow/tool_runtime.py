"""Resolve bundled executables and recover permissions lost during install.

Wheel installers and application package managers do not all preserve the
execute bit on package data files.  MfsFlow keeps its native tools as package
resources, so this module repairs the bit in place when possible and falls
back to a user-writable cache when the installed package is read-only.
"""

import logging
import os
import shutil
import stat
import tempfile
from pathlib import Path

from mfsflow import __version__


TOOL_SPECS = (
    ("samtools_exec", "samtools"),
    ("pigz_exec", "pigz"),
    ("seqkit_exec", "seqkit"),
    ("STAR_exec", "STAR"),
    ("featureCounts_exec", "featureCounts"),
)

logger = logging.getLogger(__name__)


def _is_executable(path):
    """Return whether *path* is a regular file executable by this process."""
    try:
        return path.is_file() and os.access(str(path), os.X_OK)
    except OSError:
        return False


def _cache_root(cache_dir=None):
    """Return the configured user-writable tool cache directory."""
    configured = cache_dir or os.environ.get("MFSFLOW_TOOL_CACHE")
    if configured:
        return Path(configured).expanduser()

    xdg_cache = os.environ.get("XDG_CACHE_HOME")
    if xdg_cache:
        return Path(xdg_cache).expanduser() / "mfsflow" / "tools"

    try:
        home = Path.home()
    except RuntimeError:
        home = Path(tempfile.gettempdir())
    return home / ".cache" / "mfsflow" / "tools"


def _cache_path(source, cache_dir, tool_name):
    """Build a stable cache path for a bundled tool version."""
    del source  # The source metadata is stored beside the cached file.
    return _cache_root(cache_dir) / str(__version__) / "bin" / tool_name


def _metadata_path(destination):
    return destination.with_name(destination.name + ".meta")


def _cache_is_current(destination, source):
    """Check the cheap source metadata used to invalidate cached binaries."""
    metadata = _metadata_path(destination)
    try:
        source_stat = source.stat()
        if not destination.is_file() or not _is_executable(destination):
            return False
        with metadata.open("r", encoding="ascii") as handle:
            size, mtime_ns = handle.read().strip().split(" ")
        return int(size) == source_stat.st_size and int(mtime_ns) == source_stat.st_mtime_ns
    except (OSError, ValueError):
        return False


def _write_cache_metadata(destination, source):
    source_stat = source.stat()
    metadata = _metadata_path(destination)
    temporary = metadata.with_name(f".{metadata.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="ascii") as handle:
            handle.write(f"{source_stat.st_size} {source_stat.st_mtime_ns}\n")
        os.replace(str(temporary), str(metadata))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _copy_to_cache(source, cache_dir, tool_name):
    """Copy a non-executable bundled tool atomically to the user cache."""
    destination = _cache_path(source, cache_dir, tool_name)
    if _cache_is_current(destination, source):
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{tool_name}.",
            suffix=".tmp",
            dir=str(destination.parent),
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            with source.open("rb") as source_handle:
                shutil.copyfileobj(source_handle, handle)

        source_mode = source.stat().st_mode
        os.chmod(
            str(temporary),
            source_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH,
        )
        os.replace(str(temporary), str(destination))
        _write_cache_metadata(destination, source)
    finally:
        if temporary is not None:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    return destination


def ensure_executable_file(path, tool_name=None, cache_dir=None):
    """Make a file executable, or return an executable cache copy.

    The original path is preferred.  A cache copy is only created when the
    original file exists but its mode cannot be repaired, which avoids copying
    the roughly 50 MiB bundled tool set on every run.
    """
    source = Path(path).expanduser()
    if not source.is_file():
        raise FileNotFoundError(f"Bundled tool is missing: {source}")

    if _is_executable(source):
        return str(source)

    chmod_error = None
    try:
        source.chmod(source.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except OSError as exc:
        chmod_error = exc

    if _is_executable(source):
        logger.warning("Restored execute permission for bundled tool: %s", source)
        return str(source)

    name = tool_name or source.name
    try:
        cached = _copy_to_cache(source, cache_dir, name)
    except (OSError, shutil.Error) as cache_error:
        reason = f"chmod failed: {chmod_error}; cache copy failed: {cache_error}"
        raise RuntimeError(
            f"Bundled tool is not executable and could not be repaired: {source}. "
            f"{reason}. Set MFSFLOW_TOOL_CACHE to a writable directory or run "
            f"chmod +x {source}."
        ) from cache_error

    logger.warning(
        "Bundled tool is not writable; using executable cache copy: %s -> %s",
        source,
        cached,
    )
    return str(cached)


def resolve_bundled_tool(name, toolkit_dir, fallback=None, cache_dir=None):
    """Resolve a bundled tool, preserving the existing PATH fallback."""
    bundled = Path(os.path.abspath(os.path.join(os.path.expanduser(str(toolkit_dir)), "software", name)))
    if not bundled.exists():
        return fallback or name

    try:
        return ensure_executable_file(bundled, tool_name=name, cache_dir=cache_dir)
    except RuntimeError:
        # Keep compatibility with installations that intentionally rely on a
        # system tool when the package resource cannot be repaired or cached.
        if shutil.which(name):
            logger.warning("Falling back to PATH tool %s because bundled tool could not be repaired.", name)
            return fallback or name
        raise


def resolve_config_tool_paths(config, toolkit_dir, cache_dir=None):
    """Resolve all tool paths in *config* and return changed values.

    Existing explicit paths are repaired in place.  Bare command names first
    resolve against ``<toolkit_dir>/software`` and otherwise remain PATH
    commands, matching the historical configuration behavior.
    """
    resolved = {}
    for config_key, tool_name in TOOL_SPECS:
        configured = config.get(config_key)
        configured = str(configured) if configured else ""
        is_path = os.path.isabs(configured) or os.path.sep in configured
        if is_path and os.path.isfile(os.path.expanduser(configured)):
            try:
                resolved[config_key] = ensure_executable_file(
                    configured,
                    tool_name=tool_name,
                    cache_dir=cache_dir,
                )
            except RuntimeError:
                if shutil.which(tool_name):
                    resolved[config_key] = configured
                else:
                    raise
        else:
            fallback = configured or tool_name
            resolved[config_key] = resolve_bundled_tool(
                tool_name,
                toolkit_dir,
                fallback=fallback,
                cache_dir=cache_dir,
            )
    return resolved


def persist_tool_paths(yaml_file, config, resolved):
    """Persist repaired tool paths for child processes in a resumed run."""
    changed = any(config.get(key) != value for key, value in resolved.items())
    if not changed or not yaml_file or not os.path.isfile(yaml_file):
        return

    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required to persist repaired tool paths in the run configuration."
        ) from exc

    temporary = None
    try:
        config_dir = os.path.dirname(os.path.abspath(yaml_file))
        persisted_config = dict(config)
        persisted_config.update(resolved)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".run_config.",
            suffix=".tmp",
            dir=config_dir,
            delete=False,
        ) as handle:
            temporary = handle.name
            yaml.safe_dump(persisted_config, handle, default_flow_style=False, sort_keys=False)
        os.replace(temporary, yaml_file)
    except (OSError, yaml.YAMLError) as exc:
        try:
            os.unlink(temporary)
        except (UnboundLocalError, FileNotFoundError):
            pass
        raise RuntimeError(
            f"Resolved tool paths but could not update run configuration {yaml_file}: {exc}"
        ) from exc
