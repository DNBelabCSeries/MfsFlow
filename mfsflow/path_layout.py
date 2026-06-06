"""
Directory layout helpers for organizing pipeline output directories.

This module provides functions to construct standard directory paths
for configuration, logs, barcodes, expression data, statistics,
intermediate files, and temporary merge directories within the
pipeline output structure.
"""

import os


def config_dir(out_dir):
    """Return the path to the configuration directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'config' subdirectory.
    """
    return os.path.join(out_dir, "config")


def logs_dir(out_dir):
    """Return the path to the logs directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'logs' subdirectory.
    """
    return os.path.join(out_dir, "logs")


def stage_state_dir(out_dir):
    """Return the directory containing stage manifests and success markers."""
    return os.path.join(logs_dir(out_dir), "stages")


def barcode_dir(out_dir):
    """Return the path to the barcodes directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'barcodes' subdirectory.
    """
    return os.path.join(out_dir, "barcodes")


def expression_dir(out_dir):
    """Return the path to the expression data directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'expression' subdirectory.
    """
    return os.path.join(out_dir, "expression")


def stats_dir(out_dir):
    """Return the path to the statistics directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'stats' subdirectory.
    """
    return os.path.join(out_dir, "stats")


def intermediate_dir(out_dir):
    """Return the path to the intermediate files directory within the output tree.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'intermediate' subdirectory.
    """
    return os.path.join(out_dir, "intermediate")


def tmp_merge_dir(out_dir):
    """Return the path to the temporary merge directory within the intermediate files.
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'tmp_merge' subdirectory within intermediate files.
    """
    return os.path.join(intermediate_dir(out_dir), "tmp_merge")


def outputs_dir(out_dir):
    """Return the path to the final outputs directory (sibling of out_dir).
    
    Args:
        out_dir (str): Root output directory path.
        
    Returns:
        str: Full path to the 'outs' directory at the same level as out_dir.
    """
    return os.path.join(os.path.dirname(out_dir), "outs")


def ensure_layout(out_dir):
    """Create the complete directory layout for pipeline output.
    
    Ensures all standard directories exist within the output tree,
    creating them if necessary. This includes the root output directory
    and all subdirectories for config, logs, barcodes, expression data,
    statistics, intermediate files, and final outputs.
    
    Args:
        out_dir (str): Root output directory path.
    """
    for path in (
        out_dir,
        config_dir(out_dir),
        logs_dir(out_dir),
        stage_state_dir(out_dir),
        barcode_dir(out_dir),
        expression_dir(out_dir),
        stats_dir(out_dir),
        intermediate_dir(out_dir),
        tmp_merge_dir(out_dir),
        outputs_dir(out_dir),
    ):
        os.makedirs(path, exist_ok=True)
