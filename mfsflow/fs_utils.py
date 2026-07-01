"""Filesystem helpers shared by pipeline stages."""

import os


def remove_path(path):
    """Remove a file and its associated .bai index file if they exist."""
    if not path:
        return
    if os.path.exists(path):
        os.remove(path)
    bai = path + ".bai"
    if os.path.exists(bai):
        os.remove(bai)
