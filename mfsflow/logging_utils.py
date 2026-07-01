"""Logging helpers used by pipeline entry points and stages."""

import logging
import sys


logger = logging.getLogger("mfsflow")
logger.propagate = False

_formatter = logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class _StdoutFilter(logging.Filter):
    def filter(self, record):
        return record.levelno < logging.ERROR


def _configure_logger():
    if logger.handlers:
        return

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.setFormatter(_formatter)
    stdout_handler.addFilter(_StdoutFilter())
    logger.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.ERROR)
    stderr_handler.setFormatter(_formatter)
    logger.addHandler(stderr_handler)

    logger.setLevel(logging.INFO)


_configure_logger()


def log_info(msg):
    """Log an informational message with timestamp to stdout."""
    logger.info(msg)


def log_error(msg):
    """Log an error message with timestamp to stderr."""
    logger.error(msg)


class Tee:
    """Tee-like stream multiplexer that writes to multiple streams."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, s):
        for stream in self._streams:
            stream.write(s)

    def flush(self):
        for stream in self._streams:
            stream.flush()
