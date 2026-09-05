"""Pipeline timing utilities."""

import os
import sys
import time
from contextlib import contextmanager
from datetime import datetime

try:
    import resource
except ImportError:  # pragma: no cover - unavailable on some platforms
    resource = None

from mfsflow.logging_utils import log_error, log_info


def format_duration(seconds):
    """Format a duration in seconds to a human-readable string."""
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.2f}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{sec:05.2f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes):02d}m{sec:05.2f}s"


def resource_usage_details():
    """Return process and waited-child CPU/RSS metrics for timing logs."""
    if resource is None:
        return ""
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    if sys.platform == "darwin":
        self_rss = self_usage.ru_maxrss / (1024 * 1024)
        child_rss = child_usage.ru_maxrss / (1024 * 1024)
    else:
        self_rss = self_usage.ru_maxrss / 1024
        child_rss = child_usage.ru_maxrss / 1024
    peak_rss_mb = max(self_rss, child_rss)
    cpu_user_sec = self_usage.ru_utime + child_usage.ru_utime
    cpu_sys_sec = self_usage.ru_stime + child_usage.ru_stime
    return (
        f"rss_peak_mb={peak_rss_mb:.1f};"
        f"cpu_user_sec={cpu_user_sec:.2f};"
        f"cpu_sys_sec={cpu_sys_sec:.2f}"
    )


class PipelineTimer:
    """Timer for recording pipeline stage execution times."""

    def __init__(self, timing_path, project):
        self.timing_path = timing_path
        self.project = project
        self._ensure_header()

    def _ensure_header(self):
        os.makedirs(os.path.dirname(self.timing_path), exist_ok=True)
        if not os.path.exists(self.timing_path) or os.path.getsize(self.timing_path) == 0:
            with open(self.timing_path, "w") as handle:
                handle.write("timestamp\tproject\tstage\tstatus\tduration_sec\tduration_human\tdetails\n")

    def record(self, stage, status, duration, details=""):
        safe_details = str(details or "").replace("\t", " ").replace("\n", " ")
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(self.timing_path, "a") as handle:
            handle.write(
                f"{ts}\t{self.project}\t{stage}\t{status}\t{duration:.3f}\t"
                f"{format_duration(duration)}\t{safe_details}\n"
            )

    @contextmanager
    def section(self, stage, details=""):
        start = time.perf_counter()
        try:
            yield
        except Exception:
            duration = time.perf_counter() - start
            extra = resource_usage_details()
            self.record(stage, "failed", duration, ";".join(filter(None, (details, extra))))
            log_error(f"Failed {stage} (Duration: {format_duration(duration)})")
            raise
        else:
            duration = time.perf_counter() - start
            extra = resource_usage_details()
            self.record(stage, "ok", duration, ";".join(filter(None, (details, extra))))
            log_info(f"Finished {stage} (Duration: {format_duration(duration)})")
