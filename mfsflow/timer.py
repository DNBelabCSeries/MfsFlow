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


def resource_usage_snapshot():
    """Capture cumulative usage counters used to calculate stage deltas."""
    if resource is None:
        return None
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    child_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    if sys.platform == "darwin":
        self_rss = self_usage.ru_maxrss / (1024 * 1024)
        child_rss = child_usage.ru_maxrss / (1024 * 1024)
    else:
        self_rss = self_usage.ru_maxrss / 1024
        child_rss = child_usage.ru_maxrss / 1024
    return {
        "self_user": self_usage.ru_utime,
        "self_sys": self_usage.ru_stime,
        "child_user": child_usage.ru_utime,
        "child_sys": child_usage.ru_stime,
        "self_rss": self_rss,
        "child_rss": child_rss,
    }


def resource_usage_details(baseline=None):
    """Return stage CPU deltas and process/child RSS high-water metrics.

    ``getrusage`` exposes a lifetime high-water RSS value rather than a
    per-stage peak.  Keep that limitation explicit in the field name while
    reporting CPU as a delta when a baseline is supplied.
    """
    current = resource_usage_snapshot()
    if current is None:
        return ""
    baseline = baseline or {
        "self_user": 0.0,
        "self_sys": 0.0,
        "child_user": 0.0,
        "child_sys": 0.0,
    }
    peak_rss_mb = max(current["self_rss"], current["child_rss"])
    cpu_user_sec = max(
        0.0,
        (current["self_user"] - baseline.get("self_user", 0.0))
        + (current["child_user"] - baseline.get("child_user", 0.0)),
    )
    cpu_sys_sec = max(
        0.0,
        (current["self_sys"] - baseline.get("self_sys", 0.0))
        + (current["child_sys"] - baseline.get("child_sys", 0.0)),
    )
    return (
        # Keep the historical key for consumers of pipeline_timing.tsv; the
        # explicit scope below documents that it is not a per-stage sample.
        f"rss_peak_mb={peak_rss_mb:.1f};"
        f"cpu_user_sec={cpu_user_sec:.2f};"
        f"cpu_sys_sec={cpu_sys_sec:.2f};"
        f"rss_scope=process_lifetime_highwater"
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
        usage_start = resource_usage_snapshot()
        try:
            yield
        except Exception:
            duration = time.perf_counter() - start
            extra = resource_usage_details(usage_start)
            self.record(stage, "failed", duration, ";".join(filter(None, (details, extra))))
            log_error(f"Failed {stage} (Duration: {format_duration(duration)})")
            raise
        else:
            duration = time.perf_counter() - start
            extra = resource_usage_details(usage_start)
            self.record(stage, "ok", duration, ";".join(filter(None, (details, extra))))
            log_info(f"Finished {stage} (Duration: {format_duration(duration)})")
