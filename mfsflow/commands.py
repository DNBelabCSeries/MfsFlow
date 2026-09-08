"""Subprocess execution helpers for pipeline stages."""

import subprocess
import time
from collections import deque

from mfsflow.logging_utils import log_error, log_info
from mfsflow.timer import format_duration, resource_usage_details, resource_usage_snapshot


def run_stage_cmd(cmd, stage_name, run_log, exec_env, timer, log_path, shell=False):
    """Execute a pipeline stage command with logging and timing."""
    if isinstance(cmd, list) and not shell:
        cmd_str = " ".join(map(str, cmd))
    else:
        cmd_str = str(cmd)
    start = time.perf_counter()
    usage_start = resource_usage_snapshot()
    res = subprocess.run(cmd, stdout=run_log, stderr=subprocess.STDOUT, shell=shell, env=exec_env)
    duration = time.perf_counter() - start
    resource_details = resource_usage_details(usage_start)
    if res.returncode != 0:
        timer.record(
            stage_name,
            "failed",
            duration,
            ";".join(filter(None, (cmd_str, resource_details))),
        )
        log_error(f"Failed {stage_name} (Duration: {format_duration(duration)})")
        run_log.flush()
        try:
            with open(log_path, "r") as lr:
                tail = "".join(deque(lr, maxlen=30))
                log_error(f"{stage_name} failed (rc={res.returncode}). Last 30 lines of log ({log_path}):\n{tail}")
        except Exception:
            pass
        raise RuntimeError(f"{stage_name} failed with exit code {res.returncode}.")
    timer.record(
        stage_name,
        "ok",
        duration,
        ";".join(filter(None, (cmd_str, resource_details))),
    )
    log_info(f"Finished {stage_name} (Duration: {format_duration(duration)})")
