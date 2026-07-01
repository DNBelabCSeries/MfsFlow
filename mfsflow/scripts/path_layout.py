"""Bootstrap access to mfsflow.path_layout for direct script execution."""

import os
import sys


_PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, _PACKAGE_ROOT)

from mfsflow.path_layout import (  # noqa: E402,F401
    barcode_dir,
    config_dir,
    ensure_layout,
    expression_dir,
    intermediate_dir,
    load_config,
    logs_dir,
    outputs_dir,
    stage_state_dir,
    stats_dir,
    tmp_merge_dir,
)
