from __future__ import annotations

import copy
import os
import random
from typing import Any
from typing import Optional

import numpy as np


_LOCATION_ONLY_IO_KEYS = {
    "data_root",
    "output_root",
    "st_path",
    "sc_ref_path",
    "st_file",
    "sc_ref_file",
    "gt_svc_file",
    "spatialdata_path",
}


def canonical_config_projection(config: dict[str, Any]) -> dict[str, Any]:
    """Return the resolved config without input and output locators."""
    projected = copy.deepcopy(config)
    projected["io"] = {
        key: value
        for key, value in projected.get("io", {}).items()
        if key not in _LOCATION_ONLY_IO_KEYS
    }
    return projected


def set_global_seed(seed: Optional[int], deterministic: bool = True) -> None:
    if seed is None:
        return
    random.seed(seed)
    np.random.seed(seed)
    if deterministic:
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":16:8")
