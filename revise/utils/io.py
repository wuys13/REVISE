from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Mapping
from uuid import uuid4


def benchmark_case_leaf(route_key: str, io_cfg: Mapping[str, Any]) -> str | None:
    if route_key == "sim2real:segmentation":
        seg_method = io_cfg.get("seg_method")
        return str(seg_method) if seg_method else None

    if route_key == "sim2real:bin2cell":
        seg_method = io_cfg.get("seg_method")
        return str(seg_method) if seg_method else "bin2cell"

    if route_key == "sim2real:batch_effect":
        batch_by_sc_ref = {
            "selected_xenium.h5ad": "1",
            "real_sc_ref_part.h5ad": "2",
            "real_sc_ref_all.h5ad": "3",
            "real_sc_ref_others.h5ad": "4",
        }
        sc_ref_file = io_cfg.get("sc_ref_file")
        batch_id = batch_by_sc_ref.get(str(sc_ref_file)) if sc_ref_file else None
        spot_size = io_cfg.get("spot_size")
        if batch_id is None:
            return None
        if spot_size is None:
            return batch_id
        return f"{spot_size}_{batch_id}"

    if route_key == "sim2real:spot_size":
        spot_size = io_cfg.get("spot_size")
        return str(spot_size) if spot_size is not None else None

    if route_key in {"sim2real:gene_panel", "sim2real:gene_dropout"}:
        return None

    return None


def build_task_dir(output_root: str, sample_name: str, route_key: str, io_cfg: Mapping[str, Any] | None = None) -> Path:
    if route_key.startswith("sim2real:"):
        return Path(output_root) / sample_name
    return build_run_dir(output_root=output_root, sample_name=sample_name, route_key=route_key, io_cfg=io_cfg)


def build_run_dir(output_root: str, sample_name: str, route_key: str, io_cfg: Mapping[str, Any] | None = None) -> Path:
    if route_key.startswith("sim2real:"):
        task_dir = Path(output_root) / sample_name
        leaf = benchmark_case_leaf(route_key, io_cfg or {})
        return task_dir if leaf is None else task_dir / leaf

    safe_route = route_key.replace(":", "__")
    leaf = f"{datetime.now():%Y%m%d_%H%M%S_%f}_{uuid4().hex[:8]}"
    return Path(output_root) / sample_name / safe_route / leaf
