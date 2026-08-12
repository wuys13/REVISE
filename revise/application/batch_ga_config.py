"""Minimal YAML compilation for the Batch GA entrypoint."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class BatchGAConfig:
    spatial_path: Path
    reference_path: Path
    filter_column: str
    filter_value: str
    spatial_preprocessing: dict
    reference_preprocessing: dict
    broad_column: str
    output_dir: Path
    seed: int


def _relative_path(root: Path, value: object, field: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{field} must be a non-empty path")
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"{field} must be relative to paths.root_dir")
    return root / path


def load_batch_config(config_path: str | Path) -> BatchGAConfig:
    source = Path(config_path)
    document = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    paths = document["paths"]
    root_value = paths["root_dir"]
    root = source.parent if root_value == "." else Path(root_value)
    if not root.is_absolute():
        root = (source.parent / root).resolve()
    inputs = document["inputs"]
    spatial = inputs["st"]
    reference = inputs["reference"]
    if spatial.get("format") != "h5ad" or reference.get("format") != "h5ad":
        raise ValueError("batch GA only supports h5ad inputs")
    preprocessing = document["preprocessing"]
    return BatchGAConfig(
        spatial_path=_relative_path(root, spatial["path"], "inputs.st.path"),
        reference_path=_relative_path(root, reference["path"], "inputs.reference.path"),
        filter_column=reference["filter_column"],
        filter_value=reference["filter_value"],
        spatial_preprocessing=dict(preprocessing["spatial"]),
        reference_preprocessing=dict(preprocessing["reference"]),
        broad_column=document["global_anchoring"]["broad_column"],
        output_dir=_relative_path(root, document["output"]["dir"], "output.dir"),
        seed=document["execution"]["seed"],
    )


__all__ = ["BatchGAConfig", "load_batch_config"]
