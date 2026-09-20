"""Assemble the zero-gene H2 baseline from historical spatial carriers.

Only obs metadata and spatial coordinates are read from backed H5AD handles.
Source expression matrices are never accessed or concatenated.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from revise.analysis.assembly_comparison import file_sha256, normalize_broad_labels  # noqa: E402


CARRIER_TYPES = ("T", "Mono_Macro", "Fibroblast")
DEFAULT_CARRIER_PATHS = {
    "T": ROOT / "results/sc_SVC_case/P2CRC_Xenium/T/spatial.h5ad",
    "Mono_Macro": ROOT / "results/sc_SVC_case/P2CRC_Xenium/Mono_Macro/spatial.h5ad",
    "Fibroblast": ROOT / "results/sc_SVC_case/P2CRC_Xenium/Fibroblast/spatial.h5ad",
}


def _validate_ids(ids: pd.Index, *, source: str) -> None:
    if not ids.is_unique:
        raise ValueError(f"{source}: obs_names must be unique")
    blank = np.asarray(ids.isna()) | (np.asarray(ids.astype(str).str.strip()) == "")
    if bool(np.any(blank)):
        raise ValueError(f"{source}: obs_names must not be empty")


def _close_backed(adata: AnnData) -> None:
    file = getattr(adata, "file", None)
    close = getattr(file, "close", None)
    if close is not None:
        close()


def _read_carrier_metadata(
    path: str | Path,
    *,
    broad_column: str = "Level1",
    subtype_column: str = "SVC_cluster",
    spatial_key: str = "spatial",
    source_name: str,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    """Read only obs and spatial metadata from one backed H5AD carrier."""

    source_path = Path(path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"{source_name}: missing carrier {source_path}")

    native = read_h5ad(source_path, backed="r")
    try:
        if broad_column not in native.obs:
            raise KeyError(f"{source_name}: missing obs[{broad_column!r}]")
        if subtype_column not in native.obs:
            raise KeyError(f"{source_name}: missing obs[{subtype_column!r}]")
        if spatial_key not in native.obsm:
            raise KeyError(f"{source_name}: missing obsm[{spatial_key!r}]")

        ids = pd.Index(native.obs_names)
        _validate_ids(ids, source=source_name)
        obs = native.obs.loc[:, [broad_column, subtype_column]].copy()
        obs.index = ids
        obs[broad_column] = normalize_broad_labels(obs[broad_column])
        coordinates = np.asarray(native.obsm[spatial_key]).copy()
        shape = tuple(coordinates.shape)
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError(
                f"{source_name}: obsm[{spatial_key!r}] must be 2D with at least two columns"
            )
        if coordinates.shape[0] != len(obs):
            raise ValueError(
                f"{source_name}: spatial rows {coordinates.shape[0]} != obs rows {len(obs)}"
            )
        if not np.issubdtype(coordinates.dtype, np.number):
            raise TypeError(f"{source_name}: spatial coordinates must be numeric")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError(f"{source_name}: spatial coordinates must be finite")
    finally:
        _close_backed(native)

    metadata = {
        "source_name": source_name,
        "path": str(source_path),
        "sha256": file_sha256(source_path),
        "n_obs": int(len(obs)),
        "n_vars_read": 0,
        "broad_column": broad_column,
        "subtype_column": subtype_column,
        "spatial_key": spatial_key,
        "spatial_shape": list(shape),
        "expression_loaded": False,
    }
    return obs, coordinates, metadata


def assemble_baseline(
    carrier_paths: Mapping[str, str | Path],
    output_path: str | Path,
    *,
    broad_column: str = "Level1",
    subtype_column: str = "SVC_cluster",
    spatial_key: str = "spatial",
    coordinate_unit: str = "pixel",
    microns_per_coordinate: float | None = 0.2125,
) -> dict[str, Any]:
    """Create a zero-gene baseline from three explicit carrier paths."""

    expected = set(CARRIER_TYPES)
    supplied = set(carrier_paths)
    if supplied != expected:
        raise ValueError(
            f"carrier_paths must contain exactly {CARRIER_TYPES}; "
            f"missing={sorted(expected - supplied)}, extra={sorted(supplied - expected)}"
        )

    output = Path(output_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    observations: list[pd.DataFrame] = []
    coordinates: list[np.ndarray] = []
    source_records: list[dict[str, Any]] = []
    seen: set[Any] = set()

    for cell_type in CARRIER_TYPES:
        obs, xy, record = _read_carrier_metadata(
            carrier_paths[cell_type],
            broad_column=broad_column,
            subtype_column=subtype_column,
            spatial_key=spatial_key,
            source_name=cell_type,
        )
        duplicate = sorted(set(obs.index).intersection(seen))
        if duplicate:
            raise ValueError(
                f"cross-carrier duplicate obs_names for {cell_type}: {duplicate[:5]}"
            )
        seen.update(obs.index)
        observations.append(obs)
        coordinates.append(xy)
        record["cell_type"] = cell_type
        source_records.append(record)

    combined_obs = pd.concat(observations, axis=0)
    combined_xy = np.concatenate(coordinates, axis=0)
    baseline = AnnData(
        X=np.empty((len(combined_obs), 0), dtype=np.float32),
        obs=combined_obs,
        obsm={spatial_key: combined_xy},
    )
    source_sha256 = {
        record["cell_type"]: record["sha256"] for record in source_records
    }
    source_records_by_type = {
        record["cell_type"]: record for record in source_records
    }
    baseline.uns["assembly_baseline"] = {
        "schema": "h2-historical-carrier-obs-spatial-v1",
        "zero_gene_axis": True,
        "expression_loaded": False,
        "broad_column": broad_column,
        "subtype_column": subtype_column,
        "spatial_key": spatial_key,
        "label_normalization": "replace slash with underscore in Level1 only; preserve NA",
        "coordinate_provenance": {
            "source": "formal sample.yaml and baseline provenance",
            "unit": coordinate_unit,
            "microns_per_coordinate": microns_per_coordinate,
        },
        "source_sha256": source_sha256,
        "sources": source_records_by_type,
        "n_obs": int(baseline.n_obs),
        "n_vars": int(baseline.n_vars),
    }
    baseline.write_h5ad(output)

    return {
        "output": str(output),
        "n_obs": int(baseline.n_obs),
        "n_vars": int(baseline.n_vars),
        "source_sha256": source_sha256,
        "sources": source_records_by_type,
        "coordinate_provenance": baseline.uns["assembly_baseline"][
            "coordinate_provenance"
        ],
    }


def _parse_carrier_specs(values: Sequence[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"--carrier must be TYPE=PATH, got {value!r}")
        cell_type, raw_path = value.split("=", 1)
        cell_type = cell_type.strip()
        if cell_type not in CARRIER_TYPES:
            raise ValueError(f"unknown carrier type {cell_type!r}")
        if cell_type in parsed:
            raise ValueError(f"duplicate carrier type {cell_type!r}")
        parsed[cell_type] = Path(raw_path).expanduser()
    return parsed


def prepare_baseline(
    carrier_paths: Mapping[str, str | Path],
    output_path: str | Path,
    **kwargs: Any,
) -> dict[str, Any]:
    """Named entry point for callers that treat H2 as baseline preparation."""

    return assemble_baseline(carrier_paths, output_path, **kwargs)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--carrier",
        action="append",
        default=[],
        metavar="TYPE=PATH",
        help="repeat for T, Mono_Macro, and Fibroblast",
    )
    parser.add_argument("--T", dest="t_carrier", type=Path)
    parser.add_argument("--t-carrier", dest="t_carrier_alias", type=Path)
    parser.add_argument("--mono-macro-carrier", type=Path)
    parser.add_argument("--fibroblast-carrier", type=Path)
    parser.add_argument("--coordinate-unit", default="pixel")
    parser.add_argument("--microns-per-coordinate", type=float, default=0.2125)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    carrier_paths = _parse_carrier_specs(args.carrier)
    explicit = {
        "T": args.t_carrier or args.t_carrier_alias,
        "Mono_Macro": args.mono_macro_carrier,
        "Fibroblast": args.fibroblast_carrier,
    }
    for cell_type, path in explicit.items():
        if path is not None:
            if cell_type in carrier_paths:
                raise ValueError(f"carrier {cell_type!r} supplied more than once")
            carrier_paths[cell_type] = path
    if not carrier_paths:
        carrier_paths = dict(DEFAULT_CARRIER_PATHS)
    summary = assemble_baseline(
        carrier_paths,
        args.output,
        coordinate_unit=args.coordinate_unit,
        microns_per_coordinate=args.microns_per_coordinate,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
