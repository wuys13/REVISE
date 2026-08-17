"""Application YAML loading and compilation.

This module deliberately stops at a validated, resolved application
configuration. Runtime mapping is package-owned; the source entrypoint keeps
the user-facing load, preprocess, reconstruct, and publication flow visible.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import importlib.resources as resources
import math
from pathlib import Path
from typing import Any, Mapping
import unicodedata

import yaml

from revise.config.authority import ENGINE_DEFAULTS
from revise.utils.labels import normalize_cell_type_label


class ApplicationConfigError(ValueError):
    """Raised when an application YAML cannot be used."""


@dataclass(frozen=True)
class ApplicationConfigSource:
    label: str
    payload: bytes
    path: Path | None


@dataclass(frozen=True)
class ApplicationConfig:
    config_path: Path | None
    config_sha256: str
    config_source: str
    declared_root: str
    resolved_root: Path
    cwd: Path
    svc_type: str
    mode: str | None
    st_path: Path
    reference_path: Path
    st_format: str
    spatialdata_table: str | None
    spatialdata_element: str | None
    reference_filter_column: str | None
    reference_filter_value: str | None
    spatial_min_transcript_counts: int | None
    spatial_min_counts: int | None
    spatial_min_cell_counts: int
    reference_min_transcript_counts: int | None
    reference_min_genes: int | None
    reference_min_cell_counts: int
    broad_column: str
    subtype_column: str | None
    select_cell_type: str | None
    local_refinement_strength: float | None
    local_refinement_alpha: float | None
    local_refinement_resolutions: tuple[float, ...] | None
    local_refinement_graph_method: str | None
    local_refinement_graph_alpha: float | None
    local_refinement_graph_n_neighbors: int | None
    local_refinement_graph_exp_neighbors: int | None
    local_refinement_graph_spatial_neighbors: int | None
    local_refinement_match_spot_sum: bool | None
    ot_method: str | None
    pm_on_cell_path: Path | None
    output_root: Path
    output_dir: Path
    output_name: str | None
    seed: int | None

    @property
    def source_path(self) -> str:
        return self.config_source

    @property
    def resolved_inputs(self) -> dict[str, str]:
        paths = {
            "st": str(self.st_path),
            "reference": str(self.reference_path),
            "output_dir": str(self.output_dir),
        }
        if self.pm_on_cell_path is not None:
            paths["pm_on_cell"] = str(self.pm_on_cell_path)
        return paths


def _compile_engine_config(
    config: ApplicationConfig,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Translate a compiled Application YAML into engine configuration sections."""
    runtime = {"seed": config.seed} if config.seed is not None else {}
    io = {
        "st_path": str(config.st_path),
        "sc_ref_path": str(config.reference_path),
        "pm_on_cell_path": str(config.pm_on_cell_path) if config.pm_on_cell_path else "",
        "output_root": str(config.output_dir),
        "sample_name": config.svc_type,
        "patient_key": "",
        "save_outputs": False,
        "input_format": config.st_format,
        "data_root": "",
        "st_file": "",
        "sc_ref_file": "",
    }
    if config.st_format in {"spatialdata", "auto"}:
        io["spatialdata_path"] = str(config.st_path)
        if config.spatialdata_table is not None:
            io["spatialdata_table"] = config.spatialdata_table
        if config.spatialdata_element is not None:
            io["spatialdata_spatial_element"] = config.spatialdata_element

    algorithm: dict[str, Any] = {"columns": {"cell_type_col": config.broad_column}}
    if config.subtype_column is not None:
        algorithm["columns"]["sub_cell_type_col"] = config.subtype_column
    if config.select_cell_type is not None:
        algorithm["sc"] = {"select_ct": config.select_cell_type}
    if config.local_refinement_strength is not None:
        algorithm["local_refinement"] = {"strength": config.local_refinement_strength}
    if config.local_refinement_alpha is not None:
        algorithm["graph"] = {"alpha": config.local_refinement_alpha}
    if config.local_refinement_resolutions is not None:
        algorithm.setdefault("sc", {})["resolutions"] = list(config.local_refinement_resolutions)
    if config.local_refinement_graph_method is not None:
        algorithm["graph"] = {
            "method": config.local_refinement_graph_method,
            "alpha": config.local_refinement_graph_alpha,
            "n_neighbors": config.local_refinement_graph_n_neighbors,
            "exp_neighbors": config.local_refinement_graph_exp_neighbors,
            "spatial_neighbors": config.local_refinement_graph_spatial_neighbors,
        }
    if config.local_refinement_match_spot_sum is not None:
        algorithm.setdefault("sc", {})["match_spot_sum"] = config.local_refinement_match_spot_sum
    if config.ot_method is not None:
        algorithm["ot"] = {
            "ga": {"solver": config.ot_method},
            "lr": {"solver": config.ot_method},
        }
    return runtime, io, algorithm


_OFFICIAL_TEMPLATES = frozenset(
    {
        "VisiumHD.yaml",
        "Xenium.yaml",
        "Visium.yaml",
    }
)
_TOP_LEVEL_KEYS = {
    "schema_version",
    "application",
    "paths",
    "algorithm",
    "inputs",
    "preprocessing",
    "global_anchoring",
    "local_refinement",
    "output",
    "execution",
}
_SVC_TYPES = {"sp-SVC", "sc-SVC"}
_SC_SVC_MODES = {"cluster", "sr"}
_ST_FORMATS = {"h5ad", "spatialdata", "auto"}
_ALL_CELL_TYPES = {"", "all", "*", "__all__", "all_cell_types"}
_MIGRATION_MESSAGES = {
    "application.sample_name": "use output.name",
    "inputs.mode": "Application now accepts exact inputs.*.path fields",
    "inputs.data_root": "use inputs.st.path and inputs.reference.path",
    "inputs.st.file": "use inputs.st.path",
    "inputs.reference.file": "use inputs.reference.path",
    "inputs.reference.patient_key": "prepare the reference before running",
    "output.path": "use output.dir and output.name",
}


def _mapping(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ApplicationConfigError(f"{field} must be a mapping")
    return dict(value)


def _reject_unknown(mapping: Mapping[str, Any], allowed: set[str], field: str) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        first = f"{field}.{unknown[0]}"
        migration = _MIGRATION_MESSAGES.get(first)
        if migration:
            raise ApplicationConfigError(f"{first} was removed; {migration}")
        raise ApplicationConfigError(
            f"{field} contains unknown field(s): {', '.join(unknown)}"
        )


def _required(mapping: Mapping[str, Any], key: str, field: str) -> Any:
    if key not in mapping:
        raise ApplicationConfigError(f"{field}.{key} is required")
    return mapping[key]


def _string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ApplicationConfigError(f"{field} must be a non-empty string")
    return value.strip()


def _optional_string(value: Any, field: str) -> str | None:
    if value is None:
        return None
    return _string(value, field)


def _concrete_cell_type(value: Any, field: str) -> str:
    if not isinstance(value, str):
        raise ApplicationConfigError(
            f"{field} must name one concrete broad cell type safe for an output directory"
        )
    raw_selected = value.strip()
    selected = normalize_cell_type_label(value)
    if (
        not raw_selected
        or raw_selected.lower() in _ALL_CELL_TYPES
        or "*" in raw_selected
        or raw_selected in {".", ".."}
        or any(part in {".", ".."} for part in raw_selected.split("/"))
        or "\\" in raw_selected
        or any(unicodedata.category(character) == "Cc" for character in value)
    ):
        raise ApplicationConfigError(
            f"{field} must name one concrete broad cell type safe for an output directory"
        )
    return selected


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ApplicationConfigError(f"{field} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ApplicationConfigError(f"{field} must be a finite non-negative number")
    return result


def _count(value: Any, field: str, *, optional: bool = False) -> int | None:
    if value is None and optional:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ApplicationConfigError(f"{field} must be a non-negative integer")
    return value


def _root_dir(value: Any, *, cwd: Path) -> tuple[str, Path]:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ApplicationConfigError(
            "paths.root_dir must be literal . or an existing absolute directory"
        )
    if value == ".":
        return value, cwd
    path = Path(value)
    if not path.is_absolute() or not path.exists() or not path.is_dir():
        raise ApplicationConfigError(
            "paths.root_dir must be literal . or an existing absolute directory"
        )
    return value, path.resolve()


def _relative_child(root: Path, value: Any, field: str) -> Path:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ApplicationConfigError(
            f"{field} must be a non-empty relative path under paths.root_dir"
        )
    child = Path(value)
    if child.is_absolute() or ".." in child.parts:
        raise ApplicationConfigError(f"{field} must be relative to paths.root_dir")
    return root / child


def _parse_st(inputs: Mapping[str, Any], *, root: Path):
    st = _mapping(_required(inputs, "st", "inputs"), "inputs.st")
    _reject_unknown(st, {"path", "format", "spatialdata"}, "inputs.st")
    path = _relative_child(root, _required(st, "path", "inputs.st"), "inputs.st.path")
    st_format = _string(_required(st, "format", "inputs.st"), "inputs.st.format").lower()
    if st_format not in _ST_FORMATS:
        raise ApplicationConfigError("inputs.st.format must be one of: h5ad, spatialdata, auto")
    table = element = None
    if "spatialdata" in st:
        if st_format == "h5ad":
            raise ApplicationConfigError(
                "inputs.st.spatialdata is not allowed when inputs.st.format is h5ad"
            )
        spatialdata = _mapping(st["spatialdata"], "inputs.st.spatialdata")
        _reject_unknown(spatialdata, {"table", "element"}, "inputs.st.spatialdata")
        table = _optional_string(spatialdata.get("table"), "inputs.st.spatialdata.table")
        element = _optional_string(spatialdata.get("element"), "inputs.st.spatialdata.element")
    return path, st_format, table, element


def _parse_reference(
    inputs: Mapping[str, Any],
    *,
    root: Path,
) -> tuple[Path, str | None, str | None]:
    reference = _mapping(_required(inputs, "reference", "inputs"), "inputs.reference")
    _reject_unknown(
        reference,
        {"path", "format", "filter_column", "filter_value"},
        "inputs.reference",
    )
    path = _relative_child(root, _required(reference, "path", "inputs.reference"), "inputs.reference.path")
    if _string(_required(reference, "format", "inputs.reference"), "inputs.reference.format").lower() != "h5ad":
        raise ApplicationConfigError("inputs.reference.format must be h5ad")
    filter_column = _optional_string(
        reference.get("filter_column"),
        "inputs.reference.filter_column",
    )
    filter_value = _optional_string(
        reference.get("filter_value"),
        "inputs.reference.filter_value",
    )
    if (filter_column is None) != (filter_value is None):
        raise ApplicationConfigError(
            "inputs.reference.filter_column and filter_value must be supplied together"
        )
    return path, filter_column, filter_value


def _canonical_template_name(value: str) -> str | None:
    path = Path(value)
    if path.is_absolute():
        return None
    parts = path.parts
    if len(parts) == 1 and parts[0] in _OFFICIAL_TEMPLATES:
        return parts[0]
    if len(parts) == 3 and parts[:2] == ("configs", "application"):
        if parts[2] in _OFFICIAL_TEMPLATES:
            return parts[2]
    return None


def _local_source(path: Path) -> ApplicationConfigSource:
    resolved = path.expanduser().resolve()
    try:
        payload = resolved.read_bytes()
    except OSError as exc:
        raise ApplicationConfigError(f"Cannot read application config {resolved}: {exc}") from exc
    return ApplicationConfigSource(str(resolved), payload, resolved)


def resolve_application_source(config: str | Path) -> ApplicationConfigSource:
    """Use an existing file first, then only the three official templates."""
    value = str(config)
    path = Path(value).expanduser()
    if path.exists():
        if not path.is_file():
            raise ApplicationConfigError(f"Cannot read application config {path}: not a file")
        return _local_source(path)
    template_name = _canonical_template_name(value)
    if template_name is None:
        raise ApplicationConfigError(f"Cannot read application config {path}: file does not exist")
    if Path(value).parts == (template_name,):
        mirror = Path.cwd() / "configs" / "application" / template_name
        if mirror.is_file():
            return _local_source(mirror)
    try:
        payload = resources.files("revise.application").joinpath("templates", template_name).read_bytes()
    except (FileNotFoundError, ModuleNotFoundError, OSError) as exc:
        raise ApplicationConfigError(
            f"Cannot read packaged application template {template_name}: {exc}"
        ) from exc
    return ApplicationConfigSource(
        f"package:revise.application.templates/{template_name}", payload, None
    )


def load_application_yaml(config: str | Path) -> tuple[ApplicationConfigSource, dict[str, Any]]:
    source = resolve_application_source(config)
    try:
        raw = yaml.safe_load(source.payload)
    except yaml.YAMLError as exc:
        raise ApplicationConfigError(f"Invalid application YAML {source.label}: {exc}") from exc
    document = _mapping(raw, "application config")
    if "schema_version" not in document and any(
        key in document for key in ("version", "defaults", "router", "profiles")
    ):
        raise ApplicationConfigError(
            "--config requires an application YAML; the engine config is internally managed"
        )
    return source, document


def compile_application_config(
    document: Mapping[str, Any],
    *,
    source: ApplicationConfigSource,
    cwd: Path | None = None,
) -> ApplicationConfig:
    """Validate one Application YAML document."""
    cwd = (cwd or Path.cwd()).resolve()
    document = _mapping(document, "application config")
    _reject_unknown(document, _TOP_LEVEL_KEYS, "application config")
    missing = sorted(_TOP_LEVEL_KEYS - set(document))
    if missing:
        raise ApplicationConfigError(
            f"application config is missing required field(s): {', '.join(missing)}"
        )
    if document["schema_version"] != 1:
        raise ApplicationConfigError("schema_version must be 1")

    application = _mapping(document["application"], "application")
    _reject_unknown(application, {"svc_type", "mode"}, "application")
    svc_type = _string(_required(application, "svc_type", "application"), "application.svc_type")
    if svc_type == "sc-SVC-sr":
        raise ApplicationConfigError(
            "application.svc_type sc-SVC-sr was removed; use sc-SVC with application.mode: sr"
        )
    if svc_type not in _SVC_TYPES:
        raise ApplicationConfigError("application.svc_type must be one of: sp-SVC, sc-SVC")
    mode = _optional_string(application.get("mode"), "application.mode")
    if svc_type == "sc-SVC":
        if mode is None:
            raise ApplicationConfigError(
                "application.mode is required for sc-SVC; use cluster or sr"
            )
        if mode not in _SC_SVC_MODES:
            raise ApplicationConfigError(
                "application.mode must be one of: cluster, sr"
            )
    elif mode is not None:
        raise ApplicationConfigError("application.mode is only valid for sc-SVC")

    paths = _mapping(document["paths"], "paths")
    _reject_unknown(paths, {"root_dir"}, "paths")
    declared_root, resolved_root = _root_dir(_required(paths, "root_dir", "paths"), cwd=cwd)

    algorithm = _mapping(document["algorithm"], "algorithm")
    if "base_config" in algorithm:
        raise ApplicationConfigError("algorithm.base_config is not public; the engine config is internally managed")
    _reject_unknown(algorithm, {"ot_method"}, "algorithm")
    ot_method = _optional_string(algorithm.get("ot_method"), "algorithm.ot_method")
    if ot_method is not None and ot_method not in {"pot", "tacco"}:
        raise ApplicationConfigError("algorithm.ot_method must be pot or tacco")

    inputs = _mapping(document["inputs"], "inputs")
    _reject_unknown(inputs, {"st", "reference", "pm_on_cell"}, "inputs")
    st_path, st_format, table, element = _parse_st(inputs, root=resolved_root)
    reference_path, reference_filter_column, reference_filter_value = _parse_reference(
        inputs,
        root=resolved_root,
    )
    pm_path = None
    if "pm_on_cell" in inputs:
        if mode != "sr":
            raise ApplicationConfigError("inputs.pm_on_cell is only allowed for sc-SVC sr mode")
        pm = _mapping(inputs["pm_on_cell"], "inputs.pm_on_cell")
        _reject_unknown(pm, {"path"}, "inputs.pm_on_cell")
        pm_path = _relative_child(resolved_root, _required(pm, "path", "inputs.pm_on_cell"), "inputs.pm_on_cell.path")

    preprocessing = _mapping(document["preprocessing"], "preprocessing")
    _reject_unknown(preprocessing, {"spatial", "reference"}, "preprocessing")
    spatial_preprocessing = _mapping(
        _required(preprocessing, "spatial", "preprocessing"),
        "preprocessing.spatial",
    )
    _reject_unknown(
        spatial_preprocessing,
        {"min_transcript_counts", "min_counts", "min_cell_counts"},
        "preprocessing.spatial",
    )
    spatial_min_transcript_counts = _count(
        _required(
            spatial_preprocessing,
            "min_transcript_counts",
            "preprocessing.spatial",
        ),
        "preprocessing.spatial.min_transcript_counts",
        optional=True,
    )
    spatial_min_cell_counts = _count(
        _required(spatial_preprocessing, "min_cell_counts", "preprocessing.spatial"),
        "preprocessing.spatial.min_cell_counts",
    )
    spatial_min_counts = _count(
        spatial_preprocessing.get("min_counts"),
        "preprocessing.spatial.min_counts",
        optional=True,
    )
    reference_preprocessing = _mapping(
        _required(preprocessing, "reference", "preprocessing"),
        "preprocessing.reference",
    )
    _reject_unknown(
        reference_preprocessing,
        {"min_transcript_counts", "min_genes", "min_cell_counts"},
        "preprocessing.reference",
    )
    reference_min_transcript_counts = _count(
        _required(
            reference_preprocessing,
            "min_transcript_counts",
            "preprocessing.reference",
        ),
        "preprocessing.reference.min_transcript_counts",
        optional=True,
    )
    reference_min_cell_counts = _count(
        _required(
            reference_preprocessing,
            "min_cell_counts",
            "preprocessing.reference",
        ),
        "preprocessing.reference.min_cell_counts",
    )
    reference_min_genes = _count(
        reference_preprocessing.get("min_genes"),
        "preprocessing.reference.min_genes",
        optional=True,
    )

    anchoring = _mapping(document["global_anchoring"], "global_anchoring")
    _reject_unknown(anchoring, {"broad_column"}, "global_anchoring")
    broad_column = _string(_required(anchoring, "broad_column", "global_anchoring"), "global_anchoring.broad_column")

    refinement = _mapping(document["local_refinement"], "local_refinement")
    subtype = select = strength = alpha = resolutions = None
    graph_method = graph_alpha = graph_n_neighbors = None
    graph_exp_neighbors = graph_spatial_neighbors = match_spot_sum = None
    if mode == "cluster":
        _reject_unknown(
            refinement,
            {"subtype_column", "select_cell_type", "alpha", "resolutions"},
            "local_refinement",
        )
        subtype = _string(_required(refinement, "subtype_column", "local_refinement"), "local_refinement.subtype_column")
        select = _concrete_cell_type(
            _required(refinement, "select_cell_type", "local_refinement"),
            "local_refinement.select_cell_type",
        )
        alpha = _number(
            _required(refinement, "alpha", "local_refinement"),
            "local_refinement.alpha",
        )
        raw_resolutions = _required(refinement, "resolutions", "local_refinement")
        if not isinstance(raw_resolutions, list) or not raw_resolutions:
            raise ApplicationConfigError("local_refinement.resolutions must be a non-empty list")
        resolutions = tuple(
            _number(value, "local_refinement.resolutions")
            for value in raw_resolutions
        )
    elif mode == "sr":
        _reject_unknown(
            refinement,
            {"strength", "graph", "match_spot_sum"},
            "local_refinement",
        )
        if "strength" in refinement:
            strength = _number(refinement["strength"], "local_refinement.strength")
        if "graph" in refinement:
            graph = _mapping(refinement["graph"], "local_refinement.graph")
            _reject_unknown(
                graph,
                {"method", "alpha", "n_neighbors", "exp_neighbors", "spatial_neighbors"},
                "local_refinement.graph",
            )
            graph_method = _string(
                _required(graph, "method", "local_refinement.graph"),
                "local_refinement.graph.method",
            )
            graph_alpha = _number(
                _required(graph, "alpha", "local_refinement.graph"),
                "local_refinement.graph.alpha",
            )
            graph_n_neighbors = _count(
                _required(graph, "n_neighbors", "local_refinement.graph"),
                "local_refinement.graph.n_neighbors",
            )
            graph_exp_neighbors = _count(
                _required(graph, "exp_neighbors", "local_refinement.graph"),
                "local_refinement.graph.exp_neighbors",
            )
            graph_spatial_neighbors = _count(
                _required(graph, "spatial_neighbors", "local_refinement.graph"),
                "local_refinement.graph.spatial_neighbors",
            )
        if "match_spot_sum" in refinement:
            match_spot_sum = refinement["match_spot_sum"]
            if not isinstance(match_spot_sum, bool):
                raise ApplicationConfigError(
                    "local_refinement.match_spot_sum must be a boolean"
                )
    else:
        _reject_unknown(refinement, {"strength"}, "local_refinement")
        if "strength" in refinement:
            strength = _number(refinement["strength"], "local_refinement.strength")

    output = _mapping(document["output"], "output")
    _reject_unknown(output, {"dir", "name"}, "output")
    output_root = _relative_child(resolved_root, _required(output, "dir", "output"), "output.dir")
    output_dir = output_root / select if mode == "cluster" else output_root
    output_name = _optional_string(output.get("name"), "output.name")
    if output_name is not None:
        if output_name in {".", ".."} or "/" in output_name or "\\" in output_name:
            raise ApplicationConfigError("output.name must be a filename stem without separators")
        if output_name.endswith(".h5ad"):
            raise ApplicationConfigError("output.name must not include .h5ad")

    execution = _mapping(document["execution"], "execution")
    if "action" in execution:
        raise ApplicationConfigError("execution.action is not supported")
    _reject_unknown(execution, {"seed"}, "execution")
    seed = execution.get("seed", ENGINE_DEFAULTS["runtime"]["seed"])
    if seed is None:
        seed = ENGINE_DEFAULTS["runtime"]["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ApplicationConfigError("execution.seed must be an integer")
    if not 0 <= seed <= 2**32 - 1:
        raise ApplicationConfigError("execution.seed must be between 0 and 4294967295")

    return ApplicationConfig(
        config_path=source.path,
        config_sha256=sha256(source.payload).hexdigest(),
        config_source=source.label,
        declared_root=declared_root,
        resolved_root=resolved_root,
        cwd=cwd,
        svc_type=svc_type,
        mode=mode,
        st_path=st_path,
        reference_path=reference_path,
        st_format=st_format,
        spatialdata_table=table,
        spatialdata_element=element,
        reference_filter_column=reference_filter_column,
        reference_filter_value=reference_filter_value,
        spatial_min_transcript_counts=spatial_min_transcript_counts,
        spatial_min_counts=spatial_min_counts,
        spatial_min_cell_counts=spatial_min_cell_counts,
        reference_min_transcript_counts=reference_min_transcript_counts,
        reference_min_genes=reference_min_genes,
        reference_min_cell_counts=reference_min_cell_counts,
        broad_column=broad_column,
        subtype_column=subtype,
        select_cell_type=select,
        local_refinement_strength=strength,
        local_refinement_alpha=alpha,
        local_refinement_resolutions=resolutions,
        local_refinement_graph_method=graph_method,
        local_refinement_graph_alpha=graph_alpha,
        local_refinement_graph_n_neighbors=graph_n_neighbors,
        local_refinement_graph_exp_neighbors=graph_exp_neighbors,
        local_refinement_graph_spatial_neighbors=graph_spatial_neighbors,
        local_refinement_match_spot_sum=match_spot_sum,
        ot_method=ot_method,
        pm_on_cell_path=pm_path,
        output_root=output_root,
        output_dir=output_dir,
        output_name=output_name,
        seed=seed,
    )


def override_select_cell_type(
    config: ApplicationConfig,
    select_ct: str | None,
) -> ApplicationConfig:
    """Apply a cluster override and derive its output from the configured root."""
    if select_ct is None:
        return config
    if config.svc_type != "sc-SVC" or config.mode != "cluster":
        raise ApplicationConfigError(
            "--select-ct is only valid for sc-SVC cluster mode"
        )
    selected = _concrete_cell_type(select_ct, "--select-ct")
    return replace(
        config,
        select_cell_type=selected,
        output_dir=config.output_root / selected,
    )


__all__ = [
    "ApplicationConfig",
    "ApplicationConfigError",
    "ApplicationConfigSource",
    "compile_application_config",
    "load_application_yaml",
    "override_select_cell_type",
    "resolve_application_source",
]
