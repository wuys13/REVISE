"""Application YAML loading and compilation.

This module deliberately stops at a validated, resolved application
configuration.  Runtime mapping and result publication live in the source
entrypoint so the complete user-facing flow remains visible there.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import importlib.resources as resources
import math
from pathlib import Path
from typing import Any, Mapping

import yaml


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
    st_path: Path
    reference_path: Path
    st_format: str
    spatialdata_table: str | None
    spatialdata_element: str | None
    broad_column: str
    subtype_column: str | None
    select_cell_type: str | None
    local_refinement_strength: float | None
    ot_method: str | None
    pm_on_cell_path: Path | None
    output_dir: Path
    output_name: str
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


_OFFICIAL_TEMPLATES = frozenset(
    {
        "Xenium_T.yaml",
        "Xenium_Fib.yaml",
        "Xenium_Mono.yaml",
        "VisiumHD.yaml",
        "Visium.yaml",
    }
)
_TOP_LEVEL_KEYS = {
    "schema_version",
    "application",
    "paths",
    "algorithm",
    "inputs",
    "global_anchoring",
    "local_refinement",
    "output",
    "execution",
}
_SVC_TYPES = {"sp-SVC", "sc-SVC", "sc-SVC-sr"}
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


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ApplicationConfigError(f"{field} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ApplicationConfigError(f"{field} must be a finite non-negative number")
    return result


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


def _parse_reference(inputs: Mapping[str, Any], *, root: Path) -> Path:
    reference = _mapping(_required(inputs, "reference", "inputs"), "inputs.reference")
    _reject_unknown(reference, {"path", "format"}, "inputs.reference")
    path = _relative_child(root, _required(reference, "path", "inputs.reference"), "inputs.reference.path")
    if _string(_required(reference, "format", "inputs.reference"), "inputs.reference.format").lower() != "h5ad":
        raise ApplicationConfigError("inputs.reference.format must be h5ad")
    return path


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
    """Use an existing file first, then only the five official templates."""
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
    """Validate one effective document after all CLI/Python overrides."""
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
    _reject_unknown(application, {"svc_type"}, "application")
    svc_type = _string(_required(application, "svc_type", "application"), "application.svc_type")
    if svc_type not in _SVC_TYPES:
        raise ApplicationConfigError("application.svc_type must be one of: sp-SVC, sc-SVC, sc-SVC-sr")

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
    reference_path = _parse_reference(inputs, root=resolved_root)
    pm_path = None
    if "pm_on_cell" in inputs:
        if svc_type != "sc-SVC-sr":
            raise ApplicationConfigError("inputs.pm_on_cell is only allowed for sc-SVC-sr")
        pm = _mapping(inputs["pm_on_cell"], "inputs.pm_on_cell")
        _reject_unknown(pm, {"path"}, "inputs.pm_on_cell")
        pm_path = _relative_child(resolved_root, _required(pm, "path", "inputs.pm_on_cell"), "inputs.pm_on_cell.path")

    anchoring = _mapping(document["global_anchoring"], "global_anchoring")
    _reject_unknown(anchoring, {"broad_column"}, "global_anchoring")
    broad_column = _string(_required(anchoring, "broad_column", "global_anchoring"), "global_anchoring.broad_column")

    refinement = _mapping(document["local_refinement"], "local_refinement")
    subtype = select = strength = None
    if svc_type == "sc-SVC":
        _reject_unknown(refinement, {"subtype_column", "select_cell_type"}, "local_refinement")
        subtype = _string(_required(refinement, "subtype_column", "local_refinement"), "local_refinement.subtype_column")
        raw_select = _required(refinement, "select_cell_type", "local_refinement")
        if not isinstance(raw_select, str) or raw_select.strip().lower() in _ALL_CELL_TYPES:
            raise ApplicationConfigError("local_refinement.select_cell_type must name one concrete broad cell type")
        select = raw_select.strip()
    else:
        _reject_unknown(refinement, {"strength"}, "local_refinement")
        if "strength" in refinement:
            strength = _number(refinement["strength"], "local_refinement.strength")

    output = _mapping(document["output"], "output")
    _reject_unknown(output, {"dir", "name"}, "output")
    output_dir = _relative_child(resolved_root, _required(output, "dir", "output"), "output.dir")
    output_name = _string(_required(output, "name", "output"), "output.name")
    if output_name in {".", ".."} or "/" in output_name or "\\" in output_name:
        raise ApplicationConfigError("output.name must be a filename stem without separators")
    if output_name.endswith(".h5ad"):
        raise ApplicationConfigError("output.name must not include .h5ad")

    execution = _mapping(document["execution"], "execution")
    if "action" in execution:
        raise ApplicationConfigError("execution.action was removed; use --dry-run or dry_run=True")
    _reject_unknown(execution, {"seed"}, "execution")
    seed = execution.get("seed")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ApplicationConfigError("execution.seed must be an integer")
    if seed is not None and not 0 <= seed <= 2**32 - 1:
        raise ApplicationConfigError("execution.seed must be between 0 and 4294967295")

    return ApplicationConfig(
        config_path=source.path,
        config_sha256=sha256(source.payload).hexdigest(),
        config_source=source.label,
        declared_root=declared_root,
        resolved_root=resolved_root,
        cwd=cwd,
        svc_type=svc_type,
        st_path=st_path,
        reference_path=reference_path,
        st_format=st_format,
        spatialdata_table=table,
        spatialdata_element=element,
        broad_column=broad_column,
        subtype_column=subtype,
        select_cell_type=select,
        local_refinement_strength=strength,
        ot_method=ot_method,
        pm_on_cell_path=pm_path,
        output_dir=output_dir,
        output_name=output_name,
        seed=seed,
    )


__all__ = [
    "ApplicationConfig",
    "ApplicationConfigError",
    "ApplicationConfigSource",
    "compile_application_config",
    "load_application_yaml",
    "resolve_application_source",
]
