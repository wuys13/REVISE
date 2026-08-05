"""Compile one application YAML into the request consumed by the service.

This module owns YAML parsing and request normalization only. The public CLI
lives in ``revise.application.cli`` and pipeline execution lives in
``service.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
from pathlib import Path
from typing import Any, Mapping

import yaml


class ApplicationConfigError(ValueError):
    """Raised when an application YAML does not satisfy the public contract."""


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


@dataclass(frozen=True)
class ApplicationRequest:
    config_path: Path | None
    config_sha256: str
    declared_root: str
    resolved_root: Path
    cwd: Path
    svc_type: str
    sample_name: str
    input_mode: str
    st_path: Path
    reference_path: Path
    st_format: str
    patient_key: str
    broad_column: str
    subtype_column: str | None
    select_cell_type: str | None
    local_refinement_strength: float | None
    ot_method: str | None
    output_root: Path
    action: str
    effective_action: str
    dry_run_override: bool
    seed: int | None
    data_root: Path | None = None
    st_file: str | None = None
    reference_file: str | None = None
    spatialdata_table: str | None = None
    spatialdata_element: str | None = None
    config_source: str | None = None

    @property
    def source_path(self) -> str:
        return self.config_source or str(self.config_path)

    @property
    def resolved_paths(self) -> dict[str, str]:
        return {
            "st": str(self.st_path),
            "reference": str(self.reference_path),
            "output": str(self.output_root),
        }


def _mapping(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ApplicationConfigError(f"{field} must be a mapping")
    non_string = [key for key in value if not isinstance(key, str)]
    if non_string:
        rendered = ", ".join(repr(key) for key in non_string)
        raise ApplicationConfigError(
            f"{field} contains non-string field name(s): {rendered}"
        )
    return dict(value)


def _reject_unknown(mapping: Mapping[str, Any], allowed: set[str], field: str) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
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


def _safe_output_component(value: Any, field: str) -> str:
    component = _string(value, field)
    if component in {".", ".."} or "/" in component or "\\" in component:
        raise ApplicationConfigError(f"{field} must be a safe output path component")
    if "\x00" in component:
        raise ApplicationConfigError(f"{field} must not contain NUL")
    return component


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ApplicationConfigError(
            f"{field} must be a non-negative finite number"
        )
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ApplicationConfigError(
            f"{field} must be a non-negative finite number"
        )
    return result


def _root_dir(value: Any, *, cwd: Path) -> tuple[str, Path]:
    if not isinstance(value, str) or not value:
        raise ApplicationConfigError(
            "paths.root_dir must be literal . or an existing absolute directory"
        )
    if value.startswith("~") or value != value.strip():
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
    if value.startswith("~"):
        raise ApplicationConfigError(f"{field} must not use ~ expansion")
    child = Path(value)
    if child.is_absolute():
        raise ApplicationConfigError(f"{field} must be relative to paths.root_dir")
    if ".." in child.parts:
        raise ApplicationConfigError(f"{field} escapes paths.root_dir")
    resolved = root / child
    if resolved == root:
        raise ApplicationConfigError(f"{field} must be a child of paths.root_dir")
    return resolved


def _parse_st(
    inputs: Mapping[str, Any],
    *,
    mode: str,
    root: Path,
    data_root_value: str | None,
    sample_name: str,
) -> tuple[Path, str | None, str, str | None, str | None]:
    st = _mapping(_required(inputs, "st", "inputs"), "inputs.st")
    path_key = "path" if mode == "direct" else "file"
    _reject_unknown(st, {path_key, "format", "spatialdata"}, "inputs.st")
    st_value = _string(_required(st, path_key, "inputs.st"), f"inputs.st.{path_key}")
    st_format = _string(_required(st, "format", "inputs.st"), "inputs.st.format").lower()
    if st_format not in _ST_FORMATS:
        raise ApplicationConfigError(
            "inputs.st.format must be one of: h5ad, spatialdata, auto"
        )

    spatialdata_table = None
    spatialdata_element = None
    if "spatialdata" in st:
        if st_format == "h5ad":
            raise ApplicationConfigError(
                "inputs.st.spatialdata is not allowed when inputs.st.format is h5ad"
            )
        spatialdata = _mapping(st["spatialdata"], "inputs.st.spatialdata")
        _reject_unknown(spatialdata, {"table", "element"}, "inputs.st.spatialdata")
        spatialdata_table = _optional_string(
            spatialdata.get("table"), "inputs.st.spatialdata.table"
        )
        spatialdata_element = _optional_string(
            spatialdata.get("element"), "inputs.st.spatialdata.element"
        )

    if mode == "direct":
        return (
            _relative_child(root, st_value, "inputs.st.path"),
            None,
            st_format,
            spatialdata_table,
            spatialdata_element,
        )
    assert data_root_value is not None
    _relative_child(root, st_value, "inputs.st.file")
    return (
        _relative_child(
            root,
            str(Path(data_root_value) / f"{sample_name}_{st_value}"),
            "inputs.st.file",
        ),
        st_value,
        st_format,
        spatialdata_table,
        spatialdata_element,
    )


def _parse_reference(
    inputs: Mapping[str, Any],
    *,
    mode: str,
    root: Path,
    data_root_value: str | None,
) -> tuple[Path, str | None, str]:
    reference = _mapping(
        _required(inputs, "reference", "inputs"), "inputs.reference"
    )
    path_key = "path" if mode == "direct" else "file"
    _reject_unknown(
        reference,
        {path_key, "format", "patient_key"},
        "inputs.reference",
    )
    reference_value = _string(
        _required(reference, path_key, "inputs.reference"),
        f"inputs.reference.{path_key}",
    )
    reference_format = _string(
        _required(reference, "format", "inputs.reference"),
        "inputs.reference.format",
    ).lower()
    if reference_format != "h5ad":
        raise ApplicationConfigError("inputs.reference.format must be h5ad")
    patient_key = _string(
        _required(reference, "patient_key", "inputs.reference"),
        "inputs.reference.patient_key",
    )
    if mode == "direct":
        return (
            _relative_child(root, reference_value, "inputs.reference.path"),
            None,
            patient_key,
        )
    assert data_root_value is not None
    _relative_child(root, reference_value, "inputs.reference.file")
    return (
        _relative_child(
            root,
            str(Path(data_root_value) / reference_value),
            "inputs.reference.file",
        ),
        reference_value,
        patient_key,
    )


def load_application_request(
    config_path: str | Path | None,
    *,
    dry_run: bool = False,
    _payload: bytes | None = None,
    _source_path: str | None = None,
) -> ApplicationRequest:
    """Load, validate, and resolve one application YAML."""
    cwd = Path.cwd().resolve()
    path = Path(config_path).expanduser().resolve() if config_path is not None else None
    source_path = _source_path or str(path)
    payload = _payload
    if _payload is None:
        try:
            payload = path.read_bytes()
        except OSError as exc:
            raise ApplicationConfigError(
                f"Cannot read application config {source_path}: {exc}"
            ) from exc
    try:
        raw = yaml.safe_load(payload)
    except yaml.YAMLError as exc:
        raise ApplicationConfigError(
            f"Invalid application YAML {source_path}: {exc}"
        ) from exc
    document = _mapping(raw, "application config")
    if "schema_version" not in document and any(
        key in document for key in ("version", "defaults", "router", "profiles")
    ):
        raise ApplicationConfigError(
            "--config requires an application YAML; the full engine config is "
            "internally managed by the package"
        )
    _reject_unknown(document, _TOP_LEVEL_KEYS, "application config")
    missing_sections = sorted(_TOP_LEVEL_KEYS - set(document))
    if missing_sections:
        raise ApplicationConfigError(
            f"application config is missing required field(s): {', '.join(missing_sections)}"
        )

    schema_version = document["schema_version"]
    if isinstance(schema_version, bool) or not isinstance(schema_version, int):
        raise ApplicationConfigError("schema_version must be integer 1")
    if schema_version != 1:
        raise ApplicationConfigError("schema_version must be 1")

    application = _mapping(document["application"], "application")
    _reject_unknown(application, {"svc_type", "sample_name"}, "application")
    svc_type = _string(
        _required(application, "svc_type", "application"), "application.svc_type"
    )
    if svc_type not in _SVC_TYPES:
        raise ApplicationConfigError(
            "application.svc_type must be one of: sp-SVC, sc-SVC, sc-SVC-sr"
        )
    sample_name = _safe_output_component(
        _required(application, "sample_name", "application"),
        "application.sample_name",
    )

    paths = _mapping(document["paths"], "paths")
    _reject_unknown(paths, {"root_dir"}, "paths")
    declared_root, resolved_root = _root_dir(
        _required(paths, "root_dir", "paths"), cwd=cwd
    )

    algorithm = _mapping(document["algorithm"], "algorithm")
    if "base_config" in algorithm:
        raise ApplicationConfigError(
            "algorithm.base_config is not a public field; the engine config is "
            "internally managed by the package"
        )
    _reject_unknown(algorithm, {"ot_method"}, "algorithm")
    ot_method = _optional_string(algorithm.get("ot_method"), "algorithm.ot_method")
    if ot_method is not None and ot_method not in {"pot", "tacco"}:
        raise ApplicationConfigError("algorithm.ot_method must be pot or tacco")

    inputs = _mapping(document["inputs"], "inputs")
    mode = _string(_required(inputs, "mode", "inputs"), "inputs.mode")
    if mode not in {"direct", "legacy_layout"}:
        raise ApplicationConfigError("inputs.mode must be direct or legacy_layout")
    allowed_inputs = {"mode", "st", "reference"}
    data_root_value = None
    data_root = None
    if mode == "legacy_layout":
        allowed_inputs.add("data_root")
        data_root_value = _string(
            _required(inputs, "data_root", "inputs"), "inputs.data_root"
        )
        data_root = _relative_child(
            resolved_root, data_root_value, "inputs.data_root"
        )
    _reject_unknown(inputs, allowed_inputs, "inputs")
    st_path, st_file, st_format, spatialdata_table, spatialdata_element = _parse_st(
        inputs,
        mode=mode,
        root=resolved_root,
        data_root_value=data_root_value,
        sample_name=sample_name,
    )
    reference_path, reference_file, patient_key = _parse_reference(
        inputs,
        mode=mode,
        root=resolved_root,
        data_root_value=data_root_value,
    )

    anchoring = _mapping(document["global_anchoring"], "global_anchoring")
    _reject_unknown(anchoring, {"broad_column"}, "global_anchoring")
    broad_column = _string(
        _required(anchoring, "broad_column", "global_anchoring"),
        "global_anchoring.broad_column",
    )

    refinement = _mapping(document["local_refinement"], "local_refinement")
    subtype_column = None
    select_cell_type = None
    local_refinement_strength = None
    if svc_type == "sc-SVC":
        if "strength" in refinement:
            raise ApplicationConfigError(
                "local_refinement.strength is not allowed for sc-SVC"
            )
        _reject_unknown(
            refinement,
            {"subtype_column", "select_cell_type"},
            "local_refinement",
        )
        subtype_column = _string(
            _required(refinement, "subtype_column", "local_refinement"),
            "local_refinement.subtype_column",
        )
        raw_cell_type = _required(
            refinement, "select_cell_type", "local_refinement"
        )
        normalized = raw_cell_type.strip().lower() if isinstance(raw_cell_type, str) else None
        if normalized is None or normalized in _ALL_CELL_TYPES:
            raise ApplicationConfigError(
                "local_refinement.select_cell_type must name one concrete broad cell type"
            )
        select_cell_type = _safe_output_component(
            raw_cell_type, "local_refinement.select_cell_type"
        )
    else:
        _reject_unknown(refinement, {"strength"}, "local_refinement")
        if "strength" in refinement:
            local_refinement_strength = _number(
                refinement["strength"], "local_refinement.strength"
            )

    output = _mapping(document["output"], "output")
    _reject_unknown(output, {"path"}, "output")
    output_root = _relative_child(
        resolved_root,
        _required(output, "path", "output"),
        "output.path",
    )

    execution = _mapping(document["execution"], "execution")
    _reject_unknown(execution, {"action", "seed"}, "execution")
    action = _string(
        _required(execution, "action", "execution"), "execution.action"
    )
    if action not in {"run", "preflight"}:
        raise ApplicationConfigError("execution.action must be run or preflight")
    seed = execution.get("seed")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ApplicationConfigError("execution.seed must be an integer")
    if seed is not None and not 0 <= seed <= 2**32 - 1:
        raise ApplicationConfigError(
            "execution.seed must be between 0 and 4294967295"
        )
    dry_run_override = bool(dry_run and action == "run")
    effective_action = "preflight" if dry_run_override else action

    return ApplicationRequest(
        config_path=path,
        config_sha256=sha256(payload).hexdigest(),
        declared_root=declared_root,
        resolved_root=resolved_root,
        cwd=cwd,
        svc_type=svc_type,
        sample_name=sample_name,
        input_mode=mode,
        st_path=st_path,
        reference_path=reference_path,
        st_format=st_format,
        patient_key=patient_key,
        broad_column=broad_column,
        subtype_column=subtype_column,
        select_cell_type=select_cell_type,
        local_refinement_strength=local_refinement_strength,
        ot_method=ot_method,
        output_root=output_root,
        action=action,
        effective_action=effective_action,
        dry_run_override=dry_run_override,
        seed=seed,
        data_root=data_root,
        st_file=st_file,
        reference_file=reference_file,
        spatialdata_table=spatialdata_table,
        spatialdata_element=spatialdata_element,
        config_source=_source_path,
    )
