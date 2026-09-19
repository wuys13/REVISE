"""Strict configuration contracts for reference preparation and consumption."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import yaml


class ConfigError(ValueError):
    """A preparation or resolved-reference configuration is invalid."""


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    path: Path


@dataclass(frozen=True)
class PairedConfig:
    source: Path
    pair_column: str
    pair_key: str
    output_dir: Path


@dataclass(frozen=True)
class ScreenConfig:
    reconstruction_config: Path
    candidates: tuple[Candidate, ...]
    output_dir: Path


@dataclass(frozen=True)
class ResolvedReference:
    path: Path
    format: str = "h5ad"


class _UniqueKeyLoader(yaml.SafeLoader):
    pass


def _unique_mapping(loader: yaml.SafeLoader, node: yaml.MappingNode, deep: bool = False):
    result: dict[str, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str):
            raise ConfigError("YAML mapping keys must be strings")
        if key in result:
            raise ConfigError(f"Duplicate YAML key: {key}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping
)


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a mapping")
    return value


def _keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    missing = expected - value.keys()
    unknown = value.keys() - expected
    if missing or unknown:
        raise ConfigError(
            f"{label}: missing fields {sorted(missing)}; unknown fields {sorted(unknown)}"
        )


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{label} must be a non-empty string")
    return value


def _path(value: Any, base: Path, label: str) -> Path:
    path = Path(_text(value, label)).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def _read_yaml(path: str | Path) -> tuple[Path, dict[str, Any]]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ConfigError(f"Configuration file does not exist: {source}")
    try:
        document = yaml.load(source.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except yaml.YAMLError as error:
        raise ConfigError(f"Invalid YAML in {source}: {error}") from error
    return source, _mapping(document, str(source))


def _version(document: dict[str, Any]) -> None:
    if type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise ConfigError("schema_version must be the integer 1")


def _candidate_list(source: Path) -> tuple[Candidate, ...]:
    path, document = _read_yaml(source)
    _keys(document, {"schema_version", "candidates"}, "candidate list")
    _version(document)
    entries = document["candidates"]
    if not isinstance(entries, list):
        raise ConfigError("candidates must be a list")
    candidates = []
    for entry in entries:
        item = _mapping(entry, "candidate")
        _keys(item, {"id", "path"}, "candidate")
        candidates.append(Candidate(
            _text(item["id"], "candidate.id"),
            _path(item["path"], path.parent, "candidate.path"),
        ))
    return tuple(candidates)


def _candidates(value: Any, base: Path) -> tuple[Candidate, ...]:
    source = _mapping(value, "candidates")
    if set(source) == {"directory"}:
        directory = _path(source["directory"], base, "candidates.directory")
        if not directory.is_dir():
            raise ConfigError(f"Candidate directory does not exist: {directory}")
        candidates = tuple(
            Candidate(path.name, path.resolve())
            for path in sorted(directory.iterdir())
            if path.is_file() and path.suffix.lower() == ".h5ad"
        )
    elif set(source) == {"list"}:
        candidates = _candidate_list(_path(source["list"], base, "candidates.list"))
    else:
        raise ConfigError("candidates requires exactly one of directory or list")
    ids = [candidate.candidate_id for candidate in candidates]
    paths = [candidate.path for candidate in candidates]
    if len(ids) != len(set(ids)):
        raise ConfigError("Candidate IDs must be unique")
    if len(paths) != len(set(paths)):
        raise ConfigError("Candidate paths must be unique after resolution")
    return candidates


def load_preparation_config(path: str | Path) -> PairedConfig | ScreenConfig:
    source, document = _read_yaml(path)
    _version(document)
    base = source.parent
    if document.get("mode") == "paired":
        _keys(document, {
            "schema_version", "mode", "source", "pair_column", "pair_key", "output_dir",
        }, "paired configuration")
        pair_key = document["pair_key"]
        if not isinstance(pair_key, str) or not pair_key:
            raise ConfigError("pair_key must be a non-empty string")
        return PairedConfig(
            source=_path(document["source"], base, "source"),
            pair_column=_text(document["pair_column"], "pair_column"),
            pair_key=pair_key,
            output_dir=_path(document["output_dir"], base, "output_dir"),
        )
    if document.get("mode") == "screen":
        _keys(document, {
            "schema_version", "mode", "reconstruction_config", "candidates", "output_dir",
        }, "screen configuration")
        reconstruction = _path(document["reconstruction_config"], base, "reconstruction_config")
        if not reconstruction.is_file():
            raise ConfigError(f"Reconstruction configuration does not exist: {reconstruction}")
        return ScreenConfig(
            reconstruction_config=reconstruction,
            candidates=_candidates(document["candidates"], base),
            output_dir=_path(document["output_dir"], base, "output_dir"),
        )
    raise ConfigError("mode must be paired or screen")


def _reference_file(path: Path) -> None:
    if path.suffix.lower() != ".h5ad" or not path.is_file():
        raise ConfigError(f"Reference must be an existing H5AD file: {path}")


def load_reference_config(path: str | Path) -> ResolvedReference:
    """Resolve the minimal reference-only YAML independently of cwd."""
    source, document = _read_yaml(path)
    _keys(document, {"schema_version", "reference"}, "resolved reference")
    _version(document)
    reference = _mapping(document["reference"], "reference")
    _keys(reference, {"path", "format"}, "reference")
    if reference["format"] != "h5ad":
        raise ConfigError("reference.format must be h5ad")
    resolved = _path(reference["path"], source.parent, "reference.path")
    _reference_file(resolved)
    return ResolvedReference(path=resolved)


def _atomic_write_text(path: Path, content: str, *, replace: bool = False) -> None:
    """Publish atomically, refusing to overwrite unless updating our own report."""
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        if replace:
            temporary.replace(path)
        else:
            os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_reference_config(path: Path, reference_path: Path) -> None:
    reference_path = reference_path.resolve()
    _reference_file(reference_path)
    relative = os.path.relpath(reference_path, path.resolve().parent)
    _atomic_write_text(path, yaml.safe_dump({
        "schema_version": 1,
        "reference": {"path": relative, "format": "h5ad"},
    }, sort_keys=False, allow_unicode=True))


def atomic_write_json(path: Path, payload: dict[str, Any], *, replace: bool = False) -> None:
    _atomic_write_text(path, json.dumps(
        payload, indent=2, ensure_ascii=False, allow_nan=False,
    ) + "\n", replace=replace)
