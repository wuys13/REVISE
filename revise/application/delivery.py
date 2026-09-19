"""Whole-sample Raw and analysis handoff, without changing input annotations."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

import pandas as pd
import numpy as np

from revise.utils.provenance import input_identities
from .config import ApplicationConfig
from .expression import consumer_declaration, is_known_linear


DELIVERY_VERSION = 2
LABEL_COLUMNS = {"broad": "revise_Level1", "subtype": "revise_Level2"}


def is_sample_delivery(config):
    return isinstance(config, ApplicationConfig) and (
        config.mode != "cluster" or config.select_cell_type is None
    )


def source_identity(config):
    if not config.st_path.exists():
        raise ValueError("Complete Raw requires an original spatial object or an existing input source")
    return input_identities([SimpleNamespace(role="raw", path=config.st_path)])[0]


def validate_raw(raw):
    ids = raw.obs_names
    if not ids.is_unique or ids.isna().any() or any(not str(i).strip() for i in ids):
        raise ValueError("Raw observation IDs must be unique and nonempty")
    occupied = set(LABEL_COLUMNS.values()).intersection(raw.obs.columns)
    if occupied:
        raise ValueError(f"Original Raw already contains reserved inference columns: {sorted(occupied)}")
    if "revise_delivery" in raw.uns:
        raise ValueError("Original Raw already contains reserved metadata namespace: revise_delivery")


def _inferred_labels(config, ctx):
    """Only runner spatial units are eligible for backfilling Raw labels."""
    runner = getattr(ctx, "runner", None)
    anchored = getattr(runner, "st_adata", None)
    outputs = ctx.svc.artifacts.get("outputs", {})
    refined = outputs.get("sc_svc_spatial") if config.mode == "cluster" else None
    sources = [(config.broad_column, LABEL_COLUMNS["broad"], anchored)]
    # A pre-existing subtype column on the input is not evidence of an LR inference.
    if refined is not None and config.subtype_column:
        sources.append((config.subtype_column, LABEL_COLUMNS["subtype"], refined))
    return sources


def prepare_raw(config, raw, svc, ctx, *, owned=False):
    validate_raw(raw)
    if not owned:
        raw = raw.copy()
    conflicts = {}
    for source, target, inferred in _inferred_labels(config, ctx):
        raw.obs[target] = pd.Categorical([None] * raw.n_obs)
        if inferred is None or source not in inferred.obs:
            continue
        ids = inferred.obs_names
        if not ids.is_unique or not ids.isin(raw.obs_names).all():
            raise ValueError("Inferred observation IDs do not uniquely align to original Raw")
        values = inferred.obs[source].astype("string").reindex(raw.obs_names)
        raw.obs[target] = pd.Categorical(values.to_numpy(dtype=object, na_value=None))
        if source in raw.obs:
            original = raw.obs[source].astype("string")
            conflicts[source] = int((original.notna() & values.notna() & original.ne(values)).sum())
    # sST assigns virtual-cell broad labels in cell_type, not the spot source column.
    svc_broad = "cell_type" if config.mode == "sr" and "cell_type" in svc.obs else config.broad_column
    # Publish inference aliases on SVC without replacing its existing labels.
    for source, target in ((svc_broad, LABEL_COLUMNS["broad"]),
                           (config.subtype_column, LABEL_COLUMNS["subtype"])):
        if source and source in svc.obs:
            if target in svc.obs:
                raise ValueError(f"SVC already contains reserved inference column: {target}")
            svc.obs[target] = pd.Categorical(svc.obs[source].astype(object))
    raw.uns["revise_delivery"] = {
        "version": DELIVERY_VERSION,
        "annotation_sources": {key: column for key, column in LABEL_COLUMNS.items()
                               if column in raw.obs},
        "original_label_conflicts": conflicts,
    }
    return raw


def sample_document(config, paths, raw, svc, *, source=None, sample_id=None, coordinates=None):
    original_id = sample_id or config.output_dir.name
    # Encode percent too, making hierarchical IDs reversible and collision-free.
    encoded_id = quote(original_id, safe="-_").replace(".", "%2E")
    if not encoded_id:
        raise ValueError("A complete delivery requires a nonempty sample identity")
    raw_expression = consumer_declaration(config.st_expression, raw)
    svc_source = None
    svc_identity = None
    if config.mode == "cluster":
        svc_source = config.reference_expression
        if is_known_linear(svc_source):
            svc_identity = f"reference_expression_{config.ist_mapping}"
    elif config.mode == "sr":
        svc_identity = "sst_parent_spot_corrected_expression"
        if is_known_linear(config.st_expression) and is_known_linear(config.reference_expression):
            svc_source = config.st_expression
    elif is_known_linear(config.st_expression) and is_known_linear(config.reference_expression):
        svc_source = config.st_expression
        svc_identity = "hst_reconstructed_expression"
    svc_expression = consumer_declaration(svc_source, svc, identity=svc_identity)
    columns = {key: column for key, column in LABEL_COLUMNS.items()
               if column in raw.obs and column in svc.obs}
    if "SVC_cluster" in svc.obs:
        columns["reconstruction"] = "SVC_cluster"
    spatial = {"key": "spatial", "unit": "unknown"}
    if coordinates:
        spatial.update({key: value for key, value in coordinates.items()
                        if key in {"key", "unit", "microns_per_coordinate"} and value is not None})
    if spatial["unit"] == "um":
        spatial["unit"] = "micron"
    for role, obj in (("Raw", raw), ("SVC", svc)):
        key = spatial["key"]
        if key not in obj.obsm:
            raise ValueError(f"{role} is missing declared spatial coordinates obsm[{key!r}]")
        values = np.asarray(obj.obsm[key])
        if (values.ndim != 2 or values.shape[0] != obj.n_obs or values.shape[1] < 2
                or not np.isfinite(values).all()):
            raise ValueError(f"{role} declared spatial coordinates must be finite with shape (n_obs, >=2)")
    return {
        "schema_version": 1, "sample_id": encoded_id,
        "files": {role: Path(paths[role]).name for role in ("raw", "svc")},
        "columns": columns, "expression": {"raw": raw_expression, "svc": svc_expression},
        "spatial": spatial,
        "provenance": {"delivery_version": DELIVERY_VERSION, "original_sample_id": original_id,
                       "sample_id_encoding": "percent-encoded UTF-8 including dots",
                       "raw_source": source or {"kind": "explicit_original_object"},
                       "route": config.mode or config.svc_type,
                       "svc_expression_processing": (
                           "internal normalization followed by parent-spot per-gene correction; not raw counts"
                           if config.mode == "sr" else
                           "reference carrier X as supplied; no assembly normalization"
                           if config.mode == "cluster" else "route-specific reconstruction"),
                       "expression_sources": {"spatial": config.st_expression or {"identity": "unknown"},
                                              "reference": config.reference_expression or {"identity": "unknown"}}},
    }
