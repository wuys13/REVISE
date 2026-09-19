"""Explicit source declarations for the independent linear-expression consumer.

Declarations attest to processing history; numerical checks alone cannot detect log data.
"""
from dataclasses import replace
from collections.abc import Mapping

import numpy as np
from scipy import sparse


LINEAR_SCALE = "untransformed_nonnegative"


def parse_expression_declaration(value, field="expression"):
    from .config import ApplicationConfigError

    if value is None:
        return None
    if not isinstance(value, Mapping) or set(value) - {"identity", "scale"}:
        raise ApplicationConfigError(f"{field} accepts identity and scale only")
    identity = value.get("identity", "unknown")
    scale = value.get("scale", "unknown")
    if not isinstance(identity, str) or not identity.strip():
        raise ApplicationConfigError(f"{field}.identity must be a nonempty source description")
    if not isinstance(scale, str) or scale not in {"unknown", LINEAR_SCALE}:
        raise ApplicationConfigError(
            f"{field}.scale must be unknown or {LINEAR_SCALE}; log expression is not supported"
        )
    return {"identity": identity.strip(), "scale": scale}


def is_known_linear(declaration):
    return bool(declaration and declaration["identity"] != "unknown"
                and declaration["scale"] == LINEAR_SCALE)


def validate_linear_matrix(adata, role):
    if adata.X is None:
        raise ValueError(f"{role}.X is missing despite its linear expression declaration")
    for start in range(0, adata.n_obs, 4096):
        block = adata.X[start:start + 4096]
        values = block.data if sparse.issparse(block) else np.asarray(block)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{role}.X must be finite and nonnegative")


def bind_expression_sources(config, spatial, reference):
    """Use explicit config, else the source object's explicit revise_expression record."""
    declarations = {}
    for field, obj in (("st_expression", spatial), ("reference_expression", reference)):
        declaration = getattr(config, field, None)
        if declaration is None:
            declaration = parse_expression_declaration(
                obj.uns.get("revise_expression"), f"{field}.uns.revise_expression"
            )
        if is_known_linear(declaration):
            validate_linear_matrix(obj, field)
        declarations[field] = declaration
    return replace(config, **declarations)


def consumer_declaration(declaration, adata, *, identity=None):
    if not is_known_linear(declaration):
        source_identity = declaration["identity"] if declaration else "unknown"
        return {"matrix": "X", "identity": identity or source_identity, "scale": "unknown"}
    validate_linear_matrix(adata, identity or declaration["identity"])
    # Current consumer assumes linear X for a known identity without a legacy scale key.
    return {"matrix": "X", "identity": identity or declaration["identity"]}
