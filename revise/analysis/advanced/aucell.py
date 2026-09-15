from __future__ import annotations

from collections.abc import Sequence
import importlib
from importlib import metadata

import numpy as np
from anndata import AnnData
from scipy import sparse


def _require_omicverse():
    try:
        return importlib.import_module("omicverse")
    except (ImportError, ModuleNotFoundError) as exc:
        if isinstance(exc, ModuleNotFoundError) and exc.name != "omicverse":
            raise ImportError(
                f"OmicVerse import failed because dependency {exc.name!r} is unavailable"
            ) from exc
        raise ImportError(
            "OmicVerse is required for AUCell scoring. Install "
            "`revise-svc[pathway]` to use this capability."
        ) from exc


def get_aucell_provider_metadata() -> dict[str, str]:
    """Return the loaded OmicVerse provider identity used by the scorer."""
    omicverse = _require_omicverse()
    provider = getattr(getattr(omicverse, "single", None), "geneset_aucell", None)
    if provider is None or not callable(provider):
        raise RuntimeError("OmicVerse single.geneset_aucell provider is unavailable")
    version = getattr(omicverse, "__version__", None)
    if not version:
        try:
            version = metadata.version("omicverse")
        except metadata.PackageNotFoundError:
            version = "unknown"
    return {
        "name": "omicverse",
        "version": str(version),
        "scorer": "single.geneset_aucell",
        "provider_module": str(getattr(provider, "__module__", "unknown")),
    }


def score_gene_set_aucell(
    adata: AnnData,
    genes: Sequence[str],
    *,
    score_name: str,
    AUC_threshold: float | None = None,
    seed: int | None = None,
) -> tuple[AnnData, str]:
    """Run OmicVerse AUCell on a copy and validate its result column."""
    overlapping = [gene for gene in genes if gene in adata.var_names]
    if not overlapping:
        raise ValueError("Gene set has no overlap with adata.var_names")

    omicverse = _require_omicverse()
    work = adata.copy()
    # OmicVerse 1.7.5 derives its threshold and ranks through sparse-matrix
    # methods.  Preserve the supplied values while giving that provider its
    # required carrier type.
    if not sparse.issparse(work.X):
        work.X = sparse.csr_matrix(work.X)
    provider_kwargs = {
        "adata": work,
        "geneset_name": score_name,
        "geneset": overlapping,
    }
    if AUC_threshold is not None:
        provider_kwargs["AUC_threshold"] = AUC_threshold
    if seed is not None:
        provider_kwargs["seed"] = seed
    omicverse.single.geneset_aucell(**provider_kwargs)

    score_key = f"{score_name}_aucell"
    if score_key not in work.obs:
        raise RuntimeError(f"AUCell did not create {score_key!r}")
    scores = np.asarray(work.obs[score_key], dtype=float)
    if scores.shape != (work.n_obs,) or not np.isfinite(scores).all():
        raise RuntimeError(f"AUCell returned non-finite scores in {score_key!r}")
    return work, score_key
