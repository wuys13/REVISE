from __future__ import annotations

from collections.abc import Sequence
import importlib

from anndata import AnnData


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


def score_gene_set_aucell(
    adata: AnnData,
    genes: Sequence[str],
    *,
    score_name: str,
) -> tuple[AnnData, str]:
    """Run OmicVerse AUCell on a copy and validate its result column."""
    overlapping = [gene for gene in genes if gene in adata.var_names]
    if not overlapping:
        raise ValueError("Gene set has no overlap with adata.var_names")

    omicverse = _require_omicverse()
    work = adata.copy()
    omicverse.single.geneset_aucell(
        adata=work,
        geneset_name=score_name,
        geneset=overlapping,
    )

    score_key = f"{score_name}_aucell"
    if score_key not in work.obs:
        raise RuntimeError(f"AUCell did not create {score_key!r}")
    return work, score_key
