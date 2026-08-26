from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Iterable, Mapping
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SEED = 42
MAX_OBSERVATIONS_PER_SOURCE = 10_000


def find_repo_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "reconstruct.py").is_file() and (
            candidate / "reproduce" / "case" / "tacco"
        ).is_dir():
            return candidate
    raise RuntimeError(f"Cannot locate the REVISE repository from {current}")


def data_root() -> Path:
    # A case notebook must declare its data root explicitly.  Keeping the
    # environment lookup fail-closed avoids silently reading a checkout-local
    # or machine-specific path when a notebook is moved to another host.
    return Path(os.environ["REVISE_TACCO_DATA_ROOT"]).expanduser().resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def matrix_values_are_finite_nonnegative(matrix) -> tuple[bool, bool]:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.isfinite(values).all()), bool((values >= 0).all())


def audit_adata(name: str, adata: ad.AnnData) -> dict:
    finite, nonnegative = matrix_values_are_finite_nonnegative(adata.X)
    spatial = adata.obsm.get("spatial")
    return {
        "object": name,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "unique_obs_names": bool(adata.obs_names.is_unique),
        "unique_var_names": bool(adata.var_names.is_unique),
        "finite_X": finite,
        "nonnegative_X": nonnegative,
        "has_spatial": spatial is not None,
        "finite_spatial": bool(np.isfinite(np.asarray(spatial)).all())
        if spatial is not None
        else None,
    }


def audit_table(objects: Mapping[str, ad.AnnData]) -> pd.DataFrame:
    return pd.DataFrame([audit_adata(name, obj) for name, obj in objects.items()])


def manifest_from_output(adata: ad.AnnData) -> tuple[Path, dict]:
    metadata = adata.uns.get("revise_reconstruction")
    if not isinstance(metadata, Mapping) or "run_manifest" not in metadata:
        raise KeyError("Output lacks uns['revise_reconstruction']['run_manifest']")
    path = Path(str(metadata["run_manifest"]))
    manifest = load_json(path)
    run = manifest.get("run", {})
    if (
        run.get("status") != "succeeded"
        or run.get("dry_run") is not False
        or not run.get("ended_at")
        or run.get("error") is not None
    ):
        raise RuntimeError(f"REVISE run is not terminal-success: {path}")
    artifact_by_role = {}
    for artifact in manifest.get("artifacts", []):
        role = str(artifact.get("role", ""))
        artifact_path = Path(str(artifact.get("path", "")))
        if artifact.get("status") != "completed" or not artifact_path.is_file():
            raise RuntimeError(f"Incomplete REVISE artifact in {path}: {artifact}")
        observed_size = artifact_path.stat().st_size
        observed_sha256 = sha256_file(artifact_path)
        if observed_size != int(artifact.get("size", -1)) or observed_sha256 != artifact.get(
            "sha256"
        ):
            raise RuntimeError(f"REVISE artifact identity mismatch: {artifact_path}")
        artifact_by_role[role] = artifact
    output_role = str(metadata.get("output_role", ""))
    if f"publication:{output_role}" not in artifact_by_role:
        raise RuntimeError(
            f"Published output role {output_role!r} is absent from {path}"
        )
    return path, manifest


def _sample(adata: ad.AnnData, *, max_observations: int, seed: int) -> ad.AnnData:
    if adata.n_obs <= max_observations:
        return adata.copy()
    rng = np.random.default_rng(seed)
    selected = np.sort(rng.choice(adata.n_obs, size=max_observations, replace=False))
    return adata[selected].copy()


def _total_counts(adata: ad.AnnData) -> np.ndarray:
    return np.asarray(adata.X.sum(axis=1)).ravel().astype(float)


def joint_umap(
    sources: Mapping[str, ad.AnnData],
    *,
    categorical_labels: Mapping[str, Iterable[str] | pd.Series | np.ndarray],
    max_observations: int = MAX_OBSERVATIONS_PER_SOURCE,
    seed: int = SEED,
) -> ad.AnnData:
    if set(sources) != set(categorical_labels):
        raise ValueError("Each UMAP source must have one explicit display-label vector")
    common = None
    for obj in sources.values():
        names = pd.Index(obj.var_names)
        common = names if common is None else common.intersection(names, sort=False)
    if common is None or len(common) < 2:
        raise ValueError("Joint UMAP requires at least two shared genes")

    prepared: dict[str, ad.AnnData] = {}
    for offset, (name, obj) in enumerate(sources.items()):
        subset = obj[:, common].copy()
        subset.obs["_display_label"] = pd.Series(
            np.asarray(list(categorical_labels[name]), dtype=str),
            index=subset.obs_names,
            dtype="object",
        )
        subset.obs["_total_counts"] = _total_counts(subset)
        subset = _sample(
            subset,
            max_observations=max_observations,
            seed=seed + offset,
        )
        prepared[name] = subset

    combined = ad.concat(
        prepared,
        join="inner",
        label="_source",
        index_unique="::",
        merge="same",
    )
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    n_top = min(2_000, combined.n_vars)
    if n_top < combined.n_vars:
        sc.pp.highly_variable_genes(combined, n_top_genes=n_top, subset=True)
    n_comps = min(30, combined.n_obs - 1, combined.n_vars - 1)
    if n_comps < 2:
        raise ValueError("Joint UMAP has fewer than two usable principal components")
    sc.pp.pca(combined, n_comps=n_comps, random_state=seed)
    n_neighbors = min(15, combined.n_obs - 1)
    sc.pp.neighbors(
        combined,
        n_neighbors=n_neighbors,
        n_pcs=n_comps,
        random_state=seed,
    )
    sc.tl.umap(combined, random_state=seed)
    combined.uns["joint_umap_contract"] = {
        "shared_genes_before_hvg": len(common),
        "genes_after_hvg": int(combined.n_vars),
        "max_observations_per_source": int(max_observations),
        "n_pcs": int(n_comps),
        "n_neighbors": int(n_neighbors),
        "seed": int(seed),
    }
    return combined


def independent_umap(
    adata: ad.AnnData,
    *,
    categorical_labels: Iterable[str] | pd.Series | np.ndarray,
    source_name: str,
    max_observations: int = MAX_OBSERVATIONS_PER_SOURCE,
    seed: int = SEED,
) -> ad.AnnData:
    labels = np.asarray(list(categorical_labels), dtype=str)
    if len(labels) != adata.n_obs:
        raise ValueError(
            f"UMAP labels for {source_name!r} have length {len(labels)}, "
            f"expected {adata.n_obs}"
        )

    prepared = adata.copy()
    prepared.obs["_display_label"] = pd.Series(
        labels,
        index=prepared.obs_names,
        dtype="object",
    )
    prepared.obs["_total_counts"] = _total_counts(prepared)
    input_observations = int(prepared.n_obs)
    input_genes = int(prepared.n_vars)
    prepared = _sample(
        prepared,
        max_observations=max_observations,
        seed=seed,
    )

    sc.pp.normalize_total(prepared, target_sum=1e4)
    sc.pp.log1p(prepared)
    n_top = min(2_000, prepared.n_vars)
    if n_top < prepared.n_vars:
        sc.pp.highly_variable_genes(prepared, n_top_genes=n_top, subset=True)
    n_comps = min(30, prepared.n_obs - 1, prepared.n_vars - 1)
    if n_comps < 2:
        raise ValueError(
            f"Independent UMAP for {source_name!r} has fewer than two usable "
            "principal components"
        )
    sc.pp.pca(prepared, n_comps=n_comps, random_state=seed)
    n_neighbors = min(15, prepared.n_obs - 1)
    sc.pp.neighbors(
        prepared,
        n_neighbors=n_neighbors,
        n_pcs=n_comps,
        random_state=seed,
    )
    sc.tl.umap(prepared, random_state=seed)
    prepared.uns["independent_umap_contract"] = {
        "source": str(source_name),
        "input_observations": input_observations,
        "sampled_observations": int(prepared.n_obs),
        "input_genes": input_genes,
        "genes_after_hvg": int(prepared.n_vars),
        "max_observations": int(max_observations),
        "n_pcs": int(n_comps),
        "n_neighbors": int(n_neighbors),
        "seed": int(seed),
    }
    return prepared


def independent_umaps(
    sources: Mapping[str, ad.AnnData],
    *,
    categorical_labels: Mapping[str, Iterable[str] | pd.Series | np.ndarray],
    max_observations: int = MAX_OBSERVATIONS_PER_SOURCE,
    seed: int = SEED,
) -> dict[str, ad.AnnData]:
    if set(sources) != set(categorical_labels):
        raise ValueError("Each UMAP source must have one explicit display-label vector")
    return {
        name: independent_umap(
            obj,
            categorical_labels=categorical_labels[name],
            source_name=name,
            max_observations=max_observations,
            seed=seed,
        )
        for name, obj in sources.items()
    }


def plot_independent_panels(
    embeddings: Mapping[str, ad.AnnData],
    order: Iterable[str],
    *,
    titles: Mapping[str, str] | None = None,
    raw_source: str | None = None,
    figsize: tuple[float, float] | None = None,
):
    order = list(order)
    if set(order) != set(embeddings):
        raise ValueError("Independent UMAP panel order must list every embedding once")

    category_order: list[str] = []
    for source in order:
        if source == raw_source:
            continue
        labels = embeddings[source].obs["_display_label"].astype(str)
        for label in pd.unique(labels):
            if label not in category_order:
                category_order.append(label)
    palette = plt.get_cmap("tab20").colors
    colors = {
        label: palette[index % len(palette)]
        for index, label in enumerate(category_order)
    }

    fig, axes = plt.subplots(
        1,
        len(order),
        figsize=figsize or (5.2 * len(order), 4.5),
        squeeze=False,
    )
    axes = axes.ravel()
    for ax, source in zip(axes, order):
        embedding = embeddings[source]
        coords = np.asarray(embedding.obsm["X_umap"])
        if raw_source is not None and source == raw_source:
            values = np.log1p(
                embedding.obs["_total_counts"].to_numpy(dtype=float)
            )
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                cmap="viridis",
                s=4,
                alpha=0.75,
                linewidths=0,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02, label="log1p counts")
        else:
            labels = embedding.obs["_display_label"].astype(str)
            present = set(labels)
            for label in (label for label in category_order if label in present):
                selected = labels.to_numpy() == label
                ax.scatter(
                    coords[selected, 0],
                    coords[selected, 1],
                    s=4,
                    alpha=0.72,
                    color=colors[label],
                    label=label,
                    linewidths=0,
                    rasterized=True,
                )
            ax.legend(
                bbox_to_anchor=(1.01, 1),
                loc="upper left",
                frameon=False,
                fontsize=7,
                markerscale=2,
            )
        ax.set_title((titles or {}).get(source, source))
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    plt.tight_layout()
    return fig


def plot_joint_panels(
    embedding: ad.AnnData,
    order: Iterable[str],
    *,
    titles: Mapping[str, str] | None = None,
    raw_source: str | None = None,
    shared_limits: bool = True,
    figsize: tuple[float, float] | None = None,
):
    order = list(order)
    coords = np.asarray(embedding.obsm["X_umap"])
    fig, axes = plt.subplots(
        1,
        len(order),
        figsize=figsize or (5.2 * len(order), 4.5),
        squeeze=False,
    )
    axes = axes.ravel()
    palette = plt.get_cmap("tab20").colors
    x_span = float(np.ptp(coords[:, 0]))
    y_span = float(np.ptp(coords[:, 1]))
    x_pad = 0.03 * x_span if x_span else 1.0
    y_pad = 0.03 * y_span if y_span else 1.0
    x_limits = (float(coords[:, 0].min() - x_pad), float(coords[:, 0].max() + x_pad))
    y_limits = (float(coords[:, 1].min() - y_pad), float(coords[:, 1].max() + y_pad))
    for ax, source in zip(axes, order):
        mask = embedding.obs["_source"].astype(str).to_numpy() == source
        if raw_source is not None and source == raw_source:
            values = np.log1p(
                embedding.obs.loc[mask, "_total_counts"].to_numpy(dtype=float)
            )
            scatter = ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=values,
                cmap="viridis",
                s=4,
                alpha=0.75,
                linewidths=0,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02, label="log1p counts")
        else:
            labels = embedding.obs.loc[mask, "_display_label"].astype(str)
            categories = labels.value_counts().index.tolist()
            for index, label in enumerate(categories):
                selected = mask.copy()
                selected[mask] = labels.to_numpy() == label
                ax.scatter(
                    coords[selected, 0],
                    coords[selected, 1],
                    s=4,
                    alpha=0.72,
                    color=palette[index % len(palette)],
                    label=label,
                    linewidths=0,
                    rasterized=True,
                )
            ax.legend(
                bbox_to_anchor=(1.01, 1),
                loc="upper left",
                frameon=False,
                fontsize=7,
                markerscale=2,
            )
        ax.set_title((titles or {}).get(source, source))
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        if shared_limits:
            ax.set_xlim(x_limits)
            ax.set_ylim(y_limits)
    plt.tight_layout()
    return fig


def assigned_labels(adata: ad.AnnData) -> pd.Series:
    for column in ("cell_type", "Level1", "SVC_cluster", "Level2"):
        if column in adata.obs:
            return adata.obs[column].astype(str)
    raise KeyError("No reconstructed label column found")


def broad_osmfish_label(label: str) -> str:
    text = str(label)
    lowered = text.lower()
    for prefix, broad in (
        ("pyramidal", "Pyramidal"),
        ("inhibitory", "Inhibitory"),
        ("astrocyte", "Astrocyte"),
        ("oligodendrocyte", "Oligodendrocyte"),
        ("endothelial", "Endothelial"),
    ):
        if lowered.startswith(prefix):
            return broad
    return text
