"""Bounded Raw-versus-reconstruction pathway activity analysis."""

from __future__ import annotations

from collections.abc import Mapping
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .advanced import aucell as _aucell_module
from .advanced.aucell import get_aucell_provider_metadata, score_gene_set_aucell


CODE_DEPENDENCIES = ("revise.analysis.advanced.aucell",)

_REQUIRED_PARAMETERS = (
    "resource",
    "min_coverage",
    "min_genes",
    "min_detected_genes",
    "seed",
    "cutoff_policy",
    "auc_threshold_quantile",
    "level1_column",
)
_STATUS_COMPUTED = "computed"
_VALID_SCIENCE_STATUSES = {
    "computed",
    "unmeasured",
    "low_coverage",
    "insufficient_support",
    "not_computable",
}
_PROVIDER_DETECTED_QUANTILES = (0.01, 0.05, 0.10, 0.50, 1.0)


def _json_safe(value: Any) -> Any:
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, pd.DataFrame):
        return _json_safe(value.to_dict(orient="records"))
    if isinstance(value, pd.Series):
        return _json_safe(value.to_dict())
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (pd.Timestamp, pd.Timedelta)):
        return value.isoformat()
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    serialised = json.dumps(
        _json_safe(payload),
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
        allow_nan=False,
    )
    path.write_text(serialised + "\n", encoding="utf-8")


def _write_csv(path: Path, frame: pd.DataFrame, *, compressed: bool = False) -> None:
    work = frame.copy()
    work.columns = [str(column) for column in work.columns]
    text = work.to_csv(
        index=False,
        na_rep="",
        float_format="%.17g",
        lineterminator="\n",
    )
    if not compressed:
        path.write_text(text, encoding="utf-8")
        return
    with path.open("wb") as raw_stream:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_stream, mtime=0
        ) as compressed_stream:
            compressed_stream.write(text.encode("utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_parameters(parameters: Any) -> dict[str, Any]:
    if not isinstance(parameters, Mapping):
        raise ValueError("pathway activity parameters must be a mapping")
    missing = [
        name
        for name in _REQUIRED_PARAMETERS
        if name != "resource" and name not in parameters
    ]
    if "resource" not in parameters and not all(
        f"resource_{name}" in parameters
        for name in ("name", "species", "gene_id_type", "version")
    ):
        missing.insert(0, "resource")
    if missing:
        raise ValueError("pathway activity requires explicit parameters: " + ", ".join(missing))

    resource = parameters.get("resource")
    if not isinstance(resource, Mapping):
        # Accept the equivalent flat spelling without supplying any value.
        flat_names = ("name", "species", "gene_id_type", "version")
        if all(f"resource_{name}" in parameters for name in flat_names):
            resource = {
                name: parameters[f"resource_{name}"] for name in flat_names
            }
        else:
            raise ValueError("resource must contain species, gene_id_type, and version")
    resource = dict(resource)
    for name in ("name", "species", "gene_id_type", "version"):
        value = resource.get(name)
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"resource.{name} must be an explicit nonempty string")

    min_coverage = parameters["min_coverage"]
    if isinstance(min_coverage, bool) or not isinstance(min_coverage, (int, float)):
        raise ValueError("min_coverage must be a number in (0, 1]")
    if not 0 < float(min_coverage) <= 1:
        raise ValueError("min_coverage must be a number in (0, 1]")

    min_genes = parameters["min_genes"]
    if isinstance(min_genes, bool) or not isinstance(min_genes, int) or min_genes < 1:
        raise ValueError("min_genes must be a positive integer")

    min_detected_genes = parameters["min_detected_genes"]
    if (
        isinstance(min_detected_genes, bool)
        or not isinstance(min_detected_genes, int)
        or min_detected_genes < 0
    ):
        raise ValueError("min_detected_genes must be a nonnegative integer")

    seed = parameters["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an explicit integer")

    if parameters["cutoff_policy"] != "provider_detected_gene_quantile":
        raise ValueError(
            "cutoff_policy must be 'provider_detected_gene_quantile'"
        )
    quantile = parameters["auc_threshold_quantile"]
    if isinstance(quantile, bool) or not isinstance(quantile, (int, float)):
        raise ValueError("auc_threshold_quantile must be a number in [0, 1]")
    if float(quantile) not in _PROVIDER_DETECTED_QUANTILES:
        allowed = ", ".join(str(value) for value in _PROVIDER_DETECTED_QUANTILES)
        raise ValueError(
            "auc_threshold_quantile must be one of the provider quantiles: " + allowed
        )

    level1_column = parameters["level1_column"]
    if not isinstance(level1_column, str) or not level1_column.strip():
        raise ValueError("level1_column must be an explicit nonempty string")

    normalised = dict(parameters)
    normalised["resource"] = resource
    normalised["min_coverage"] = float(min_coverage)
    normalised["auc_threshold_quantile"] = float(quantile)
    return normalised


def _load_gene_sets(path: Path) -> dict[str, list[str]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"gene_sets resource is not valid JSON: {path}") from exc
    if not isinstance(payload, Mapping) or not payload:
        raise ValueError("gene_sets resource must be a nonempty pathway-to-gene mapping")
    gene_sets: dict[str, list[str]] = {}
    for pathway, genes in payload.items():
        if not isinstance(pathway, str) or not pathway.strip():
            raise ValueError("gene_sets pathway names must be nonempty strings")
        if not isinstance(genes, list):
            raise ValueError(f"gene_sets[{pathway!r}] must be a list")
        if any(
            not isinstance(gene, str) or not gene.strip() or gene != gene.strip()
            for gene in genes
        ):
            raise ValueError(f"gene_sets[{pathway!r}] contains an invalid gene ID")
        if len(set(genes)) != len(genes):
            raise ValueError(f"gene_sets[{pathway!r}] must contain exact unique genes")
        gene_sets[pathway] = list(genes)
    return gene_sets


def _detected_count(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        return np.asarray((values != 0).sum(axis=1)).ravel().astype(np.int64, copy=False)
    return np.count_nonzero(np.asarray(values), axis=1).astype(np.int64, copy=False)


def _signature_detected_observations(values: Any, indices: list[int]) -> int:
    if not indices:
        return 0
    if sparse.issparse(values):
        row_counts = np.asarray((values[:, indices] != 0).sum(axis=1)).ravel()
    else:
        row_counts = np.count_nonzero(np.asarray(values)[:, indices], axis=1)
    return int(np.count_nonzero(row_counts))


def _signature_detected_gene_count(values: Any, indices: list[int]) -> int:
    if not indices:
        return 0
    if sparse.issparse(values):
        return int(np.asarray((values[:, indices] != 0).sum(axis=0)).ravel().astype(bool).sum())
    return int(np.count_nonzero(np.any(np.asarray(values)[:, indices] != 0, axis=0)))


def _stored_matrix_bytes(values: Any) -> int:
    if sparse.issparse(values):
        return int(
            values.data.nbytes
            + values.indices.nbytes
            + values.indptr.nbytes
        )
    return int(np.asarray(values).nbytes)


def _view_profile(view: Any, quantile: float) -> dict[str, Any]:
    adata = getattr(view, "adata", view)
    var_names = list(adata.var_names)
    if len(set(var_names)) != len(var_names) or any(
        not isinstance(gene, str) or not gene.strip() for gene in var_names
    ):
        raise ValueError("expression gene axis must contain unique nonempty IDs")
    if int(adata.n_obs) < 1 or int(adata.n_vars) < 1:
        raise ValueError("expression view must contain observations and genes")
    detected = _detected_count(adata.X)
    count_quantile = float(pd.Series(detected).quantile(quantile))
    return {
        "adata": adata,
        "var_names": var_names,
        "var_index": {gene: index for index, gene in enumerate(var_names)},
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "detected_counts": detected,
        "detected_count_quantile": count_quantile,
        "n_vars_cutoff": int(adata.n_vars),
        "auc_rank_cutoff": count_quantile / float(adata.n_vars),
        "all_expression_zero": bool(np.count_nonzero(detected) == 0),
        "stored_matrix_bytes": _stored_matrix_bytes(adata.X),
        "estimated_dense_matrix_bytes": int(adata.n_obs * adata.n_vars * 8),
    }


def _side_results(
    view: Any,
    profile: dict[str, Any],
    gene_sets: Mapping[str, list[str]],
    parameters: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    values = profile["adata"].X
    for pathway, genes in gene_sets.items():
        indices = [profile["var_index"][gene] for gene in genes if gene in profile["var_index"]]
        available_genes = [gene for gene in genes if gene in profile["var_index"]]
        coverage = len(available_genes) / float(len(genes)) if genes else 0.0
        reasons = []
        if not genes:
            reasons.append("resource_has_no_genes")
        if coverage < float(parameters["min_coverage"]):
            reasons.append("coverage_below_minimum")
        if len(available_genes) < int(parameters["min_genes"]):
            reasons.append("min_genes_below_minimum")
        result: dict[str, Any] = {
            "pathway": pathway,
            "resource_genes": list(genes),
            "available_genes": available_genes,
            "available_gene_count": len(available_genes),
            "resource_gene_count": len(genes),
            "coverage": coverage,
            "detectable_observations": 0,
            "detected_signature_gene_count": 0,
            "min_detected_genes": int(parameters["min_detected_genes"]),
            "status": None,
            "reason": None,
            "values": None,
        }
        result["detectable_observations"] = _signature_detected_observations(values, indices)
        detected_signature_gene_count = _signature_detected_gene_count(values, indices)
        result["detected_signature_gene_count"] = detected_signature_gene_count
        if not available_genes and genes:
            result["status"] = "unmeasured"
            result["reason"] = "resource_genes_unmeasured"
        elif reasons:
            result["status"] = "low_coverage"
            result["reason"] = ";".join(reasons)
        elif profile["all_expression_zero"]:
            result["status"] = "not_computable"
            result["reason"] = "zero_expression_no_detectable_signature"
        elif profile["auc_rank_cutoff"] <= 0:
            result["status"] = "not_computable"
            result["reason"] = "zero_detected_gene_rank_cutoff"
        elif detected_signature_gene_count < int(parameters["min_detected_genes"]):
            result["status"] = "insufficient_support"
            result["reason"] = "insufficient_detected_signature_genes"
        else:
            # Provider errors intentionally propagate as execution errors.
            scored, score_key = score_gene_set_aucell(
                profile["adata"],
                available_genes,
                score_name=pathway,
                AUC_threshold=float(parameters["auc_threshold_quantile"]),
                seed=int(parameters["seed"]),
            )
            scores = np.asarray(scored.obs[score_key], dtype=float)
            if scores.shape != (profile["n_obs"],) or not np.isfinite(scores).all():
                raise RuntimeError(
                    f"AUCell provider returned invalid scores for pathway {pathway!r}"
                )
            result["status"] = _STATUS_COMPUTED
            result["values"] = scores
        if result["status"] not in _VALID_SCIENCE_STATUSES:
            raise RuntimeError(f"invalid pathway status for {pathway!r}")
        results[pathway] = result
    return results


def _coordinates(view: Any) -> np.ndarray | None:
    coordinates = getattr(view, "coordinates", None)
    if coordinates is None:
        adata = getattr(view, "adata", view)
        if "spatial" in adata.obsm:
            coordinates = adata.obsm["spatial"]
    if coordinates is None:
        return None
    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 2 or values.shape[1] < 2 or not np.isfinite(values[:, :2]).all():
        raise ValueError("pathway analysis coordinates must be finite two-dimensional values")
    return values[:, :2]


def _labels(raw_view: Any, level1_column: str, observation_ids: list[str]) -> list[str]:
    raw_adata = getattr(raw_view, "adata", raw_view)
    if level1_column not in raw_adata.obs:
        raise ValueError(f"aligned Raw view is missing {level1_column!r}")
    labels = raw_adata.obs[level1_column].tolist()
    if any(pd.isna(label) or not str(label).strip() for label in labels):
        raise ValueError(f"aligned Raw view contains empty {level1_column!r} labels")
    if len(labels) != len(observation_ids):
        raise ValueError("Raw Level1 labels do not match the observation axis")
    return [str(label) for label in labels]


def _validate_parent_labels(
    parameters: Mapping[str, Any], reconstruction: Mapping[str, Any], scopes: list[str]
) -> dict[str, str]:
    """Validate the shared task-to-Raw-Level1 mapping without normalizing labels."""
    mapping = parameters.get("parent_labels")
    if mapping is None:
        if reconstruction.get("modality") == "iST":
            raise ValueError("iST pathway analysis requires explicit parent_labels")
        return {}
    if not isinstance(mapping, Mapping) or not mapping:
        raise ValueError("parent_labels must be a nonempty task-to-Raw-label mapping")
    validated: dict[str, str] = {}
    for task_id, raw_label in mapping.items():
        if (
            not isinstance(task_id, str)
            or not task_id.strip()
            or task_id != task_id.strip()
            or task_id in {".", ".."}
            or "/" in task_id
            or "\\" in task_id
            or any(ord(char) < 32 for char in task_id)
        ):
            raise ValueError("parent_labels keys must be safe task IDs")
        if not isinstance(raw_label, str) or not raw_label.strip():
            raise ValueError("parent_labels values must be exact nonempty Raw labels")
        validated[task_id] = raw_label
    if len(set(validated.values())) != len(validated):
        raise ValueError("parent_labels must not map multiple tasks to one raw label")
    task_id = reconstruction.get("cell_type")
    if reconstruction.get("modality") == "iST":
        if not isinstance(task_id, str) or task_id not in validated:
            raise ValueError("parent_labels must contain the exact iST task cell_type")
        if validated[task_id] not in scopes:
            raise ValueError(
                f"parent_labels[{task_id!r}] does not match an exact Raw Level1 label"
            )
    return validated


def _comparison_rows(
    context: Any,
    raw_view: Any,
    reconstruction_view: Any,
    scopes: list[str],
    level1_column: str,
) -> list[dict[str, Any]]:
    reconstruction = getattr(context, "reconstruction", {})
    sample_id = str(reconstruction.get("sample_id", ""))
    task_cell_type = str(reconstruction.get("cell_type", ""))
    raw_mapping = str(getattr(raw_view, "mapping", "aligned_raw"))
    reconstruction_mapping = str(
        getattr(reconstruction_view, "mapping", "reconstruction_expression")
    )
    provenance = getattr(reconstruction_view, "provenance", {})
    normalization = reconstruction.get("expression_semantics")
    if normalization is None and isinstance(provenance, Mapping):
        normalization = provenance.get("expression_semantics")
    normalization = "unspecified" if normalization is None else str(normalization)
    rows = []
    for scope in scopes:
        rows.append(
            {
                "comparison_id": json.dumps(
                    {
                        "sample_id": sample_id,
                        "task_cell_type": task_cell_type,
                        "raw_level1": scope,
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=False,
                ),
                "sample_id": sample_id,
                "task_cell_type": task_cell_type,
                "scope": scope,
                "comparison_edge": "raw_vs_reconstruction",
                "raw_view": raw_mapping,
                "reconstruction_view": reconstruction_mapping,
                "label_source": f"raw.{level1_column}",
                "observation_basis": "spatial_unit",
                "gene_rule": "independent_full_gene_space",
                "normalization": normalization,
            }
        )
    return rows


def _summary(values: list[float]) -> dict[str, Any]:
    if not values:
        return {
            "n_valid": 0,
            "mean": None,
            "median": None,
            "q1": None,
            "q3": None,
        }
    series = pd.Series(values, dtype=float)
    return {
        "n_valid": int(series.size),
        "mean": float(series.mean()),
        "median": float(series.median()),
        "q1": float(series.quantile(0.25)),
        "q3": float(series.quantile(0.75)),
    }


def _input_tables(
    raw_view: Any,
    reconstruction_view: Any,
    labels: list[str],
    comparisons: list[dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_adata = getattr(raw_view, "adata", raw_view)
    reconstruction_adata = getattr(reconstruction_view, "adata", reconstruction_view)
    observation_ids = [str(value) for value in raw_adata.obs_names]
    raw_coordinates = _coordinates(raw_view)
    reconstruction_coordinates = _coordinates(reconstruction_view)
    if raw_coordinates is not None and raw_coordinates.shape[0] != len(observation_ids):
        raise ValueError("Raw coordinates do not match the observation axis")
    if reconstruction_coordinates is not None and reconstruction_coordinates.shape[0] != len(observation_ids):
        raise ValueError("reconstruction coordinates do not match the observation axis")
    comparison_by_scope = {row["scope"]: row["comparison_id"] for row in comparisons}
    observations = []
    for index, (unit_id, scope) in enumerate(zip(observation_ids, labels)):
        row = {
            "observation_order": index,
            "unit_id": unit_id,
            "comparison_id": comparison_by_scope[scope],
            "raw_level1": scope,
            "included_raw": True,
            "included_reconstruction": True,
            "raw_x": None,
            "raw_y": None,
            "reconstruction_x": None,
            "reconstruction_y": None,
        }
        if raw_coordinates is not None:
            row["raw_x"], row["raw_y"] = raw_coordinates[index]
        if reconstruction_coordinates is not None:
            row["reconstruction_x"], row["reconstruction_y"] = reconstruction_coordinates[index]
        observations.append(row)

    raw_genes = [str(value) for value in raw_adata.var_names]
    reconstruction_genes = [str(value) for value in reconstruction_adata.var_names]
    genes = list(dict.fromkeys(raw_genes + reconstruction_genes))
    genes_table = pd.DataFrame(
        [
            {
                "gene_id": gene,
                "raw_provided": gene in raw_genes,
                "reconstruction_provided": gene in reconstruction_genes,
                "feature_role": "expression",
            }
            for gene in genes
        ]
    )

    provenance = getattr(reconstruction_view, "provenance", {})
    cluster_labels = []
    donor_counts = {}
    if isinstance(provenance, Mapping):
        cluster_labels = [str(value) for value in provenance.get("cluster_labels", [])]
        counts = provenance.get("donor_cells_per_cluster", [])
        donor_counts = {
            label: int(count)
            for label, count in zip(cluster_labels, counts)
        }
    if "SVC_cluster" in reconstruction_adata.obs:
        unit_clusters = [str(value) for value in reconstruction_adata.obs["SVC_cluster"].tolist()]
    else:
        unit_clusters = [None] * len(observation_ids)
    mapping_name = str(getattr(reconstruction_view, "mapping", "reconstruction_expression"))
    mapping_table = pd.DataFrame(
        [
            {
                "unit_id": unit_id,
                "comparison_id": comparison_by_scope[scope],
                "raw_level1": scope,
                "mapping": mapping_name,
                "cluster_label": cluster,
                "donor_cells_in_cluster": donor_counts.get(cluster) if cluster is not None else None,
                "mapping_source": str(getattr(reconstruction_view, "source", "")),
            }
            for unit_id, scope, cluster in zip(observation_ids, labels, unit_clusters)
        ]
    )
    return pd.DataFrame(observations), genes_table, mapping_table


def _side_cutoff(profile: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "n_observations": int(profile["n_obs"]),
        "n_vars": int(profile["n_vars"]),
        "detected_count_quantile": float(profile["detected_count_quantile"]),
        "detected_count_cutoff": float(profile["detected_count_quantile"]),
        "n_vars_cutoff": int(profile["n_vars_cutoff"]),
        "auc_rank_cutoff": float(profile["auc_rank_cutoff"]),
        "provider_detected_count_quantile": float(profile["detected_count_quantile"]),
        "provider_n_vars": int(profile["n_vars"]),
        "provider_auc_threshold": float(profile["auc_rank_cutoff"]),
        "all_expression_zero": bool(profile["all_expression_zero"]),
        "stored_matrix_bytes": int(profile["stored_matrix_bytes"]),
        "estimated_dense_matrix_bytes": int(profile["estimated_dense_matrix_bytes"]),
    }


def _source_identity(view: Any) -> dict[str, Any]:
    source = str(getattr(view, "source", ""))
    identity: dict[str, Any] = {"path": source}
    path = Path(source) if source else None
    if path is not None and path.is_file():
        identity["sha256"] = _sha256(path)
    return identity


def run(context: Any) -> dict[str, Any]:
    """Run the bounded P4 pathway activity aspect for one verified task."""
    parameters = _validate_parameters(getattr(context, "parameters", None))
    resources = getattr(context, "resources", None)
    if not isinstance(resources, Mapping) or "gene_sets" not in resources:
        raise ValueError("pathway activity requires context.resources['gene_sets']")
    resource_path = Path(resources["gene_sets"]).resolve()
    if not resource_path.is_file():
        raise FileNotFoundError(f"gene_sets resource does not exist: {resource_path}")
    gene_sets = _load_gene_sets(resource_path)

    # Resolve the provider before science statuses are emitted. Provider failures
    # are execution failures and must not be represented as missing Raw biology.
    provider = get_aucell_provider_metadata()
    raw_view = context.inputs.aligned_raw()
    reconstruction = getattr(context, "reconstruction", {})
    expression_strategy = (
        "cluster_mean_projection"
        if reconstruction.get("modality") == "iST"
        and reconstruction.get("ist_mapping", "paired") == "paired"
        else "native"
    )
    reconstruction_view = context.inputs.reconstruction_expression(
        strategy=expression_strategy
    )
    raw_adata = getattr(raw_view, "adata", raw_view)
    reconstruction_adata = getattr(reconstruction_view, "adata", reconstruction_view)
    raw_axis = list(raw_adata.obs_names)
    reconstruction_axis = list(reconstruction_adata.obs_names)
    if raw_axis != reconstruction_axis:
        raise ValueError("aligned Raw and reconstruction expression axes must match exactly")
    raw_ids = [str(value) for value in raw_axis]
    if len(set(raw_ids)) != len(raw_ids):
        raise ValueError("observation IDs must remain unique after serialization")

    labels = _labels(raw_view, parameters["level1_column"], raw_ids)
    scopes = list(dict.fromkeys(labels))
    parent_labels = _validate_parent_labels(
        parameters,
        getattr(context, "reconstruction", {}),
        scopes,
    )
    comparisons = _comparison_rows(
        context,
        raw_view,
        reconstruction_view,
        scopes,
        parameters["level1_column"],
    )
    raw_profile = _view_profile(raw_view, parameters["auc_threshold_quantile"])
    reconstruction_profile = _view_profile(
        reconstruction_view, parameters["auc_threshold_quantile"]
    )
    raw_results = _side_results(raw_view, raw_profile, gene_sets, parameters)
    reconstruction_results = _side_results(
        reconstruction_view, reconstruction_profile, gene_sets, parameters
    )

    output_dir = Path(context.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = output_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    observations, genes, mapping = _input_tables(
        raw_view, reconstruction_view, labels, comparisons
    )

    comparison_by_scope = {row["scope"]: row for row in comparisons}
    scope_indices = {
        scope: [index for index, label in enumerate(labels) if label == scope]
        for scope in scopes
    }
    score_rows = []
    availability_rows = []
    summary_rows = []
    resource_gene_rows = []
    for pathway, resource_genes in gene_sets.items():
        raw_result = raw_results[pathway]
        reconstruction_result = reconstruction_results[pathway]
        raw_available = set(raw_result["available_genes"])
        reconstruction_available = set(reconstruction_result["available_genes"])
        raw_used = raw_result["status"] == _STATUS_COMPUTED
        reconstruction_used = reconstruction_result["status"] == _STATUS_COMPUTED
        for order, gene in enumerate(resource_genes):
            resource_gene_rows.append(
                {
                    "pathway": pathway,
                    "resource_gene_order": order,
                    "resource_gene": gene,
                    "raw_available": gene in raw_available,
                    "reconstruction_available": gene in reconstruction_available,
                    "raw_used": raw_used and gene in raw_available,
                    "reconstruction_used": reconstruction_used and gene in reconstruction_available,
                }
            )
        for scope in scopes:
            comparison = comparison_by_scope[scope]
            comparison_id = comparison["comparison_id"]
            indices = scope_indices[scope]
            raw_values = raw_result["values"]
            reconstruction_values = reconstruction_result["values"]
            paired = []
            for index in indices:
                raw_score = (
                    float(raw_values[index])
                    if raw_values is not None
                    else None
                )
                reconstruction_score = (
                    float(reconstruction_values[index])
                    if reconstruction_values is not None
                    else None
                )
                delta = (
                    reconstruction_score - raw_score
                    if raw_score is not None and reconstruction_score is not None
                    else None
                )
                if delta is not None:
                    paired.append(delta)
                score_rows.append(
                    {
                        "comparison_id": comparison_id,
                        "sample_id": comparison["sample_id"],
                        "task_cell_type": comparison["task_cell_type"],
                        "scope": scope,
                        "unit_id": raw_ids[index],
                        "pathway": pathway,
                        "raw_score": raw_score,
                        "reconstruction_score": reconstruction_score,
                        "paired_delta": delta,
                        "raw_status": raw_result["status"],
                        "reconstruction_status": reconstruction_result["status"],
                        "raw_reason": raw_result["reason"],
                        "reconstruction_reason": reconstruction_result["reason"],
                        "comparison_status": "computed" if delta is not None else "unavailable",
                    }
                )

            raw_summary = _summary(
                [float(raw_values[index]) for index in indices]
                if raw_values is not None
                else []
            )
            reconstruction_summary = _summary(
                [float(reconstruction_values[index]) for index in indices]
                if reconstruction_values is not None
                else []
            )
            delta_summary = _summary(paired)
            comparison_status = "computed" if paired else "unavailable"
            summary_rows.append(
                {
                    "comparison_id": comparison_id,
                    "sample_id": comparison["sample_id"],
                    "task_cell_type": comparison["task_cell_type"],
                    "scope": scope,
                    "pathway": pathway,
                    "raw_status": raw_result["status"],
                    "reconstruction_status": reconstruction_result["status"],
                    "comparison_status": comparison_status,
                    "raw_reason": raw_result["reason"],
                    "reconstruction_reason": reconstruction_result["reason"],
                    "raw_n_valid": raw_summary["n_valid"],
                    "raw_mean": raw_summary["mean"],
                    "raw_median": raw_summary["median"],
                    "raw_q1": raw_summary["q1"],
                    "raw_q3": raw_summary["q3"],
                    "reconstruction_n_valid": reconstruction_summary["n_valid"],
                    "reconstruction_mean": reconstruction_summary["mean"],
                    "reconstruction_median": reconstruction_summary["median"],
                    "reconstruction_q1": reconstruction_summary["q1"],
                    "reconstruction_q3": reconstruction_summary["q3"],
                    "paired_n": delta_summary["n_valid"],
                    "paired_delta_mean": delta_summary["mean"],
                    "paired_delta_median": delta_summary["median"],
                    "resource_gene_count": raw_result["resource_gene_count"],
                    "raw_available_gene_count": raw_result["available_gene_count"],
                    "reconstruction_available_gene_count": reconstruction_result["available_gene_count"],
                    "raw_used_gene_count": len(raw_available) if raw_used else 0,
                    "reconstruction_used_gene_count": (
                        len(reconstruction_available) if reconstruction_used else 0
                    ),
                    "raw_coverage": raw_result["coverage"],
                    "reconstruction_coverage": reconstruction_result["coverage"],
                }
            )
            for side, side_result, profile in (
                ("raw", raw_result, raw_profile),
                ("reconstruction", reconstruction_result, reconstruction_profile),
            ):
                availability_rows.append(
                    {
                        "comparison_id": comparison_id,
                        "scope": scope,
                        "object_type": "pathway",
                        "object_id": pathway,
                        "pathway": pathway,
                        "side": side,
                        "status": side_result["status"],
                        "reason": side_result["reason"],
                        "resource_gene_count": side_result["resource_gene_count"],
                        "available_gene_count": side_result["available_gene_count"],
                        "used_gene_count": (
                            side_result["available_gene_count"]
                            if side_result["status"] == _STATUS_COMPUTED
                            else 0
                        ),
                        "coverage": side_result["coverage"],
                        "detectable_observations": side_result["detectable_observations"],
                        "detected_signature_gene_count": side_result["detected_signature_gene_count"],
                        "min_detected_genes": side_result["min_detected_genes"],
                        "n_observations": profile["n_obs"],
                        "n_vars": profile["n_vars"],
                        "detected_count_quantile": profile["detected_count_quantile"],
                        "n_vars_cutoff": profile["n_vars_cutoff"],
                        "auc_rank_cutoff": profile["auc_rank_cutoff"],
                    }
                )

    resource_identity = {
        "path": str(resource_path),
        "sha256": _sha256(resource_path),
    }
    resource_payload = {
        "resource": dict(parameters["resource"]),
        **dict(parameters["resource"]),
        "path": str(resource_path),
        "sha256": resource_identity["sha256"],
        "provider": provider,
        "provider_name": provider["name"],
        "provider_version": provider["version"],
        "scorer": "AUCell",
        "cutoff_policy": parameters["cutoff_policy"],
        "auc_threshold_quantile": parameters["auc_threshold_quantile"],
        "seed": parameters["seed"],
        "min_coverage": parameters["min_coverage"],
        "min_genes": parameters["min_genes"],
        "min_detected_genes": parameters["min_detected_genes"],
        "pathways": list(gene_sets),
        "side_cutoffs": {
            "raw": _side_cutoff(raw_profile),
            "reconstruction": _side_cutoff(reconstruction_profile),
        },
    }
    code_paths = {
        str(Path(__file__).resolve()): _sha256(Path(__file__).resolve()),
    }
    aucell_path = Path(_aucell_module.__file__).resolve()
    if aucell_path.is_file():
        code_paths[str(aucell_path)] = _sha256(aucell_path)
    audit = {
        "schema_version": 1,
        "aspect": "pathway_activity",
        "method": "AUCell",
        "provider": provider,
        "code": code_paths,
        "resource": resource_identity,
        "resource_metadata": dict(parameters["resource"]),
        "parameters": parameters,
        "parent_labels": parent_labels,
        "input_views": {
            "raw": {
                "mapping": str(getattr(raw_view, "mapping", "aligned_raw")),
                "source": _source_identity(raw_view),
                "n_observations": raw_profile["n_obs"],
                "n_vars": raw_profile["n_vars"],
            },
            "reconstruction": {
                "mapping": str(getattr(reconstruction_view, "mapping", "reconstruction_expression")),
                "source": _source_identity(reconstruction_view),
                "provenance": getattr(reconstruction_view, "provenance", {}),
                "n_observations": reconstruction_profile["n_obs"],
                "n_vars": reconstruction_profile["n_vars"],
            },
            "observation_axis": raw_ids,
            "raw_level1_column": parameters["level1_column"],
            "raw_level1_labels": scopes,
        },
        "side_cutoffs": {
            "raw": _side_cutoff(raw_profile),
            "reconstruction": _side_cutoff(reconstruction_profile),
        },
        "output_files": {
            "audit": "audit.json",
            "comparisons": "comparisons.csv",
            "observations": "inputs/observations.csv.gz",
            "genes": "inputs/genes.csv",
            "mapping": "inputs/mapping.csv.gz",
            "availability": "availability.csv",
            "scores": "scores.csv.gz",
            "summary": "summary.csv",
            "resource": "resource.json",
            "resource_genes": "resource_genes.csv",
        },
    }

    _write_json(output_dir / "audit.json", audit)
    _write_csv(output_dir / "comparisons.csv", pd.DataFrame(comparisons))
    _write_csv(inputs_dir / "observations.csv.gz", observations, compressed=True)
    _write_csv(inputs_dir / "genes.csv", genes, compressed=False)
    _write_csv(inputs_dir / "mapping.csv.gz", mapping, compressed=True)
    _write_csv(output_dir / "availability.csv", pd.DataFrame(availability_rows))
    _write_csv(output_dir / "scores.csv.gz", pd.DataFrame(score_rows), compressed=True)
    _write_csv(output_dir / "summary.csv", pd.DataFrame(summary_rows))
    _write_json(output_dir / "resource.json", resource_payload)
    _write_csv(
        output_dir / "resource_genes.csv",
        pd.DataFrame(
            resource_gene_rows,
            columns=[
                "pathway",
                "resource_gene_order",
                "resource_gene",
                "raw_available",
                "reconstruction_available",
                "raw_used",
                "reconstruction_used",
            ],
        ),
    )

    artifacts = {
        "audit": {"path": "audit.json", "description": "Pathway analysis provenance and parameters"},
        "comparisons": {"path": "comparisons.csv", "description": "Raw Level1 comparison definitions"},
        "observations": {"path": "inputs/observations.csv.gz", "description": "Stable spatial observation axis"},
        "genes": {"path": "inputs/genes.csv", "description": "Independent Raw and reconstruction gene axes"},
        "mapping": {"path": "inputs/mapping.csv.gz", "description": "Expression mapping and projection provenance"},
        "availability": {"path": "availability.csv", "description": "Per-side pathway availability states"},
        "scores": {"path": "scores.csv.gz", "description": "Unit-level Raw and reconstruction pathway scores"},
        "summary": {"path": "summary.csv", "description": "Raw Level1 pathway score summaries"},
        "resource": {"path": "resource.json", "description": "Resource and provider identity"},
        "resource_genes": {"path": "resource_genes.csv", "description": "Resource genes and side availability"},
    }
    return {
        "artifacts": artifacts,
        "calculation": {
            "input_view": "aligned_raw + reconstruction_expression",
            "parameters": parameters,
            "comparison_basis": "independent_full_gene_space_aucell_with_aligned_unit_delta",
        },
    }
