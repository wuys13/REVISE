"""Lazy, read-only analysis views for a verified batch reconstruction."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from revise.application.ist_assembly import cluster_means


__all__ = [
    "AnalysisInputs",
    "InputView",
    "PairedView",
    "SSTBaseline",
    "load_analysis_inputs",
]


def _typed_key(value: Any) -> tuple[type, Any]:
    return type(value), value


def _axis_values(adata, axis: str, source: str) -> pd.Index:
    values = pd.Index(getattr(adata, axis))
    if values.empty or not values.is_unique or values.isna().any() or any(not str(value).strip() for value in values):
        raise ValueError(f"{source}.{axis} must be unique and nonempty")
    return values


def _coordinate_storage(adata, *, allow_obs: bool = False) -> str | None:
    if "spatial" in adata.obsm:
        return "obsm['spatial']"
    if allow_obs and {"x", "y"}.issubset(set(adata.obs.columns)):
        return "obs[x,y]"
    return None


def _coordinates(
    adata,
    source: str,
    *,
    allow_obs: bool = False,
    return_source: bool = False,
) -> np.ndarray | tuple[np.ndarray, str]:
    storage = _coordinate_storage(adata, allow_obs=allow_obs)
    if storage is None:
        raise ValueError(f"{source} is missing spatial coordinates")
    if storage == "obsm['spatial']":
        values = np.asarray(adata.obsm["spatial"])
    else:
        try:
            values = adata.obs.loc[:, ["x", "y"]].to_numpy(dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{source} spatial coordinates must be numeric") from exc
    if values.ndim != 2 or values.shape[0] != adata.n_obs or values.shape[1] < 2:
        raise ValueError(f"{source} spatial coordinates must have shape (n_obs, >=2)")
    if not np.isfinite(values).all():
        raise ValueError(f"{source} spatial coordinates must be finite")
    return (values, storage) if return_source else values


def _path_from_entry(entry: Any, label: str) -> Path:
    if isinstance(entry, Mapping):
        entry = entry.get("path")
    if not isinstance(entry, (str, Path)) or not str(entry).strip():
        raise ValueError(f"{label} must contain a path")
    return Path(entry)


def _source_path(record: Mapping[str, Any]) -> Path:
    inputs = record.get("inputs", {})
    sources = inputs.get("sources", {}) if isinstance(inputs, Mapping) else {}
    if isinstance(sources, Mapping) and "spatial" in sources:
        return _path_from_entry(sources["spatial"], "inputs.sources.spatial")
    raise ValueError("verified reconstruction is missing its raw spatial source")


def _output_path(record: Mapping[str, Any], role: str) -> Path:
    outputs = record.get("outputs", {})
    if not isinstance(outputs, Mapping) or role not in outputs:
        raise ValueError(f"verified reconstruction is missing output role {role!r}")
    return _path_from_entry(outputs[role], f"outputs.{role}")


def _labels(adata, source: str) -> list[tuple[type, Any]]:
    if "SVC_cluster" not in adata.obs:
        raise ValueError(f"{source} is missing SVC_cluster")
    labels = adata.obs["SVC_cluster"]
    if labels.isna().any() or any(not str(value).strip() for value in labels.to_numpy(dtype=object)):
        raise ValueError(f"{source} contains null or empty SVC_cluster labels")
    return [_typed_key(value) for value in labels.to_numpy(dtype=object)]


def _copy_axis_values(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        return values.copy()
    return np.asarray(values).copy()


@dataclass
class InputView:
    """One materialized analysis carrier, with explicit source and mapping."""

    adata: ad.AnnData
    mapping: str
    source: str
    provenance: dict[str, Any] = field(default_factory=dict)
    uncovered_observations: tuple[Any, ...] = ()
    coordinate_metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def X(self):
        return self.adata.X

    @property
    def observation_ids(self) -> pd.Index:
        return self.adata.obs_names

    @property
    def gene_ids(self) -> pd.Index:
        return self.adata.var_names

    @property
    def coordinates(self) -> np.ndarray | None:
        if "spatial" not in self.adata.obsm:
            if (self.mapping != "generated_spatial_units"
                    or self.provenance.get("coordinate_source") != "output.obs[x,y]"):
                return None
            return _coordinates(self.adata, self.source, allow_obs=True)
        return _coordinates(self.adata, self.source)

    @property
    def donor_ids(self) -> pd.Index:
        column = self.provenance.get("donor_column")
        if column is None or column not in self.adata.obs:
            return pd.Index([], dtype=object)
        return pd.Index(self.adata.obs[column].to_numpy(dtype=object))

    def close(self) -> None:
        """Close an accidentally backed carrier; normal views are in memory."""
        file_handle = getattr(self.adata, "file", None)
        if file_handle is not None:
            file_handle.close()


@dataclass
class SSTBaseline(InputView):
    """Equal-split raw baseline on the generated-cell axis."""

    parent_spot_ids: tuple[Any, ...] = ()
    covered_raw_spots: tuple[Any, ...] = ()
    uncovered_raw_spots: tuple[Any, ...] = ()


class PairedView:
    """Lazy native iST spatial/donor carriers and cluster helpers."""

    def __init__(self, inputs: "AnalysisInputs") -> None:
        self._inputs = inputs

    @property
    def spatial(self) -> InputView:
        return self._inputs._load_output("spatial", mapping="cluster")

    @property
    def donor(self) -> InputView:
        return self._inputs._load_output("expression", mapping="cluster")

    @property
    def raw(self) -> InputView:
        return self._inputs.aligned_raw()

    def _validate(self) -> tuple[list[tuple[type, Any]], list[tuple[type, Any]]]:
        spatial_keys = _labels(self.spatial.adata, "spatial carrier")
        donor_keys = _labels(self.donor.adata, "expression carrier")
        if not spatial_keys or set(spatial_keys) != set(donor_keys):
            raise ValueError("spatial and donor SVC cluster sets must match exactly and be nonempty")
        return spatial_keys, donor_keys

    @property
    def spatial_cluster_index(self) -> np.ndarray:
        spatial_keys, donor_keys = self._validate()
        order = {key: index for index, key in enumerate(dict.fromkeys(donor_keys))}
        return np.asarray([order[key] for key in spatial_keys], dtype=np.intp)

    def cluster_means(self):
        self._validate()
        means, values = cluster_means(self.donor.adata)
        return means, pd.Index(values)

    def map_to_spatial(self, cluster_values):
        indexer = self.spatial_cluster_index
        if getattr(cluster_values, "ndim", None) not in {1, 2}:
            raise ValueError("cluster values must be one- or two-dimensional")
        if cluster_values.shape[0] != len(set(_labels(self.donor.adata, "expression carrier"))):
            raise ValueError("cluster values do not match the donor cluster axis")
        return cluster_values[indexer].copy()


class AnalysisInputs:
    """Lazy access to raw and reconstruction carriers from a verified handoff.

    Constructing this object performs no H5AD reads. Each view is loaded only
    when its corresponding method or property is requested.
    """

    def __init__(self, reconstruction: Mapping[str, Any]) -> None:
        if not isinstance(reconstruction, Mapping):
            raise TypeError("reconstruction must be a mapping")
        if reconstruction.get("status") != "succeeded":
            raise ValueError("reconstruction handoff is not successful")
        self.reconstruction = reconstruction
        self.modality = reconstruction.get("modality")
        if self.modality not in {"hST", "iST", "sST"}:
            raise ValueError("reconstruction modality must be hST, iST, or sST")
        self.ist_mapping = reconstruction.get("ist_mapping", "paired") if self.modality == "iST" else None
        if self.modality == "iST" and self.ist_mapping not in {"paired", "mean", "random"}:
            raise ValueError("iST reconstruction mapping must be paired, mean, or random")
        self.coordinate_metadata = dict(reconstruction.get("coordinates", {}))
        self._raw_path = _source_path(reconstruction)
        self._cache: dict[str, InputView] = {}

    @classmethod
    def from_reconstruction(cls, reconstruction: Mapping[str, Any]) -> "AnalysisInputs":
        return cls(reconstruction)

    def raw(self) -> InputView:
        if "raw" not in self._cache:
            self._cache["raw"] = self._load(self._raw_path, role="raw", mapping="raw")
        return self._cache["raw"]

    def _output_role(self) -> str:
        if self.modality == "iST" and self.ist_mapping == "paired":
            return "spatial"
        return "svc"

    def _output_mapping(self) -> str:
        if self.modality == "iST":
            return self.ist_mapping
        return "native_id" if self.modality == "hST" else "generated_spatial_units"

    def _load(self, path: Path, *, role: str, mapping: str) -> InputView:
        data = ad.read_h5ad(path)
        _axis_values(data, "obs_names", role)
        _axis_values(data, "var_names", role)
        provenance = {}
        metadata = data.uns.get("revise_reconstruction", {})
        if isinstance(metadata, Mapping):
            provenance.update(dict(metadata))
        if role == "svc" and self.modality == "sST":
            storage = _coordinate_storage(data, allow_obs=True)
            if storage is not None:
                provenance["coordinate_source"] = f"output.{storage}"
                provenance["coordinate_semantics"] = (
                    "generated spatial-unit coordinates; may be parent-spot centers "
                    "or explicit cell centers"
                )
        if role == "svc" and self.modality == "iST" and self.ist_mapping == "random":
            column = "revise_ist_donor_id"
            if column not in data.obs:
                data.file.close()
                raise ValueError("random iST SVC output is missing donor provenance")
            values = data.obs[column]
            if values.isna().any() or any(not str(value).strip() for value in values.to_numpy(dtype=object)):
                data.file.close()
                raise ValueError("random iST donor provenance must be nonempty")
            provenance["donor_column"] = column
            provenance.setdefault("donor_id_source", "reference_cell")
        return InputView(data, mapping=mapping, source=str(path), provenance=provenance,
                         coordinate_metadata=dict(self.coordinate_metadata))

    def _load_output(self, role: str, *, mapping: str) -> InputView:
        cache_key = f"output:{role}:{mapping}"
        if cache_key not in self._cache:
            self._cache[cache_key] = self._load(_output_path(self.reconstruction, role), role=role,
                                                mapping=mapping)
        return self._cache[cache_key]

    def spatial(self) -> InputView:
        return self._load_output(self._output_role(), mapping=self._output_mapping())

    def reconstructed(self) -> InputView:
        if self.modality == "iST" and self.ist_mapping == "paired":
            raise ValueError("paired iST has no reconstructed expression carrier; use paired() donor/cluster means")
        return self.spatial()

    def reconstruction_expression(self, strategy: str = "auto") -> InputView:
        """Return the expression carrier permitted by this reconstruction route."""
        allowed = {"auto", "native", "cluster_mean_projection"}
        if not isinstance(strategy, str) or strategy not in allowed:
            raise ValueError(
                "strategy must be one of: auto, native, cluster_mean_projection"
            )
        if self.modality == "sST":
            raise ValueError("sST generated units do not provide reconstruction expression")

        paired = self.modality == "iST" and self.ist_mapping == "paired"
        if strategy == "cluster_mean_projection" and not paired:
            raise ValueError("strategy 'cluster_mean_projection' is only available for paired iST")
        if strategy == "native" and paired:
            raise ValueError("strategy 'native' is unavailable for paired iST")
        if not paired:
            return self.reconstructed()

        aligned = self.aligned_raw()
        paired_view = self.paired()
        spatial = paired_view.spatial
        donor = paired_view.donor
        means, cluster_labels = paired_view.cluster_means()
        projected_x = paired_view.map_to_spatial(means)

        donor_keys = _labels(donor.adata, "expression carrier")
        donor_cluster_keys = list(dict.fromkeys(donor_keys))
        donor_cells_per_cluster = [donor_keys.count(key) for key in donor_cluster_keys]

        outputs = self.reconstruction.get("outputs", {})
        source_identity = {}
        for role, view in (("spatial", spatial), ("expression", donor)):
            entry = outputs.get(role) if isinstance(outputs, Mapping) else None
            source_identity[role] = dict(entry) if isinstance(entry, Mapping) else {"path": view.source}
        inputs = self.reconstruction.get("inputs", {})
        sources = inputs.get("sources", {}) if isinstance(inputs, Mapping) else {}
        raw_entry = sources.get("spatial") if isinstance(sources, Mapping) else None
        source_identity["raw"] = (
            dict(raw_entry) if isinstance(raw_entry, Mapping) else {"path": aligned.source}
        )
        record_identity = {
            key: self.reconstruction[key]
            for key in ("fingerprint", "sample_id", "cell_type", "modality", "ist_mapping")
            if key in self.reconstruction
        }
        provenance = {
            "expression_view": "cluster_mean_projection",
            "expression_source": "expression_carrier.X_as_is",
            "donor_id_source": "reference_cell",
            "cluster_labels": list(cluster_labels),
            "donor_cells_per_cluster": donor_cells_per_cluster,
            "source_paths": {
                "raw": aligned.source,
                "spatial": spatial.source,
                "expression": donor.source,
            },
            "source_identity": source_identity,
            "record_identity": record_identity,
            "carrier_provenance": {
                "spatial": dict(spatial.provenance),
                "expression": dict(donor.provenance),
            },
            "uncovered_raw_observation_count": len(aligned.uncovered_observations),
        }
        expression_semantics = self.reconstruction.get("expression_semantics")
        if expression_semantics is not None:
            provenance["expression_semantics"] = expression_semantics

        projected = ad.AnnData(
            X=projected_x,
            obs=spatial.adata.obs.copy(deep=True),
            var=donor.adata.var.copy(deep=True),
            uns={"revise_analysis_input": dict(provenance)},
            obsm={key: _copy_axis_values(value) for key, value in spatial.adata.obsm.items()},
            obsp={key: _copy_axis_values(value) for key, value in spatial.adata.obsp.items()},
            varm={key: _copy_axis_values(value) for key, value in donor.adata.varm.items()},
            varp={key: _copy_axis_values(value) for key, value in donor.adata.varp.items()},
        )
        return InputView(
            projected,
            mapping="cluster_mean_projection",
            source=spatial.source,
            provenance=provenance,
            uncovered_observations=aligned.uncovered_observations,
            coordinate_metadata=dict(self.coordinate_metadata),
        )

    def paired(self) -> PairedView:
        if self.modality != "iST" or self.ist_mapping != "paired":
            raise ValueError("native paired view is unavailable; generated sST units require spot_name baseline")
        pairing = self.reconstruction.get("pairing", {})
        if isinstance(pairing, Mapping) and pairing.get("status") not in (None, "available"):
            raise ValueError(f"native paired view unavailable: {pairing.get('reason', 'unverified pairing')}")
        return PairedView(self)

    def aligned_raw(self) -> InputView:
        if self.modality == "sST":
            raise ValueError("sST generated units do not have native raw observation pairing")
        raw = self.raw()
        output = self.spatial()
        raw_ids = _axis_values(raw.adata, "obs_names", "raw")
        output_ids = _axis_values(output.adata, "obs_names", "reconstructed output")
        indexer = raw_ids.get_indexer(output_ids)
        if (indexer < 0).any():
            raise ValueError("reconstructed output contains raw observation IDs that are missing from raw")
        raw_xy = _coordinates(raw.adata, "raw")
        output_xy = _coordinates(output.adata, "reconstructed output")
        expected_xy = raw_xy[indexer]
        if output_xy.shape != expected_xy.shape or not np.array_equal(output_xy, expected_xy):
            raise ValueError("reconstructed output coordinates do not match raw observation IDs")
        covered = {_typed_key(value) for value in output_ids.to_numpy(dtype=object)}
        uncovered = tuple(value for value in raw_ids if _typed_key(value) not in covered)
        aligned = raw.adata[indexer, :].copy()
        return InputView(aligned, mapping="native_id", source=raw.source,
                         provenance={"output_source": output.source, "alignment": "observation_id"},
                         uncovered_observations=uncovered,
                         coordinate_metadata=dict(self.coordinate_metadata))

    def sst_baseline(self) -> SSTBaseline:
        if self.modality != "sST":
            raise ValueError("parent_spot_equal_split baseline is only available for sST")
        raw = self.raw()
        output = self.spatial()
        raw_ids = _axis_values(raw.adata, "obs_names", "raw")
        _axis_values(output.adata, "obs_names", "generated output")
        if "spot_name" not in output.adata.obs:
            raise ValueError("generated sST output is missing spot_name parent mapping")
        parent_values = output.adata.obs["spot_name"].to_numpy(dtype=object)
        if any(pd.isna(value) or not str(value).strip() for value in parent_values):
            raise ValueError("generated sST output spot_name parent mapping must be nonempty")
        parent_index = raw_ids.get_indexer(parent_values)
        if (parent_index < 0).any():
            raise ValueError("generated sST output contains an unknown raw spot_name")
        coordinates, coordinate_storage = _coordinates(
            output.adata, "generated output", allow_obs=True, return_source=True
        )
        counts = np.bincount(parent_index, minlength=raw.adata.n_obs)
        weights = 1.0 / counts[parent_index]
        parent_map = sparse.csr_matrix(
            (weights, (np.arange(output.adata.n_obs), parent_index)),
            shape=(output.adata.n_obs, raw.adata.n_obs),
        )
        raw_x = raw.adata.X
        if sparse.issparse(raw_x):
            baseline_x = (parent_map @ raw_x).tocsr()
        else:
            baseline_x = sparse.csr_matrix(parent_map @ np.asarray(raw_x, dtype=float))
        _check_sst_conservation(baseline_x, raw_x, parent_index)
        if not np.isfinite(baseline_x.data).all():
            raise ValueError("sST raw baseline must contain only finite values")
        baseline_obsm = {
            key: _copy_axis_values(value) for key, value in output.adata.obsm.items()
        }
        if "spatial" not in baseline_obsm:
            baseline_obsm["spatial"] = coordinates.copy()
        baseline = ad.AnnData(
            X=baseline_x,
            obs=output.adata.obs.copy(deep=True),
            var=raw.adata.var.copy(deep=True),
            uns={"revise_analysis_input": {
                "mapping": "parent_spot_equal_split",
                "parent_column": "spot_name",
                "raw_source": raw.source,
                "output_source": output.source,
            }},
            obsm=baseline_obsm,
        )
        parent_ids = tuple(parent_values.tolist())
        covered = tuple(dict.fromkeys(parent_ids))
        covered_keys = {_typed_key(value) for value in covered}
        uncovered = tuple(value for value in raw_ids if _typed_key(value) not in covered_keys)
        return SSTBaseline(
            baseline,
            mapping="parent_spot_equal_split",
            source=raw.source,
            provenance={
                "parent_column": "spot_name",
                "raw_source": raw.source,
                "output_source": output.source,
                "coordinate_source": f"output.{coordinate_storage}",
                "coordinate_semantics": (
                    "generated spatial-unit coordinates; may be parent-spot centers "
                    "or explicit cell centers"
                ),
            },
            uncovered_observations=uncovered,
            coordinate_metadata=dict(self.coordinate_metadata),
            parent_spot_ids=parent_ids,
            covered_raw_spots=covered,
            uncovered_raw_spots=uncovered,
        )

    def close(self) -> None:
        for view in self._cache.values():
            view.close()
        self._cache.clear()


def _check_sst_conservation(baseline_x, raw_x, parent_index: np.ndarray) -> None:
    n_generated = len(parent_index)
    n_raw = raw_x.shape[0]
    grouping = sparse.csr_matrix(
        (np.ones(n_generated, dtype=float), (parent_index, np.arange(n_generated))),
        shape=(n_raw, n_generated),
    )
    summed = grouping @ baseline_x
    covered_rows = np.unique(parent_index)
    if sparse.issparse(raw_x):
        delta = (summed[covered_rows] - raw_x[covered_rows]).tocsr()
        differences = delta.data
        scale_values = raw_x[covered_rows].data
        scale = max(float(np.max(np.abs(scale_values))) if len(scale_values) else 0.0, 1.0)
        if len(differences) and np.max(np.abs(differences)) > 1e-8 + 1e-6 * scale:
            raise ValueError("sST raw baseline does not conserve raw spot expression")
        return
    # Keep dense raw sources row-wise so a large source is never copied into a
    # second dense matrix merely for the conservation assertion.
    for raw_row in covered_rows:
        expected = np.asarray(raw_x[raw_row]).ravel()
        actual = np.asarray(summed.getrow(raw_row).toarray()).ravel()
        scale = max(float(np.max(np.abs(expected))) if len(expected) else 0.0, 1.0)
        if np.max(np.abs(actual - expected)) > 1e-8 + 1e-6 * scale:
            raise ValueError("sST raw baseline does not conserve raw spot expression")


def load_analysis_inputs(reconstruction: Mapping[str, Any]) -> AnalysisInputs:
    """Return a lazy helper for a verified reconstruction handoff."""

    return AnalysisInputs.from_reconstruction(reconstruction)
