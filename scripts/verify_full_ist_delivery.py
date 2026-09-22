#!/usr/bin/env python3
"""Audit one full iST random delivery without changing its three-file payload."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import traceback

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from revise.utils.provenance import hash_jsonable
from revise.utils.labels import normalize_cell_type_label


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray, pd.Index)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _equal_values(left, right) -> bool:
    if sparse.issparse(left) or sparse.issparse(right):
        left = sparse.csr_matrix(left)
        right = sparse.csr_matrix(right)
        return left.shape == right.shape and (left != right).nnz == 0
    return np.array_equal(np.asarray(left), np.asarray(right), equal_nan=True)


def _compare_matrix(left, right, *, label: str, chunk_rows: int = 4096) -> None:
    if tuple(left.shape) != tuple(right.shape):
        raise AssertionError(f"{label} shape changed: {left.shape} != {right.shape}")
    for start in range(0, left.shape[0], chunk_rows):
        stop = min(left.shape[0], start + chunk_rows)
        if not _equal_values(left[start:stop], right[start:stop]):
            raise AssertionError(f"{label} changed in rows {start}:{stop}")


def _validate_nonnegative_finite(matrix, *, label: str, chunk_rows: int = 4096) -> None:
    for start in range(0, matrix.shape[0], chunk_rows):
        stop = min(matrix.shape[0], start + chunk_rows)
        block = matrix[start:stop]
        values = block.data if sparse.issparse(block) else np.asarray(block)
        if not np.isfinite(values).all() or (values < 0).any():
            raise AssertionError(f"{label} contains negative or non-finite values in rows {start}:{stop}")


def _atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(_jsonable(value), indent=2, ensure_ascii=False) + "\n")
    os.replace(temporary, path)


def _write_donor_sidecar(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as stream:
        frame.to_csv(stream, sep="\t", index=False)
    os.replace(temporary, path)


def _scientific_config(document: dict) -> dict:
    normalized = copy.deepcopy(document)
    normalized["inputs"]["st"].pop("path")
    normalized["inputs"]["reference"].pop("path")
    normalized["output"].pop("dir")
    return normalized


def verify(
    *,
    repo: Path,
    config_path: Path,
    expected_st_sha256: str,
    expected_reference_sha256: str,
    audit_path: Path,
    donor_sidecar: Path,
) -> dict:
    document = yaml.safe_load(config_path.read_text())
    base_config_path = repo / "configs/acceptance/full/P2CRC_Xenium-random.yaml"
    base_document = yaml.safe_load(base_config_path.read_text())
    scientific_config = _scientific_config(document)
    if scientific_config != _scientific_config(base_document):
        raise AssertionError("run-scoped config changes frozen scientific parameters")
    if document.get("application") != {"svc_type": "sc-SVC", "mode": "cluster"}:
        raise AssertionError("this audit requires the full sc-SVC cluster route")
    if document.get("algorithm", {}).get("ot_method") != "tacco":
        raise AssertionError("full P2 Xenium audit requires TACCO")
    if document.get("output", {}).get("ist_mapping") != "random":
        raise AssertionError("full P2 Xenium audit requires output.ist_mapping=random")
    if document.get("execution", {}).get("seed") != 42:
        raise AssertionError("full P2 Xenium audit requires execution.seed=42")
    reference_request = document.get("inputs", {}).get("reference", {})
    if (reference_request.get("filter_column"), reference_request.get("filter_value")) != ("Patient", "P2CRC"):
        raise AssertionError("full P2 Xenium audit requires reference filter Patient=P2CRC")
    coordinate_request = document.get("delivery", {}).get("coordinates", {})
    if coordinate_request != {"key": "spatial", "unit": "pixel", "microns_per_coordinate": 0.2125}:
        raise AssertionError("full P2 Xenium audit requires the frozen Xenium coordinate declaration")
    root = repo if document["paths"]["root_dir"] == "." else Path(document["paths"]["root_dir"])
    st_path = (root / document["inputs"]["st"]["path"]).resolve()
    reference_path = (root / document["inputs"]["reference"]["path"]).resolve()
    output_dir = (root / document["output"]["dir"]).resolve()
    sample_path = output_dir / "sample.yaml"
    if not sample_path.is_file():
        raise AssertionError(f"missing delivery contract: {sample_path}")
    sample = yaml.safe_load(sample_path.read_text())
    if sample.get("files") != {"raw": "raw.h5ad", "svc": "SVC.h5ad"}:
        raise AssertionError("delivery must contain the canonical raw.h5ad and SVC.h5ad roles")
    raw_path = output_dir / sample["files"]["raw"]
    svc_path = output_dir / sample["files"]["svc"]
    if not raw_path.is_file() or not svc_path.is_file():
        raise AssertionError("delivery is missing raw.h5ad or SVC.h5ad")

    hashes = {
        "st_input": digest(st_path),
        "reference_input": digest(reference_path),
        "raw": digest(raw_path),
        "svc": digest(svc_path),
        "sample": digest(sample_path),
    }
    if hashes["st_input"] != expected_st_sha256:
        raise AssertionError("spatial input SHA256 does not match the frozen identity")
    if hashes["reference_input"] != expected_reference_sha256:
        raise AssertionError("reference input SHA256 does not match the frozen identity")

    original = ad.read_h5ad(st_path, backed="r")
    reference = ad.read_h5ad(reference_path, backed="r")
    raw = ad.read_h5ad(raw_path, backed="r")
    svc = ad.read_h5ad(svc_path, backed="r")
    try:
        if sample.get("sample_id") != "P2CRC_Xenium":
            raise AssertionError("delivery sample_id must remain P2CRC_Xenium")
        if sample.get("spatial") != coordinate_request:
            raise AssertionError("sample.yaml spatial declaration differs from the frozen config")
        if not original.obs_names.equals(raw.obs_names):
            raise AssertionError("published Raw observation IDs differ from the original input")
        if not original.var_names.equals(raw.var_names):
            raise AssertionError("published Raw genes differ from the original input")
        pd.testing.assert_frame_equal(
            raw.obs.loc[:, original.obs.columns], original.obs, check_like=False
        )
        pd.testing.assert_frame_equal(raw.var, original.var, check_like=False)
        _compare_matrix(original.X, raw.X, label="Raw.X")
        if set(original.layers) != set(raw.layers):
            raise AssertionError("published Raw layers differ from the original input")
        for key in original.layers:
            _compare_matrix(original.layers[key], raw.layers[key], label=f"Raw.layers[{key!r}]")
        if not set(original.obsm).issubset(raw.obsm):
            raise AssertionError("published Raw is missing original obsm entries")
        for key in original.obsm:
            if not _equal_values(original.obsm[key], raw.obsm[key]):
                raise AssertionError(f"Raw.obsm[{key!r}] changed")

        if not svc.obs_names.is_unique or not svc.var_names.is_unique:
            raise AssertionError("SVC observation and gene IDs must be unique")
        if svc.n_obs == 0 or svc.n_vars == 0:
            raise AssertionError("SVC must contain observations and genes")
        positions = raw.obs_names.get_indexer(svc.obs_names)
        if (positions < 0).any():
            raise AssertionError("SVC observation IDs are not a subset of complete Raw")
        spatial_key = sample["spatial"]["key"]
        if spatial_key not in raw.obsm or spatial_key not in svc.obsm:
            raise AssertionError(f"delivery is missing obsm[{spatial_key!r}]")
        if not _equal_values(svc.obsm[spatial_key], raw.obsm[spatial_key][positions]):
            raise AssertionError("SVC coordinates do not exactly match complete Raw by observation ID")
        _validate_nonnegative_finite(raw.X, label="Raw.X")
        _validate_nonnegative_finite(svc.X, label="SVC.X")

        donor_column = "revise_ist_donor_id"
        donor_cluster_column = "revise_ist_donor_cluster"
        for column in ("SVC_cluster", donor_column, donor_cluster_column):
            if column not in svc.obs:
                raise AssertionError(f"SVC is missing required random-donor provenance column {column!r}")
            if svc.obs[column].isna().any():
                raise AssertionError(f"SVC random-donor provenance column {column!r} contains null values")
        donor_ids = svc.obs[donor_column].astype(str)
        if donor_ids.str.strip().eq("").any():
            raise AssertionError("SVC random-donor provenance contains blank donor IDs")
        donor_clusters = svc.obs[donor_cluster_column].astype(str)
        svc_clusters = svc.obs["SVC_cluster"].astype(str)
        if not np.array_equal(donor_clusters.to_numpy(), svc_clusters.to_numpy()):
            raise AssertionError("sampled donor clusters do not match the spatial carrier clusters")
        donor_cluster_counts = (
            pd.DataFrame({"donor_id": donor_ids, "donor_cluster": donor_clusters})
            .drop_duplicates()
            .groupby("donor_id", observed=True)["donor_cluster"]
            .nunique()
        )
        if (donor_cluster_counts != 1).any():
            raise AssertionError("one donor ID maps to multiple expression-carrier clusters")

        sidecar = pd.DataFrame(
            {
                "output_id": svc.obs_names.astype(str),
                "donor_id": donor_ids.to_numpy(),
                "donor_cluster": donor_clusters.to_numpy(),
                "svc_cluster": svc_clusters.to_numpy(),
            }
        )
        _write_donor_sidecar(donor_sidecar, sidecar)
        reuse = donor_ids.value_counts()
        donor_frame = pd.DataFrame({"donor_id": donor_ids, "donor_cluster": donor_clusters})
        reuse_by_cluster = []
        for cluster, frame in donor_frame.groupby("donor_cluster", observed=True):
            cluster_reuse = frame["donor_id"].value_counts()
            assignments = int(len(frame))
            unique_donors = int(cluster_reuse.size)
            reuse_by_cluster.append(
                {
                    "donor_cluster": str(cluster),
                    "assignments": assignments,
                    "unique_donors": unique_donors,
                    "max_reuse": int(cluster_reuse.max()),
                    "repeated_assignment_fraction": (
                        (assignments - unique_donors) / assignments if assignments else 0.0
                    ),
                    "assignments_to_reused_donors_fraction": (
                        float(cluster_reuse[cluster_reuse > 1].sum()) / assignments
                        if assignments else 0.0
                    ),
                }
            )
        provenance = _jsonable(svc.uns.get("revise_reconstruction", {}))
        for key in ("ist_mapping", "effective_seed", "reference_filter", "processed_cell_types", "skipped_cell_types"):
            if key not in provenance:
                raise AssertionError(f"SVC reconstruction provenance is missing {key!r}")
        if provenance["ist_mapping"] != "random" or provenance["effective_seed"] != 42:
            raise AssertionError("actual SVC provenance does not record random mapping with seed 42")
        if provenance.get("donor_column") != donor_column:
            raise AssertionError("actual SVC provenance does not identify the donor ID column")
        if provenance.get("donor_cluster_column") != donor_cluster_column:
            raise AssertionError("actual SVC provenance does not identify the donor cluster column")
        if provenance.get("donor_sha256") != hash_jsonable(donor_ids.tolist()):
            raise AssertionError("actual donor ID hash disagrees with the delivered assignment")
        if provenance.get("donor_cluster_sha256") != hash_jsonable(donor_clusters.tolist()):
            raise AssertionError("actual donor cluster hash disagrees with the delivered assignment")
        if provenance["reference_filter"] != {"column": "Patient", "value": "P2CRC"}:
            raise AssertionError("actual SVC provenance does not record reference filter Patient=P2CRC")
        processed = provenance["processed_cell_types"]
        skipped = provenance["skipped_cell_types"]
        if not isinstance(processed, list) or not processed:
            raise AssertionError("processed_cell_types must be a nonempty explicit list")
        if not isinstance(skipped, (dict, list)):
            raise AssertionError("skipped_cell_types must be an explicit mapping or list")
        for label_column in ("revise_Level1", "revise_Level2"):
            if label_column not in raw.obs or label_column not in svc.obs:
                raise AssertionError(f"full delivery is missing required inferred label {label_column!r}")
        broad_counts = {
            "raw": {str(key): int(value) for key, value in raw.obs["revise_Level1"].value_counts(dropna=False).items()},
            "svc": {str(key): int(value) for key, value in svc.obs["revise_Level1"].value_counts(dropna=False).items()},
        }
        gene_intersections = {
            "raw_svc": int(raw.var_names.intersection(svc.var_names).size),
            "reference_svc": int(reference.var_names.intersection(svc.var_names).size),
            "raw_reference": int(raw.var_names.intersection(reference.var_names).size),
        }
        if gene_intersections["reference_svc"] != svc.n_vars:
            raise AssertionError("every SVC gene must be present in the reference")

        if not reference.obs_names.is_unique:
            raise AssertionError("reference observation IDs must be unique for donor provenance")
        donor_reference_positions = reference.obs_names.get_indexer(donor_ids)
        if (donor_reference_positions < 0).any():
            missing = donor_ids.iloc[np.flatnonzero(donor_reference_positions < 0)[:10]].tolist()
            raise AssertionError(f"sampled donor IDs are absent from the original reference: {missing}")
        donor_patients = reference.obs.iloc[donor_reference_positions]["Patient"].astype(str)
        if not donor_patients.eq("P2CRC").all():
            raise AssertionError("sampled donors include rows outside Patient=P2CRC")
        donor_reference_broad = (
            reference.obs.iloc[donor_reference_positions]["Level1"]
            .astype(str)
            .map(normalize_cell_type_label)
            .reset_index(drop=True)
        )
        donor_output_broad = (
            svc.obs["revise_Level1"].astype(str).reset_index(drop=True)
        )
        if not np.array_equal(donor_reference_broad.to_numpy(), donor_output_broad.to_numpy()):
            raise AssertionError("sampled donor broad labels disagree with normalized reference Level1")
        reference_gene_positions = reference.var_names.get_indexer(svc.var_names)
        if (reference_gene_positions < 0).any():
            raise AssertionError("SVC genes cannot be aligned to the original reference")
        for start in range(0, svc.n_obs, 1024):
            stop = min(svc.n_obs, start + 1024)
            block_ids = pd.Index(donor_ids.iloc[start:stop].to_numpy(), dtype="object")
            block_unique = pd.Index(pd.unique(block_ids), dtype="object")
            block_reference_positions = reference.obs_names.get_indexer(block_unique)
            reference_order = np.argsort(block_reference_positions)
            reference_unique_sorted = reference.X[block_reference_positions[reference_order]][
                :, reference_gene_positions
            ]
            reference_unique = reference_unique_sorted[np.argsort(reference_order)]
            reference_block = reference_unique[block_unique.get_indexer(block_ids)]
            svc_block = svc.X[start:stop]
            if not _equal_values(reference_block, svc_block):
                raise AssertionError(f"sampled donor expression differs from source reference rows {start}:{stop}")
        record = {
            "schema_version": 1,
            "kind": "full real-data iST random delivery audit",
            "status": "technical_verified",
            "scientific_acceptance": False,
            "config": str(config_path.relative_to(repo)),
            "scientific_config_sha256": hash_jsonable(scientific_config),
            "sample_id": sample.get("sample_id"),
            "input_shape": list(original.shape),
            "raw_shape": list(raw.shape),
            "svc_shape": list(svc.shape),
            "raw_original_matrix_axes_labels_coordinates_preserved": True,
            "svc_native_ids_and_coordinates_preserved": True,
            "processed_cell_types": processed,
            "skipped_cell_types": skipped,
            "broad_counts": broad_counts,
            "gene_intersections": gene_intersections,
            "expression": sample.get("expression"),
            "spatial": sample.get("spatial"),
            "random_donor_audit": {
                "assignment_count": int(len(donor_ids)),
                "unique_donor_count": int(donor_ids.nunique()),
                "source_reference_rows_verified": int(donor_ids.nunique()),
                "source_assignments_verified": int(len(donor_ids)),
                "source_expression_identity": "exact original-reference rows restricted to the SVC gene axis",
                "source_broad_identity": "normalized original-reference Level1 equals delivered revise_Level1",
                "donor_cluster_evidence": "expression-carrier assembly provenance; the original reference does not contain SVC_cluster",
                "repeated_assignment_fraction": float((len(donor_ids) - donor_ids.nunique()) / len(donor_ids)),
                "assignments_to_reused_donors_fraction": float(reuse[reuse > 1].sum() / len(donor_ids)),
                "reuse_min": int(reuse.min()),
                "reuse_median": float(reuse.median()),
                "reuse_max": int(reuse.max()),
                "reuse_quantiles": {
                    str(key): float(value)
                    for key, value in reuse.quantile([0.5, 0.9, 0.95, 0.99, 1.0]).items()
                },
                "by_cluster": reuse_by_cluster,
                "sidecar": str(donor_sidecar),
            },
            "files": {
                "spatial_input": {"path": str(st_path), "size": st_path.stat().st_size, "sha256": hashes["st_input"]},
                "reference_input": {"path": str(reference_path), "size": reference_path.stat().st_size, "sha256": hashes["reference_input"]},
                "raw": {"path": str(raw_path), "size": raw_path.stat().st_size, "sha256": hashes["raw"]},
                "svc": {"path": str(svc_path), "size": svc_path.stat().st_size, "sha256": hashes["svc"]},
                "sample": {"path": str(sample_path), "size": sample_path.stat().st_size, "sha256": hashes["sample"]},
                "donor_sidecar": {"path": str(donor_sidecar), "size": donor_sidecar.stat().st_size, "sha256": digest(donor_sidecar)},
            },
        }
    finally:
        original.file.close()
        reference.file.close()
        raw.file.close()
        svc.file.close()

    unchanged = {"st_input": digest(st_path), "reference_input": digest(reference_path)}
    if unchanged != {"st_input": hashes["st_input"], "reference_input": hashes["reference_input"]}:
        raise AssertionError("an input changed during delivery verification")
    _atomic_json(audit_path, record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--expected-st-sha256", required=True)
    parser.add_argument("--expected-reference-sha256", required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--donor-sidecar", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    audit = args.audit if args.audit.is_absolute() else repo / args.audit
    sidecar = args.donor_sidecar if args.donor_sidecar.is_absolute() else repo / args.donor_sidecar
    config = args.config if args.config.is_absolute() else repo / args.config
    try:
        record = verify(
            repo=repo,
            config_path=config.resolve(),
            expected_st_sha256=args.expected_st_sha256,
            expected_reference_sha256=args.expected_reference_sha256,
            audit_path=audit.resolve(),
            donor_sidecar=sidecar.resolve(),
        )
    except Exception as exc:
        _atomic_json(
            audit.resolve(),
            {
                "schema_version": 1,
                "kind": "full real-data iST random delivery audit",
                "status": "failed",
                "scientific_acceptance": False,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            },
        )
        raise
    print(json.dumps({key: value for key, value in record.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
