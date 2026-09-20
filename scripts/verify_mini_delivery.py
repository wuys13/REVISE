"""Verify a real mini delivery against its immutable ROI input, without rewriting it."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def equal_matrix(left, right):
    assert left.shape == right.shape
    if sparse.issparse(left) or sparse.issparse(right):
        assert (sparse.csr_matrix(left) != sparse.csr_matrix(right)).nnz == 0
    else:
        np.testing.assert_array_equal(left, right)


def _dense_matrix(value):
    return value.toarray() if sparse.issparse(value) else np.asarray(value)


def _sst_true_zero_mask(diagnostic, *, active_parent_ids, gene_count, expected):
    if not isinstance(diagnostic, dict):
        raise AssertionError("sST parent-gene correction diagnostic must be a mapping")
    if diagnostic.get("schema_version") != 1:
        raise AssertionError("sST parent-gene correction diagnostic schema is unsupported")
    parent_ids = np.asarray(diagnostic.get("true_zero_support_parent_ids", ()), dtype=str)
    pairs = np.asarray(diagnostic.get("true_zero_support_pairs", ()), dtype=np.int64)
    if pairs.size == 0:
        pairs = np.empty((0, 2), dtype=np.int64)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise AssertionError("sST true-zero support pairs must have two columns")
    if (pairs[:, 0] < 0).any() or (pairs[:, 0] >= parent_ids.size).any():
        raise AssertionError("sST true-zero support parent index is invalid")
    if (pairs[:, 1] < 0).any() or (pairs[:, 1] >= gene_count).any():
        raise AssertionError("sST true-zero support gene index is invalid")
    if np.unique(pairs, axis=0).shape[0] != pairs.shape[0]:
        raise AssertionError("sST true-zero support pairs must be unique")
    parent_lookup = {str(parent_id): row for row, parent_id in enumerate(active_parent_ids)}
    if any(str(parent_id) not in parent_lookup for parent_id in parent_ids):
        raise AssertionError("sST true-zero support parent is absent from the delivery")
    mask = np.zeros_like(expected, dtype=bool)
    if pairs.size:
        rows = np.fromiter(
            (parent_lookup[str(parent_ids[index])] for index in pairs[:, 0]),
            dtype=np.intp,
            count=pairs.shape[0],
        )
        mask[rows, pairs[:, 1]] = True
    if np.any(mask & ~(expected > 0)):
        raise AssertionError("sST true-zero support must have a positive parent target")
    if diagnostic.get("true_zero_support_entries") != int(mask.sum()):
        raise AssertionError("sST true-zero support count disagrees with its sparse mask")
    if not np.isclose(
        diagnostic.get("true_zero_support_target_mass"), expected[mask].sum()
    ):
        raise AssertionError("sST true-zero support mass disagrees with Raw targets")
    return mask


def verify_sst_parent_gene_conservation(raw, svc) -> dict:
    """Compare virtual-cell totals with the runner's normalized parent targets.

    The sST runner first restricts Raw to the SVC gene axis, normalizes each
    parent to 10,000, then rescales each gene across that parent's virtual
    cells. A positive target with a strictly zero pre-rescale allocation
    cannot be restored by multiplication. Read that support from the saved
    diagnostic; a legacy artifact's final effective-zero values do not prove
    its pre-correction support.
    """
    parents = svc.obs["spot_name"].astype(str)
    if not parents.isin(raw.obs_names).all():
        raise AssertionError("Every virtual cell must name a Raw parent")
    gene_positions = raw.var_names.get_indexer(svc.var_names)
    if (gene_positions < 0).any():
        raise AssertionError("SVC genes must be present in Raw")
    parent_positions = raw.obs_names.get_indexer(parents)
    if (parent_positions < 0).any():
        raise AssertionError("Virtual-cell parents must align to Raw")

    raw_values = _dense_matrix(raw.X[:, gene_positions]).astype(np.float64, copy=False)
    raw_totals = raw_values.sum(axis=1)
    expected_by_raw = np.divide(
        raw_values * 1e4,
        raw_totals[:, None],
        out=np.zeros_like(raw_values),
        where=raw_totals[:, None] > 0,
    )
    svc_values = _dense_matrix(svc.X).astype(np.float64, copy=False)
    aggregated_by_raw = np.zeros_like(expected_by_raw)
    np.add.at(aggregated_by_raw, parent_positions, svc_values)
    active_parent_positions = pd.unique(parent_positions)
    active_parent_ids = raw.obs_names[active_parent_positions]
    expected = expected_by_raw[active_parent_positions]
    aggregated = aggregated_by_raw[active_parent_positions]

    positive_target = expected > 0
    effective_zero = positive_target & (aggregated <= 1e-12)
    residual = aggregated - expected
    parent_total_residual = aggregated.sum(axis=1) - expected.sum(axis=1)
    cell_totals = svc_values.sum(axis=1)
    expected_total = float(expected.sum())
    common = {
        "target": "Raw restricted to the SVC gene axis then per-parent normalize_total(10000)",
        "parent_count": int(parents.nunique()),
        "gene_count": int(svc.n_vars),
        "positive_target_entries": int(positive_target.sum()),
        "effective_zero_entries": int(effective_zero.sum()),
        "parent_total_max_absolute_error": float(np.abs(parent_total_residual).max()),
        "cell_total": {
            "min": float(cell_totals.min()),
            "median": float(np.median(cell_totals)),
            "max": float(cell_totals.max()),
        },
    }
    diagnostic = svc.uns.get("sst_parent_gene_correction")
    if diagnostic is None:
        unexpected_output = ~positive_target & (aggregated > 1e-12)
        close = np.isclose(aggregated, expected, rtol=1e-5, atol=1e-6)
        return {
            **common,
            "status": "failed" if unexpected_output.any() else "unobserved",
            "pre_correction_support_status": "unobserved",
            "effective_zero_target_mass": float(expected[effective_zero].sum()),
            "unclassified_positive_target_nonconserving_entries": int(
                (positive_target & ~close).sum()
            ),
            "effective_zero_examples": [
                {
                    "parent_id": str(active_parent_ids[row]),
                    "gene": str(svc.var_names[column]),
                    "normalized_parent_target": float(expected[row, column]),
                    "aggregated_svc": float(aggregated[row, column]),
                }
                for row, column in np.argwhere(effective_zero)[:10]
            ],
            "unexpected_positive_output_entries": int(unexpected_output.sum()),
        }

    true_zero_support = _sst_true_zero_mask(
        diagnostic,
        active_parent_ids=active_parent_ids,
        gene_count=svc.n_vars,
        expected=expected,
    )
    supported = positive_target & ~true_zero_support
    if diagnostic.get("positive_target_entries") != int(positive_target.sum()):
        raise AssertionError("sST positive target count disagrees with Raw targets")
    if diagnostic.get("positive_support_entries") != int(supported.sum()):
        raise AssertionError("sST positive support count disagrees with its sparse mask")
    close = np.isclose(aggregated, expected, rtol=1e-5, atol=1e-6)
    supported_mismatch = supported & ~close
    unexpected_output = ~positive_target & (aggregated > 1e-12)
    true_zero_nonzero_output = true_zero_support & (aggregated > 0)
    true_zero_target_mass = float(expected[true_zero_support].sum())
    status = (
        "failed"
        if unexpected_output.any() or supported_mismatch.any() or true_zero_nonzero_output.any()
        else "partial"
        if true_zero_support.any()
        else "passed"
    )
    return {
        **common,
        "status": status,
        "pre_correction_support_status": "observed",
        "correction_operator": diagnostic.get("operator"),
        "true_zero_support_entries": int(true_zero_support.sum()),
        "true_zero_support_target_mass": true_zero_target_mass,
        "true_zero_support_target_mass_fraction": (
            true_zero_target_mass / expected_total if expected_total else 0.0
        ),
        "true_zero_support_target_max": float(
            expected[true_zero_support].max() if true_zero_support.any() else 0.0
        ),
        "parent_count_with_true_zero_support": int(true_zero_support.any(axis=1).sum()),
        "parent_max_missing_mass_fraction": float(
            np.max(np.divide(expected * true_zero_support, expected.sum(axis=1)[:, None],
                             out=np.zeros_like(expected),
                             where=expected.sum(axis=1)[:, None] > 0).sum(axis=1))
        ),
        "true_zero_support_examples": [
            {
                "parent_id": str(active_parent_ids[row]),
                "gene": str(svc.var_names[column]),
                "normalized_parent_target": float(expected[row, column]),
                "aggregated_svc": float(aggregated[row, column]),
            }
            for row, column in np.argwhere(true_zero_support)[:10]
        ],
        "positive_support_entries": int(supported.sum()),
        "positive_support_nonconserving_entries": int(supported_mismatch.sum()),
        "positive_support_max_absolute_error": float(
            np.abs(residual[supported]).max() if supported.any() else 0.0
        ),
        "positive_support_max_relative_error": float(
            (np.abs(residual[supported]) / expected[supported]).max()
            if supported.any()
            else 0.0
        ),
        "true_zero_nonzero_output_entries": int(true_zero_nonzero_output.sum()),
        "unexpected_positive_output_entries": int(unexpected_output.sum()),
    }


def verify(run: Path, sample: str, mapping: str) -> dict:
    manifest = json.loads((run / "input_manifest.json").read_text())
    source_record = manifest["samples"][sample]
    input_path = ROOT / source_record["mini_input"]
    assert digest(input_path) == source_record["mini_sha256"]
    delivery = run / "delivery" / f"{sample}_mini" / mapping
    config_path = delivery / "sample.yaml"
    document = yaml.safe_load(config_path.read_text())
    raw_path = delivery / document["files"]["raw"]
    svc_path = delivery / document["files"]["svc"]
    hashes = {str(p): digest(p) for p in (input_path, raw_path, svc_path, config_path)}
    original, raw, svc = (ad.read_h5ad(p) for p in (input_path, raw_path, svc_path))
    assert original.obs_names.equals(raw.obs_names)
    pd.testing.assert_frame_equal(raw.obs.loc[:, original.obs.columns], original.obs)
    pd.testing.assert_frame_equal(raw.var, original.var)
    equal_matrix(raw.X, original.X)
    assert set(original.obsm) <= set(raw.obsm)
    for key in original.obsm:
        if isinstance(original.obsm[key], pd.DataFrame):
            pd.testing.assert_frame_equal(raw.obsm[key], original.obsm[key])
        else:
            equal_matrix(raw.obsm[key], original.obsm[key])
    for key in original.layers:
        equal_matrix(raw.layers[key], original.layers[key])
    assert svc.obs_names.is_unique and svc.var_names.is_unique
    assert document["files"] == {"raw": "raw.h5ad", "svc": "SVC.h5ad"}
    assert document["sample_id"] == f"{sample}_mini"
    if sample == "P2CRC_Visium":
        parents = svc.obs["spot_name"].astype(str)
        assert parents.isin(raw.obs_names).all()
        coordinates = raw.obsm["spatial"][raw.obs_names.get_indexer(parents)]
        assert "SVC_cluster" not in svc.obs
        relationship = "virtual-cell coordinates equal parent spot coordinates"
        parent_gene_conservation = verify_sst_parent_gene_conservation(raw, svc)
    else:
        assert svc.obs_names.isin(raw.obs_names).all()
        coordinates = raw.obsm["spatial"][raw.obs_names.get_indexer(svc.obs_names)]
        relationship = "retained native IDs and exact original coordinates"
        parent_gene_conservation = None
    np.testing.assert_array_equal(svc.obsm["spatial"], coordinates)
    for target in (raw, svc):
        for column in ("revise_Level1", "revise_Level2"):
            if column in target.obs:
                assert not target.obs[column].dropna().astype(str).str.contains("/", regex=False).any()
    for role, obj in (("raw", raw), ("svc", svc)):
        values = obj.X.data if sparse.issparse(obj.X) else np.asarray(obj.X)
        assert np.isfinite(values).all() and (values >= 0).all(), role
    assert hashes == {str(p): digest(p) for p in (input_path, raw_path, svc_path, config_path)}
    counts = lambda obj, column: ({str(k): int(v) for k, v in obj.obs[column].value_counts().items()}
                                  if column in obj.obs else {})
    status = parent_gene_conservation["status"] if parent_gene_conservation else "passed"
    record = {"kind": "real-data mini delivery verification", "sample": sample, "mapping": mapping,
              "status": status, "input_shape": list(original.shape), "raw_shape": list(raw.shape),
              "svc_shape": list(svc.shape), "raw_original_matrix_axes_labels_coordinates_preserved": True,
              "inputs_unchanged": True, "relationship": relationship,
              "raw_inferred_broad_counts": counts(raw, "revise_Level1"),
              "svc_broad_counts": counts(svc, "revise_Level1"),
              "raw_missing_subtype": int(raw.obs["revise_Level2"].isna().sum()) if "revise_Level2" in raw.obs else raw.n_obs,
              "expression": document["expression"], "coordinates": document["spatial"], "hashes": hashes}
    if parent_gene_conservation is not None:
        record["parent_gene_conservation"] = parent_gene_conservation
    path = run / "evidence" / f"delivery-{sample}-{mapping}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=ROOT / "output/mini-acceptance/20260920")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mapping", default="default")
    args = parser.parse_args()
    record = verify(args.run.resolve(), args.sample, args.mapping)
    print(json.dumps({k: v for k, v in record.items() if k != "hashes"}, indent=2))
    if record["status"] == "failed":
        raise SystemExit(1)
    if record["status"] == "unobserved":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
