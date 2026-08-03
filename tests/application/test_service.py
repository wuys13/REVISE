"""Application result assembly and publication contract."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy import sparse

from revise.svc import SVC


ROOT = Path(__file__).resolve().parents[2]


class PublicationContext(SimpleNamespace):
    def set_pending_publication(self, *, commit, rollback):
        self._publication_commit = commit
        self._publication_rollback = rollback

    def rollback_pending_publication(self):
        rollback = getattr(self, "_publication_rollback", None)
        self._publication_commit = None
        self._publication_rollback = None
        if rollback is not None:
            rollback()


def _ist_outputs(*, sparse_x: bool = False):
    spatial = AnnData(
        X=np.zeros((3, 1)),
        obs=pd.DataFrame(
            {"SVC_cluster": ["A", "B", "A"], "position_owner": [1, 2, 3]},
            index=["position-1", "position-2", "position-3"],
        ),
        var=pd.DataFrame(index=["spatial-placeholder"]),
    )
    spatial.obsm["spatial"] = np.array(
        [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
    )
    spatial.obsm["position_embedding"] = np.arange(6).reshape(3, 2)
    spatial.obsp["position_graph"] = sparse.eye(3, format="csr")
    spatial.uns["revise_reconstruction"] = {"spatial_note": "keep"}
    spatial.layers["spatial_layer"] = np.ones((3, 1))
    spatial.raw = spatial.copy()

    values = np.array(
        [[2.0, 4.0], [6.0, 8.0], [10.0, 12.0], [14.0, 16.0]]
    )
    expression = AnnData(
        X=sparse.csr_matrix(values) if sparse_x else values,
        obs=pd.DataFrame(
            {"SVC_cluster": ["A", "A", "B", "B"]},
            index=["cell-2", "cell-1", "cell-4", "cell-3"],
        ),
        var=pd.DataFrame(
            {"gene_owner": ["x", "y"]},
            index=["g1", "g2"],
        ),
    )
    expression.varm["gene_embedding"] = np.arange(4).reshape(2, 2)
    expression.varp["gene_graph"] = sparse.eye(2, format="csr")
    expression.layers["expression_layer"] = expression.X.copy()
    expression.raw = expression.copy()
    return spatial, expression


def _context(tmp_path, spatial, expression, *, seed=17):
    records = []
    return PublicationContext(
        svc=SVC(
            expr=expression,
            spatial=spatial,
            svc_kind="sc",
            artifacts={
                "outputs": {
                    "sc_svc_spatial": spatial,
                    "sc_svc_expr": expression,
                }
            },
        ),
        run_dir=tmp_path / "run",
        runtime={"seed": seed},
        merged_config={
            "runtime": {"seed": seed},
            "ot": {"ga": {"solver": "tacco"}, "lr": {"solver": "tacco"}},
        },
        provenance={},
        artifact_records=records,
        record_artifact=records.append,
    )


def _args(tmp_path, *, mapping="mean", seed=None):
    return SimpleNamespace(
        svc_type="iST-SVC",
        ist_mapping=mapping,
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=seed,
        ot_method=None,
    )


def _publish(tmp_path, *, sparse_x=False, mapping="mean", seed=None):
    from revise.application.service import _build_public_result

    spatial, expression = _ist_outputs(sparse_x=sparse_x)
    ctx = _context(tmp_path, spatial, expression)
    result, path = _build_public_result(
        _args(tmp_path, mapping=mapping, seed=seed),
        "application_sc",
        None,
        ctx,
    )
    return result, path, ctx, spatial, expression


@pytest.mark.parametrize("sparse_x", (False, True))
def test_ist_mean_assembles_single_result_with_explicit_carrier_ownership(
    tmp_path,
    sparse_x,
):
    result, path, ctx, spatial, expression = _publish(
        tmp_path,
        sparse_x=sparse_x,
    )

    assert path == tmp_path / "out" / "sample" / "iST-SVC" / "SVC.h5ad"
    assert list(path.parent.glob("*.h5ad")) == [path]
    assert result.obs_names.equals(spatial.obs_names)
    assert result.obs.equals(spatial.obs)
    assert result.var_names.equals(expression.var_names)
    assert result.var.equals(expression.var)
    assert np.allclose(
        result.X.toarray() if sparse.issparse(result.X) else result.X,
        [[4.0, 6.0], [12.0, 14.0], [4.0, 6.0]],
    )
    assert np.array_equal(result.obsm["spatial"], spatial.obsm["spatial"])
    assert np.array_equal(
        result.obsm["position_embedding"], spatial.obsm["position_embedding"]
    )
    assert np.array_equal(
        result.varm["gene_embedding"], expression.varm["gene_embedding"]
    )
    assert set(result.obsp) == {"position_graph"}
    assert set(result.varp) == {"gene_graph"}
    assert not result.layers
    assert result.raw is None
    assert "revise_ist_donor_id" not in result.obs
    assert spatial.obs.columns.tolist() == ["SVC_cluster", "position_owner"]
    assert expression.obs.columns.tolist() == ["SVC_cluster"]
    assert ctx.provenance["result"]["filename"] == "SVC.h5ad"


def test_ist_random_is_invariant_to_expression_row_order_and_records_donors(
    tmp_path,
):
    from revise.application.service import _build_public_result

    spatial, expression = _ist_outputs()
    first, _, _, _, _ = _publish(tmp_path / "first", mapping="random", seed=731)
    reordered = expression[["cell-3", "cell-1", "cell-4", "cell-2"], :].copy()
    second, _ = _build_public_result(
        _args(tmp_path / "second", mapping="random", seed=731),
        "application_sc",
        None,
        _context(tmp_path / "second", spatial, reordered),
    )

    assert first.obs["revise_ist_donor_id"].tolist() == second.obs[
        "revise_ist_donor_id"
    ].tolist()
    assert np.array_equal(first.X, second.X)
    expected_json = json.dumps(
        first.obs["revise_ist_donor_id"].tolist(),
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    metadata = first.uns["revise_reconstruction"]
    assert metadata["effective_seed"] == 731
    assert metadata["donor_column"] == "revise_ist_donor_id"
    assert metadata["donor_sha256"] == hashlib.sha256(expected_json).hexdigest()


def test_ist_random_donor_selection_is_stable_across_fresh_processes():
    probe = textwrap.dedent(
        """
        import json
        from types import SimpleNamespace

        import numpy as np

        from revise.application.service import _random_expression

        expression = SimpleNamespace(
            X=np.array([[2, 4], [6, 8], [10, 12], [14, 16]]),
            obs_names=["cell-2", "cell-1", "cell-4", "cell-3"],
        )
        spatial_keys = [(str, "A"), (str, "B"), (str, "A")]
        expression_keys = [
            (str, "A"),
            (str, "A"),
            (str, "B"),
            (str, "B"),
        ]
        values, donors = _random_expression(
            expression, spatial_keys, expression_keys, 731
        )
        print(json.dumps({"donors": donors, "values": values.tolist()}))
        """
    )
    outputs = []
    for hash_seed in ("1", "987654"):
        env = os.environ.copy()
        env["PYTHONHASHSEED"] = hash_seed
        completed = subprocess.run(
            [sys.executable, "-c", probe],
            cwd=ROOT,
            env=env,
            text=True,
            capture_output=True,
            check=True,
        )
        outputs.append(json.loads(completed.stdout))

    assert outputs[0] == outputs[1]


@pytest.mark.parametrize("mapping", ("mean", "random"))
def test_ist_metadata_uses_the_minimal_2x_contract(tmp_path, mapping):
    result, _, ctx, _, _ = _publish(tmp_path, mapping=mapping)

    metadata = result.uns["revise_reconstruction"]
    assert metadata == {
        "spatial_note": "keep",
        "schema_version": 2,
        "svc_type": "iST-SVC",
        "ist_mapping": mapping,
        "effective_seed": 17 if mapping == "random" else None,
        "expression_source": "expression_carrier.X_as_is",
        "donor_column": "revise_ist_donor_id" if mapping == "random" else None,
        "donor_sha256": metadata["donor_sha256"] if mapping == "random" else None,
    }
    assert ctx.provenance["result"] == {
        "filename": "SVC.h5ad",
        "type": "iST-SVC",
    }
    assembly = ctx.provenance["assembly"]
    if mapping == "random":
        assert assembly == {
            "ist_mapping": "random",
            "effective_seed": 17,
            "donor_column": "revise_ist_donor_id",
            "donor_sha256": metadata["donor_sha256"],
            "donor_count": result.n_obs,
        }
    else:
        assert assembly == {
            "ist_mapping": "mean",
            "effective_seed": None,
            "donor_column": None,
            "donor_sha256": None,
            "donor_count": None,
        }


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda spatial, expression: expression.obs.__setitem__(
                "SVC_cluster", ["A", "A", "A", "A"]
            ),
            "cluster sets",
        ),
        (
            lambda spatial, expression: spatial.obs.__setitem__(
                "SVC_cluster", ["A", None, "A"]
            ),
            "null",
        ),
        (
            lambda spatial, expression: expression.obs.__setitem__(
                "SVC_cluster", [1, 1, "B", "B"]
            ),
            "cluster sets",
        ),
        (
            lambda spatial, expression: expression.__setattr__(
                "var_names", ["g1", "g1"]
            ),
            "var_names.*unique",
        ),
        (
            lambda spatial, expression: expression.__setattr__(
                "obs_names", ["cell-1", "cell-1", "cell-3", "cell-4"]
            ),
            "obs_names.*unique",
        ),
        (
            lambda spatial, expression: expression.__setattr__(
                "obs_names", ["cell-1", "", "cell-3", "cell-4"]
            ),
            "obs_names.*non-empty",
        ),
        (
            lambda spatial, expression: expression.X.__setitem__((0, 0), np.inf),
            "finite",
        ),
    ],
)
def test_ist_validation_fails_before_replacing_the_public_result(
    tmp_path,
    mutate,
    message,
):
    from revise.application.service import _build_public_result

    spatial, expression = _ist_outputs()
    mutate(spatial, expression)
    path = tmp_path / "out" / "sample" / "iST-SVC" / "SVC.h5ad"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"previous-result")

    with pytest.raises((KeyError, ValueError), match=message):
        _build_public_result(
            _args(tmp_path, mapping="random"),
            "application_sc",
            None,
            _context(tmp_path, spatial, expression),
        )

    assert path.read_bytes() == b"previous-result"


def test_ist_missing_carrier_fails_before_replacing_the_public_result(tmp_path):
    from revise.application.service import _build_public_result

    spatial, expression = _ist_outputs()
    ctx = _context(tmp_path, spatial, expression)
    del ctx.svc.artifacts["outputs"]["sc_svc_expr"]
    path = tmp_path / "out" / "sample" / "iST-SVC" / "SVC.h5ad"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"previous-result")

    with pytest.raises(RuntimeError, match="sc_svc_expr"):
        _build_public_result(
            _args(tmp_path),
            "application_sc",
            None,
            ctx,
        )

    assert path.read_bytes() == b"previous-result"


def test_contract_metadata_conflict_fails_without_silent_overwrite(tmp_path):
    from revise.application.service import _build_public_result

    spatial, expression = _ist_outputs()
    spatial.uns["revise_reconstruction"]["schema_version"] = 1
    path = tmp_path / "out" / "sample" / "iST-SVC" / "SVC.h5ad"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"previous-result")

    with pytest.raises(ValueError, match="schema_version"):
        _build_public_result(
            _args(tmp_path),
            "application_sc",
            None,
            _context(tmp_path, spatial, expression),
        )

    assert path.read_bytes() == b"previous-result"


def test_sc_svc_publication_requires_selected_cell_type(tmp_path):
    from revise.application.service import _build_sc_public_results

    spatial, expression = _sc_outputs()
    ctx = PublicationContext(
        svc=SVC(
            expr=expression,
            spatial=spatial,
            svc_kind="sc",
            artifacts={
                "outputs": {
                    "sc_svc_spatial": spatial,
                    "sc_svc_expr": expression,
                }
            },
        ),
        run_dir=tmp_path / "run",
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        artifact_records=[],
        record_artifact=lambda artifact: None,
    )
    args = SimpleNamespace(
        svc_type="sc-SVC",
        output_root=str(tmp_path / "out"),
        sample_name="sample",
        seed=17,
        ot_method="pot",
        select_ct=None,
    )

    with pytest.raises(RuntimeError, match="missing its selected cell type"):
        _build_sc_public_results(args, "application_sc", ctx)


def test_single_file_publication_reloads_staged_h5ad_before_replace(
    monkeypatch,
    tmp_path,
):
    from revise.application import service

    spatial, expression = _ist_outputs()
    path = tmp_path / "out" / "sample" / "iST-SVC" / "SVC.h5ad"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"previous-result")

    def corrupt_write(self, target, *args, **kwargs):
        Path(target).write_bytes(b"not-an-h5ad")

    monkeypatch.setattr(AnnData, "write_h5ad", corrupt_write)

    with pytest.raises(OSError):
        service._build_public_result(
            _args(tmp_path),
            "application_sc",
            None,
            _context(tmp_path, spatial, expression),
        )

    assert path.read_bytes() == b"previous-result"
    assert list(path.parent.iterdir()) == [path]


def test_written_result_matches_the_returned_result(tmp_path):
    result, path, _, _, _ = _publish(tmp_path, mapping="random")

    written = read_h5ad(path)
    assert written.obs_names.equals(result.obs_names)
    assert written.var_names.equals(result.var_names)
    assert np.array_equal(written.X, result.X)
    assert written.uns["revise_reconstruction"] == result.uns[
        "revise_reconstruction"
    ]
