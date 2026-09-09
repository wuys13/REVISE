from __future__ import annotations

import builtins
import sys
from types import ModuleType, SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from revise.analysis.advanced.cell_communication import (
    normalize_cellphonedb_tables,
    run_cellphonedb_v5,
)


PAIR_COLUMNS = ["A|A", "A|B", "B|A", "B|B"]


def _table(rows: int, pair_columns=PAIR_COLUMNS) -> pd.DataFrame:
    data: dict[str, object] = {
        "gene_a": [f"L{i}" for i in range(rows)],
        "gene_b": [f"R{i}" for i in range(rows)],
    }
    for offset, column in enumerate(pair_columns):
        data[column] = np.arange(rows, dtype=float) + offset
    return pd.DataFrame(data)


def _results() -> dict[str, pd.DataFrame]:
    return {
        "pvalues": _table(2),
        "significant_means": _table(1),
        "interaction_scores": _table(3),
        "deconvoluted_percents": pd.DataFrame(),
    }


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"cluster": ["A", "B"]}, index=["c1", "c2"]),
    )


def test_normalize_cellphonedb_tables_preserves_pair_order_and_unequal_row_counts():
    results = _results()

    normalized = normalize_cellphonedb_tables(results, ["A", "B"])

    assert list(normalized) == [
        "pvalues",
        "significant_means",
        "interaction_scores",
    ]
    assert [len(table) for table in normalized.values()] == [2, 1, 3]
    assert all(table.columns.tolist() == PAIR_COLUMNS for table in normalized.values())
    assert normalized["pvalues"].index.tolist() == ["L0_R0", "L1_R1"]
    assert "LR_pair" not in results["pvalues"].columns
    assert results["pvalues"].index.tolist() == [0, 1]


def test_normalize_cellphonedb_tables_rejects_mismatched_pair_column_sequence():
    results = _results()
    results["interaction_scores"] = _table(2, list(reversed(PAIR_COLUMNS)))

    with pytest.raises(ValueError, match="pair-column sequence"):
        normalize_cellphonedb_tables(results, ["A", "B"])


def test_normalize_cellphonedb_tables_accepts_valid_empty_tables():
    empty = _table(0)
    results = {
        "pvalues": empty,
        "significant_means": empty.copy(),
        "interaction_scores": empty.copy(),
    }

    normalized = normalize_cellphonedb_tables(results, ["A", "B"])

    assert all(table.empty for table in normalized.values())
    assert all(table.columns.tolist() == PAIR_COLUMNS for table in normalized.values())


def test_normalize_cellphonedb_tables_accepts_surviving_group_pairs():
    surviving_pairs = ["A|A"]
    results = {
        "pvalues": _table(2, surviving_pairs),
        "significant_means": _table(1, surviving_pairs),
        "interaction_scores": _table(3, surviving_pairs),
    }

    normalized = normalize_cellphonedb_tables(results, ["A", "B"])

    assert all(table.columns.tolist() == surviving_pairs for table in normalized.values())


def test_run_cellphonedb_v5_forwards_explicit_parameters_and_validates_results(
    monkeypatch, tmp_path
):
    calls = []
    results = _results()
    returned_adata = _adata()

    def fake_run(adata, **kwargs):
        calls.append((adata, kwargs))
        return results, returned_adata

    omicverse = ModuleType("omicverse")
    omicverse.single = SimpleNamespace(run_cellphonedb_v5=fake_run)
    monkeypatch.setitem(sys.modules, "omicverse", omicverse)
    monkeypatch.setitem(sys.modules, "cellphonedb", ModuleType("cellphonedb"))
    database_path = tmp_path / "cellphonedb.zip"
    database_path.write_bytes(b"database")
    adata = _adata()

    actual_results, actual_adata = run_cellphonedb_v5(
        adata,
        database_path=database_path,
        celltype_key="cluster",
        min_cell_fraction=0.005,
        min_genes=200,
        min_cells=3,
        iterations=1000,
        threshold=0.1,
        pvalue=0.05,
        threads=10,
        output_dir=tmp_path / "output",
        cleanup_temp=True,
    )

    assert actual_results is results
    assert actual_adata is returned_adata
    assert calls == [
        (
            adata,
            {
                "cpdb_file_path": str(database_path),
                "celltype_key": "cluster",
                "min_cell_fraction": 0.005,
                "min_genes": 200,
                "min_cells": 3,
                "iterations": 1000,
                "threshold": 0.1,
                "pvalue": 0.05,
                "threads": 10,
                "output_dir": str(tmp_path / "output"),
                "cleanup_temp": True,
            },
        )
    ]


def test_run_cellphonedb_v5_distinguishes_missing_resource_and_dependency(
    monkeypatch, tmp_path
):
    adata = _adata()
    missing_path = tmp_path / "missing.zip"

    with pytest.raises(FileNotFoundError, match="CellPhoneDB data asset"):
        run_cellphonedb_v5(
            adata,
            database_path=missing_path,
            celltype_key="cluster",
            min_cell_fraction=0.005,
            min_genes=200,
            min_cells=3,
            iterations=1000,
            threshold=0.1,
            pvalue=0.05,
            threads=10,
            output_dir=tmp_path,
            cleanup_temp=True,
        )

    database_path = tmp_path / "cellphonedb.zip"
    database_path.write_bytes(b"database")
    real_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name in {"omicverse", "cellphonedb"}:
            raise ModuleNotFoundError(f"No module named {name!r}", name=name)
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "omicverse", raising=False)
    monkeypatch.delitem(sys.modules, "cellphonedb", raising=False)
    monkeypatch.setattr(builtins, "__import__", blocked_import)

    with pytest.raises(ImportError, match=r"revise-svc\[cci\]"):
        run_cellphonedb_v5(
            adata,
            database_path=database_path,
            celltype_key="cluster",
            min_cell_fraction=0.005,
            min_genes=200,
            min_cells=3,
            iterations=1000,
            threshold=0.1,
            pvalue=0.05,
            threads=10,
            output_dir=tmp_path,
            cleanup_temp=True,
        )
