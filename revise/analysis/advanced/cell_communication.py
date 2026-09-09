from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd


_RESULT_TABLES = ("pvalues", "significant_means", "interaction_scores")


def _require_cellphonedb_dependencies():
    try:
        import cellphonedb  # noqa: F401
        import omicverse
    except ImportError as exc:
        raise ImportError(
            "CellPhoneDB analysis requires OmicVerse and CellPhoneDB; install "
            'them with `python -m pip install "revise-svc[cci]"`.'
        ) from exc
    return omicverse


def normalize_cellphonedb_tables(
    results: Mapping[str, pd.DataFrame], cell_groups: Sequence[str]
) -> dict[str, pd.DataFrame]:
    """Return the three consumed CellPhoneDB tables on a shared cell-pair axis."""
    groups = [str(group) for group in cell_groups]
    if not groups or len(set(groups)) != len(groups):
        raise ValueError("cell_groups must contain unique group names")

    allowed_groups = set(groups)
    normalized: dict[str, pd.DataFrame] = {}
    shared_pair_columns: list[str] | None = None

    for table_name in _RESULT_TABLES:
        if table_name not in results:
            raise KeyError(f"CellPhoneDB results are missing '{table_name}'")
        table = results[table_name]
        if not isinstance(table, pd.DataFrame):
            raise ValueError(f"CellPhoneDB result '{table_name}' must be a DataFrame")
        missing_gene_columns = {"gene_a", "gene_b"}.difference(table.columns)
        if missing_gene_columns:
            missing = ", ".join(sorted(missing_gene_columns))
            raise ValueError(
                f"CellPhoneDB result '{table_name}' is missing columns: {missing}"
            )

        pair_columns = [
            column
            for column in table.columns
            if isinstance(column, str) and column.count("|") == 1
        ]
        if not pair_columns or any(
            source not in allowed_groups or target not in allowed_groups
            for source, target in (column.split("|") for column in pair_columns)
        ):
            raise ValueError(
                f"CellPhoneDB result '{table_name}' has an incompatible pair-column axis"
            )
        if shared_pair_columns is None:
            shared_pair_columns = pair_columns
        elif pair_columns != shared_pair_columns:
            raise ValueError(
                "CellPhoneDB result tables must use the same pair-column sequence"
            )

        complete = table.dropna(subset=["gene_a", "gene_b"])
        normalized_table = complete.loc[:, pair_columns].copy()
        normalized_table.index = (
            complete["gene_a"].astype(str) + "_" + complete["gene_b"].astype(str)
        )
        normalized_table.index.name = "LR_pair"
        normalized[table_name] = normalized_table

    return normalized


def run_cellphonedb_v5(
    adata,
    *,
    database_path: str | Path,
    celltype_key: str,
    min_cell_fraction: float,
    min_genes: int,
    min_cells: int,
    iterations: int,
    threshold: float,
    pvalue: float,
    threads: int,
    output_dir: str | Path,
    cleanup_temp: bool,
) -> tuple[Mapping[str, pd.DataFrame], Any]:
    """Run the notebook's OmicVerse/CellPhoneDB-v5 path and validate its tables."""
    database = Path(database_path)
    if not database.is_file():
        raise FileNotFoundError(f"CellPhoneDB data asset is missing: {database}")
    if celltype_key not in adata.obs:
        raise ValueError(f"celltype_key '{celltype_key}' not found in adata.obs")

    omicverse = _require_cellphonedb_dependencies()
    output = omicverse.single.run_cellphonedb_v5(
        adata,
        cpdb_file_path=str(database),
        celltype_key=celltype_key,
        min_cell_fraction=min_cell_fraction,
        min_genes=min_genes,
        min_cells=min_cells,
        iterations=iterations,
        threshold=threshold,
        pvalue=pvalue,
        threads=threads,
        output_dir=str(output_dir),
        cleanup_temp=cleanup_temp,
    )

    if not isinstance(output, tuple) or len(output) != 2:
        raise RuntimeError(
            "OmicVerse CellPhoneDB analysis did not return (results, adata)"
        )
    results, result_adata = output
    if not isinstance(results, Mapping):
        raise ValueError("CellPhoneDB results must be a mapping of result tables")

    return results, result_adata


__all__ = ["normalize_cellphonedb_tables", "run_cellphonedb_v5"]
