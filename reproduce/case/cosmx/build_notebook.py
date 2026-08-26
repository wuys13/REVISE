from __future__ import annotations

import argparse
import os
import tempfile
import textwrap
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
NOTEBOOK_NAME = "CosMx_SMI_267T_not_sp_SVC.ipynb"


def md(source: str):
    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip() + "\n")


def code(source: str):
    return nbf.v4.new_code_cell(textwrap.dedent(source).strip() + "\n")


SETUP = r'''
import os
from pathlib import Path

for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
    "BLIS_NUM_THREADS", "TBB_NUM_THREADS",
):
    os.environ[key] = "1"
os.environ["PYTHONHASHSEED"] = "42"
os.environ["MPLBACKEND"] = "Agg"

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc

SEED = 42

CASE_ROOT = Path(os.environ["REVISE_COSMX_CASE_ROOT"]).expanduser().resolve()
ST_PATH = CASE_ROOT / "prepared" / "CosMx_SMI_267T_not_ST.h5ad"
REFERENCE_PATH = CASE_ROOT / "prepared" / "sc_downsampled_deconvST_REVISE_panel6042.h5ad"
SVC_PATH = CASE_ROOT / "results" / "CosMx_SMI_267T_not" / "CosMx_SMI_267T_not_sp_SVC.h5ad"
'''


LOAD_AND_ALIGN = r'''
st = ad.read_h5ad(ST_PATH)
reference = ad.read_h5ad(REFERENCE_PATH)
svc = ad.read_h5ad(SVC_PATH)

shared_genes = svc.var_names[
    svc.var_names.isin(st.var_names)
    & svc.var_names.isin(reference.var_names)
]
raw_st = st[svc.obs_names, shared_genes].copy()
reference = reference[:, shared_genes].copy()
svc = svc[:, shared_genes].copy()

reference_labels = reference.obs["Level1"].astype(str).to_numpy()
formal_labels = svc.obs["Level1"].astype(str).to_numpy()
raw_st_labels = formal_labels.copy()
'''


PLOT = r'''
umap_sources = (
    ("CRC single-cell reference", reference.copy(), reference_labels),
    ("Raw CosMx ST | GA Level1", raw_st.copy(), raw_st_labels),
    ("REVISE sp-SVC | strength 0, graph α 0.8", svc.copy(), formal_labels),
)
prepared_sources = []
for title, prepared, labels in umap_sources:
    prepared.obs["Level1"] = labels
    sc.pp.normalize_total(prepared, target_sum=1e4)
    sc.pp.log1p(prepared)
    sc.pp.pca(prepared, n_comps=30, random_state=SEED)
    sc.pp.neighbors(
        prepared,
        n_neighbors=15,
        n_pcs=30,
        random_state=SEED,
    )
    sc.tl.umap(prepared, random_state=SEED)
    prepared_sources.append((title, prepared))

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (title, prepared) in zip(axes, prepared_sources):
    sc.pl.umap(
        prepared,
        color="Level1",
        ax=ax,
        show=False,
        title=title,
    )
plt.tight_layout()
plt.show()
'''


def build_notebook():
    notebook = nbf.v4.new_notebook()
    notebook["cells"] = [
        md(
            r'''
            # Case study: CosMx SMI 267T_not REVISE sp-SVC

            ## Notebook Guide

            **Purpose.** Load the prepared CosMx segmented-cell data, its matched
            single-cell reference, and the released formal sp-SVC object, then show
            three independently fitted expression-space UMAP panels.

            **Question.** How do the reference, prepared CosMx observations, and
            formal sp-SVC expression occupy their own expression-space views?

            **Method.** Reindex prepared ST rows to formal `svc.obs_names`, retain
            shared genes in formal order, and fit one independent UMAP per source.

            **Inputs.** Set `REVISE_COSMX_CASE_ROOT` to the case data directory. It
            must contain the prepared ST object, the prepared reference, and the
            formal sp-SVC object at the paths used below. CosMx rows are segmented
            cells; this notebook keeps the existing sp-SVC case interpretation.

            **Execution boundary.** The image in the plotting cell is carried over
            from the Phase 1 snapshot. The streamlined code is intentionally
            unexecuted here and has not been rerun.
            '''
        ),
        code(SETUP),
        md(
            r'''
            ## Load and align the three AnnData objects

            The formal output defines the observation order. The raw prepared ST
            object is reindexed to those `obs_names`, and all three objects are
            restricted to the shared genes in formal-output order before plotting.
            The raw-ST panel uses the formal Global Anchoring `Level1` labels so the
            three panels retain the original case-study comparison.
            '''
        ),
        code(LOAD_AND_ALIGN),
        md(
            r'''
            ## Independently fitted expression-space UMAPs

            Each source is normalized, log-transformed, reduced, and embedded on its
            own. Independent panel coordinates are descriptive and are not shared
            axes.
            '''
        ),
        code(PLOT),
        md(
            r'''
            ## Interpretation boundary

            **Direct observation.** The three panels show independently fitted
            embeddings with their displayed `Level1` labels.

            **Interpretation boundary.**

            The reference, raw-ST, and formal sp-SVC coordinates are independent:
            orientation, scale, and distances are not comparable across panels.
            Colors in the raw-ST panel are formal Global Anchoring assignments, not
            an independent raw annotation. The display is descriptive and does not
            establish newly resolved cell subtypes; it presents spatial-neighborhood
            expression reconstruction over segmented CosMx cells.
            '''
        ),
    ]
    notebook["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.10.14"},
    }
    return notebook


def write_notebook(notebook, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            nbf.write(notebook, handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, output)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def execute_notebook(notebook, *, timeout: int, kernel_name: str):
    from nbclient import NotebookClient

    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name=kernel_name,
        allow_errors=False,
        resources={"metadata": {"path": str(HERE.parents[2])}},
    )
    return client.execute()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the compact CosMx SMI sp-SVC showcase notebook."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE.parent / NOTEBOOK_NAME,
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute with nbclient before publishing the notebook.",
    )
    parser.add_argument("--kernel-name", default="python3")
    parser.add_argument("--timeout", type=int, default=1_800)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    notebook = build_notebook()
    if args.execute:
        notebook = execute_notebook(
            notebook,
            timeout=args.timeout,
            kernel_name=args.kernel_name,
        )
    write_notebook(notebook, args.output.expanduser().resolve())
    state = "executed" if args.execute else "unexecuted"
    print(f"Wrote {state} notebook: {args.output}")


if __name__ == "__main__":
    main()
