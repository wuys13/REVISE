# REVISE

<p align="center">
  <img src="./logo/REVISE.png" alt="REVISE logo" height="240" />
  <img src="./logo/Sim2Real-ST.png" alt="Sim2Real-ST logo" height="240" />
  <img src="./logo/SVC.png" alt="SVC logo" height="240" />
</p>

[![PyPI](https://img.shields.io/pypi/v/revise-svc.svg)](https://pypi.org/project/revise-svc/)
[![Documentation Status](https://readthedocs.org/projects/revise-svc/badge/?version=latest)](https://revise-svc.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

REVISE (REconstruction via Vision-integrated Spatial Estimation) reconstructs
**Spatially-inferred Virtual Cells (SVCs)** from spatial transcriptomics data
and a matched single-cell RNA-seq reference. Spatial or morphology-derived
priors can be used when available.

## Start Here: Reconstruct an SVC

Provide one spatial-transcriptomics H5AD and one matched single-cell reference,
then run [`reconstruct.py`](reconstruct.py). Choose the 1.x reconstruction type
from your input rows: **hST-like bins/pseudo-cells → `sp-SVC`; iST segmented
cells → `sc-SVC`; sST spots → `sc-SVC-sr`.** The labels hST, iST, and sST are
data-selection guidance, not current CLI values.

First validate the resolved inputs without running reconstruction:

This example reads `data/sample_st.h5ad` and `data/sc_ref.h5ad`.

```bash
python reconstruct.py \
  --svc-type sp-SVC \
  --sample-name sample \
  --data-root data \
  --st-file st.h5ad \
  --sc-ref-file sc_ref.h5ad \
  --output-root output \
  --dry-run
```

Replace `sp-SVC` using the guidance below. Remove `--dry-run` to reconstruct.
The equivalent installed command is `revise-reconstruct`.

Every successful full reconstruction publishes:

```text
output/sample/SVC.h5ad
```

`--dry-run` does not publish `SVC.h5ad`. It checks the metadata-level input
contract and required arrays; it does not fully scan expression values.

The associated `provenance.json` records the selected type, resolved route,
configuration, inputs, stages, and artifacts. Passing preflight does not prove
that the complete reconstruction or a biological interpretation will succeed.

<details>
<summary><strong>Which --svc-type should I choose?</strong></summary>

| Your ST data | Typical platform and input rows | `--svc-type` | What REVISE reconstructs |
| --- | --- | --- | --- |
| **hST-like** high-resolution sequencing data | Visium HD; each row is a bin or pseudo-cell | `sp-SVC` | Reference-annotated expression with spatial refinement for the retained high-resolution units |
| **iST** imaging-based data | iST-like segmented-cell inputs, such as Xenium, CosMx, or MERFISH | `sc-SVC` | Segmented-cell positions combined with reference-informed cell-state refinement and gene completion |
| **sST** spot-based data | Visium; each row is a multi-cell spot | `sc-SVC-sr` | Reference-informed virtual-cell expression and cell-type composition within each spot |

`sc-SVC-sr` does not by itself infer true sub-spot cell positions. When the ST
H5AD contains segmentation-derived cell centers, REVISE uses them. Virtual
cells without a supplied center remain at their source spot center.

</details>

<details>
<summary><strong>Input file layout and AnnData requirements</strong></summary>

The example command resolves these paths:

```text
data/
├── sample_st.h5ad
└── sc_ref.h5ad
```

`--sample-name sample` and `--st-file st.h5ad` resolve to
`data/sample_st.h5ad`; the reference resolves directly to
`data/sc_ref.h5ad`.

- Both inputs must have non-empty `X`, unique `obs_names`, unique `var_names`,
  and at least one shared gene.
- The ST input must contain finite two-dimensional coordinates in
  `obsm["spatial"]`.
- The reference must contain the broad annotation selected by
  `--cell-type-col` (default `obs["Level1"]`) on every route. `sc-SVC` and
  `sc-SVC-sr` also require the refined annotation selected by
  `--sub-cell-type-col` (default `obs["Level2"]`).
- If the reference contains the default `Patient` column, at least one row must
  match `--sample-name`; use `--patient-key` to select another column.
- For sST inputs with segmentation-derived centers, the optional
  `uns["revise_cell_locations"]` table uses unique `cell_id` values as its
  index and contains `spot_name`, `x`, and `y`. Its cell IDs must agree with
  `uns["all_cells_in_spot"]`; `x/y` must use the same coordinate system and
  scale as `obsm["spatial"]`. Missing centers fall back to the spot center.
  Without a `PM_on_cell.csv`, these coordinates are retained while cell types
  are assigned to the existing rows by a seeded random permutation of each
  spot's inferred quota.

</details>

<details>
<summary><strong>Frequently used reconstruction parameters</strong></summary>

- `--seed`: controls deterministic random choices; default `42`.
- `--ot-method pot|tacco`: selects one OT implementation for both Global
  Anchoring and Local Refinement. TACCO requires the optional installation
  group described below.
- `--select-ct`: for `sc-SVC`, reconstruct one broad cell type or use the
  default `all`.
- `--cell-type-col`: selects the broad reference annotation column for all
  three routes.
- `--sub-cell-type-col`: selects the refined annotation required by `sc-SVC`
  and `sc-SVC-sr`. `sc-SVC-sr` validates this column, but its current
  composition assignment is driven by the broad column.
- `--sc-mapping mean|random`: for `sc-SVC`, map each reconstructed spatial row
  to its cluster mean expression or to a seeded same-cluster reference row.
- `--set KEY=VALUE`: advanced configuration override. High-level CLI options
  cannot be contradicted through `--set`.

Run `python reconstruct.py --help` for the complete command contract.

</details>

REVISE has two public workflows:

| Workflow | Goal | Entry points | Main results |
| --- | --- | --- | --- |
| **Benchmark** | Reproduce Sim2Real-ST evaluations across six confounding factors | [`reproduce/benchmark_main.py`](reproduce/benchmark_main.py), [`reproduce/benchmark_main.sh`](reproduce/benchmark_main.sh), [`reproduce/benchmark/`](reproduce/benchmark/) | Per-gene PCC, SSIM, MSE, and NRMSE |
| **Application** | Reconstruct SVCs and perform real-data analyses | [`reconstruct.py`](reconstruct.py), `revise-reconstruct`, [`reproduce/case/`](reproduce/case/) | `SVC.h5ad`, run provenance, and notebook figures |

Documentation: <https://revise-svc.readthedocs.io/en/latest/>

Dataset and reproduced results: <https://zenodo.org/records/17705737>

## What REVISE Covers

Sim2Real-ST benchmarks six confounding factors across spatial transcriptomics
platform types:

- Spatially heterogeneous factors: image segmentation artifacts and bin-to-cell
  assignment errors.
- Spatially homogeneous factors: spot size, batch effect, gene panel
  limitation, and gene dropout.

![Spatial transcriptomics limitations](png/ST_limitations.png)

<p align="center">Confounding factors that limit current spatial transcriptomics technologies</p>

REVISE reconstructs three complementary SVC types:

- `sp-SVC`: spatial refinement for high-resolution spatial transcriptomics.
- `sc-SVC`: molecular completion and cell-state refinement for imaging-based
  spatial transcriptomics.
- `sc-SVC-sr`: spot-level super-resolution reconstruction.

![REVISE overview](png/REVISE_overview.png)

<p align="center">Overview of the REVISE framework</p>

Every application reconstruction publishes the same stable filename:

```text
<output-root>/<sample-name>/SVC.h5ad
```

The run's `provenance.json` records whether that result is an `sp-SVC`,
`sc-SVC`, or `sc-SVC-sr`, together with the resolved route, configuration,
inputs, stages, and artifacts.

## Installation

The unified CLI and optional groups described below are the current source
contract. Install that exact code from the repository:

```bash
git clone https://github.com/wuys13/REVISE.git
cd REVISE
python -m pip install .
```

The published package can be installed with `python -m pip install revise-svc`,
but releases can lag the repository; check the installed version before
expecting the current CLI or optional groups.

The base source installation contains reconstruction, the default OT
implementation, clustering, and the core scientific stack. Install additional
capabilities only when needed:

| Capability | Installation | Purpose |
| --- | --- | --- |
| Additional OT implementation | `python -m pip install ".[tacco]"` | Adds another selectable OT implementation, such as TACCO |
| Pathway analysis | `python -m pip install ".[pathway]"` | Dependencies used by pathway analysis notebooks |
| Cell-cell interaction analysis | `python -m pip install ".[cci]"` | Dependencies used by CCI notebooks; databases and reference resources are prepared separately |
| Trajectory analysis | `python -m pip install ".[trajectory]"` | Dependencies used by trajectory analysis notebooks |
| SpatialData input | `python -m pip install ".[spatialdata]"` | SpatialData/Zarr input support |

After a matching package version is published, replace `.` with `revise-svc`
in those commands. For development:

```bash
python -m pip install -e ".[dev]"
```

Installation does not download research data or external analysis databases.

## Quick Start

### Benchmark

From a source checkout, run one Sim2Real-ST confounding family:

```bash
python reproduce/benchmark_main.py \
  --confounding segmentation \
  --data-root raw_data/Sim2Real-ST \
  --sample-name P2CRC/cut_part1 \
  --dataset-task segmentation \
  --output-root output/benchmark
```

Supported values are `segmentation`, `bin2cell`, `batch_effect`, `spot_size`,
`gene_panel`, and `gene_dropout`. Run the bounded multi-family launcher with:

```bash
bash reproduce/benchmark_main.sh
```

The analysis notebooks are under [`reproduce/benchmark/`](reproduce/benchmark/).

### Application examples and notebooks

![SVC applications](png/SVC_applications.png)

<p align="center">Biological insights enabled by SVC reconstruction</p>

Application reconstruction and downstream analysis notebooks are under
[`reproduce/case/`](reproduce/case/).

## Python API

```python
from revise.framework import REVISEPipeline

pipeline = REVISEPipeline()
svc = pipeline.run(
    profile="application_sc",
    runtime_overrides={"platform": "sc_svc", "confounding": "segmentation"},
    io_overrides={
        "data_root": "raw_data/Real_application",
        "output_root": "output/sc_SVC_case",
        "sample_name": "P2CRC",
        "st_file": "Xenium.h5ad",
        "sc_ref_file": "adata_sc_all_reanno.h5ad",
        "patient_key": "Patient",
    },
    set_overrides=["sc.select_ct=T"],
)
```

## Reproduction Notebooks

The curated paper notebooks are tracked under [`reproduce/`](reproduce/).
They require the corresponding data and, for some downstream analyses, the
optional installation groups described above.

Their presence preserves the paper workflows; it does not mean the current
source checkout has rerun every real-data analysis. Real-data end-to-end
validation remains a separate release step.

## Repository Layout

- `revise/`: installable reconstruction and analysis package.
- `reproduce/`: benchmark launchers and benchmark/application notebooks.
- `docs/`: detailed user and method documentation.
- `logo/`, `png/`: public project and scientific figures.
- `tests/`: behavioral, scientific-contract, packaging, and CLI tests.
- `.github/`: continuous integration.
- `constraints/`: tested Python 3.10/3.11 dependency constraints.

## License

REVISE is released under the [MIT License](LICENSE). Third-party datasets and
notebook resources may have separate terms.
