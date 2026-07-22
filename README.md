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

REVISE has two public workflows:

| Workflow | Goal | Entry points | Main results |
| --- | --- | --- | --- |
| **Benchmark** | Reproduce Sim2Real-ST evaluations across six confounding factors | [`benchmark_main.py`](benchmark_main.py), [`benchmark_main.sh`](benchmark_main.sh), [`reproduce/benchmark/`](reproduce/benchmark/) | Per-gene PCC, SSIM, MSE, and NRMSE |
| **Application** | Reconstruct SVCs and perform real-data analyses | `revise-reconstruct`, [`reconstruct.py`](reconstruct.py), compatibility wrappers, [`reproduce/case/`](reproduce/case/) | `SVC.h5ad`, run provenance, and notebook figures |

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

REVISE reconstructs two complementary SVC types:

- `sp-SVC`: spatial refinement for high-resolution spatial transcriptomics.
- `sc-SVC`: molecular completion and cell-state refinement for imaging-based
  and spot-based spatial transcriptomics.

![REVISE overview](png/REVISE_overview.png)

<p align="center">Overview of the REVISE framework</p>

Every application reconstruction publishes the same stable filename:

```text
<output-root>/<sample-name>/SVC.h5ad
```

The run's `provenance.json` records whether that result comes from the `hST`,
`iST`, or `sST` route, together with the resolved configuration, inputs,
stages, and artifacts.

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
python benchmark_main.py \
  --confounding segmentation \
  --data-root raw_data/Sim2Real-ST \
  --sample-name P2CRC/cut_part1 \
  --dataset-task segmentation \
  --output-root output/benchmark
```

Supported values are `segmentation`, `bin2cell`, `batch_effect`, `spot_size`,
`gene_panel`, and `gene_dropout`. Run the bounded multi-family launcher with:

```bash
bash benchmark_main.sh
```

The analysis notebooks are under [`reproduce/benchmark/`](reproduce/benchmark/).

### Application

![SVC applications](png/SVC_applications.png)

<p align="center">Biological insights enabled by SVC reconstruction</p>

Prepare an ST AnnData file and a matched single-cell reference, then run the
installed command:

```bash
revise-reconstruct \
  --platform hST \
  --sample-name sample \
  --data-root data \
  --st-file st.h5ad \
  --sc-ref-file sc_ref.h5ad \
  --output-root output
```

`--platform` selects the reconstruction route. The public result is always
`output/sample/SVC.h5ad`; its `hST`, `iST`, or `sST` classification is recorded
in `provenance.json`. From a source checkout, `python reconstruct.py` delegates
to the same command implementation.

The paper-compatible wrappers are retained only for reproducing historical
notebooks. They keep their legacy output names and are not part of the
canonical `SVC.h5ad` output contract:

- [`application_sp_SVC_recon.py`](application_sp_SVC_recon.py)
- [`application_sc_SVC_recon.py`](application_sc_SVC_recon.py)

Application reconstruction and downstream analysis notebooks are under
[`reproduce/case/`](reproduce/case/).

Optional morphology priors can be prepared through the installed package
command:

```bash
revise-build-histology-priors \
  --st-h5ad st.h5ad \
  --mask segmented_cells.tif \
  --out-h5ad st_with_histology_priors.h5ad
```

Biology-facing post-reconstruction metrics use the package-owned analysis
implementation:

```bash
revise-compute-biological-metrics \
  --input-h5ad output/sample/SVC.h5ad \
  --output-dir output/sample/biological_metrics
```

## Python API

```python
from revise.framework import REVISEPipeline

pipeline = REVISEPipeline()
svc = pipeline.run(
    profile="application_sc",
    runtime_overrides={"platform": "iST", "confounding": "segmentation"},
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
- `reproduce/`: benchmark and application notebooks.
- `docs/`: detailed user and method documentation.
- `logo/`, `png/`: public project and scientific figures.
- `tests/`: behavioral, scientific-contract, packaging, and CLI tests.
- `.github/`: continuous integration.
- `constraints/`: tested Python 3.10/3.11 dependency constraints.

## License

REVISE is released under the [MIT License](LICENSE). Third-party datasets and
notebook resources may have separate terms.
