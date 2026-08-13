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

## Quick start

REVISE supports Python 3.10 and 3.11. A clean environment is recommended.

<details>
<summary><strong>Create an environment with uv or conda</strong></summary>

With uv:

```bash
uv venv --python 3.11
source .venv/bin/activate
```

Or with conda:

```bash
conda create -n revise python=3.11 -y
conda activate revise
```

</details>

Install REVISE:

```bash
pip install revise-svc
```

Choose the application template that matches the ST data:

| ST data | Template | REVISE route |
| --- | --- | --- |
| Xenium T-cell data | [`configs/application/Xenium_T.yaml`](configs/application/Xenium_T.yaml) | `sc-SVC` |
| Xenium fibroblast data | [`configs/application/Xenium_Fib.yaml`](configs/application/Xenium_Fib.yaml) | `sc-SVC` |
| Xenium mono/macro data | [`configs/application/Xenium_Mono.yaml`](configs/application/Xenium_Mono.yaml) | `sc-SVC` |
| Visium HD bins or pseudo-cells | [`configs/application/VisiumHD.yaml`](configs/application/VisiumHD.yaml) | `sp-SVC` |
| Visium or another multi-cell spot assay | [`configs/application/Visium.yaml`](configs/application/Visium.yaml) | `sc-SVC-sr` |

These files are present in a source checkout. With a package-only installation,
copy the selected file from the packaged `revise.application/templates`
resource into your project with Python's standard `importlib.resources`, edit
the copy, and pass its local path to `--config` (see the
[installation guide](docs/source/installation.rst)). A bare official template
name first uses `configs/application/<name>.yaml` in the launch directory and
then falls back to the packaged resource.

Edit the selected YAML before running it. The fields you normally change are
`paths.root_dir`, the exact ST/reference paths, optional reference
`filter_column`/`filter_value`, preprocessing thresholds, annotation columns,
and `output.dir`/`output.name`. The three Xenium files describe the same Xenium
assay with different reconstruction targets: `T`, `Fibroblast`, or
`Mono_Macro`. Confirm that the selected label exists in the full reference
`obs` column. Visium HD and Visium expose only `local_refinement.strength`.
Visium's `0.0` still runs LR; it only disables assignment-posterior
strengthening of the LR cost. The official Xenium templates filter the shared
reference with `Patient == P2CRC`; the key names are generic and may describe a
different reference cohort field in a custom YAML.

From a source checkout, run the YAML directly:

```bash
python reconstruct.py --config configs/application/VisiumHD.yaml
```

After installation, the installed command `revise-reconstruct` is the
equivalent entry name:

```bash
revise-reconstruct --config configs/application/VisiumHD.yaml
```

The three routes publish:

```text
output/P1CRC_HD_sp-SVC.h5ad
output/visium_mouse_brain_revise/REVISEVisiumMouseBrain_sc-SVC-sr.h5ad
```

Standard `sc-SVC` publishes its two reconstruction carriers separately:

```text
output/P2CRC_Xenium/T/spatial.h5ad
output/P2CRC_Xenium/T/expr.h5ad
```

The associated `provenance.json` keeps the application YAML identity under
`application_config` (`source_path`, `source_sha256`, root and resolved paths,
effective request/hash, and action). Top-level `engine_defaults_hash`,
`authority_hash`, `algorithm_config_hash`, and `effective_config_hash` identify
the package engine authority and resolved run. [`reconstruct.py`](reconstruct.py)
is the canonical Application implementation, and the installed
`revise-reconstruct` command points to its same `main()` implementation.

<details>
<summary><strong>What does each template reconstruct?</strong></summary>

| Your ST data | Typical platform and input rows | Route | What REVISE reconstructs |
| --- | --- | --- | --- |
| **hST-like** high-resolution sequencing data | Visium HD; each row is a bin or pseudo-cell | `sp-SVC` | Reference-annotated expression with spatial refinement for the retained high-resolution units |
| **iST** imaging-based data | iST-like segmented-cell inputs, such as Xenium, CosMx, or MERFISH | `sc-SVC` | Segmented-cell positions combined with reference-informed cell-state refinement and gene completion |
| **sST** spot-based data | Visium; each row is a multi-cell spot | `sc-SVC-sr` | Reference-informed virtual-cell expression and cell-type composition within each spot |

`sc-SVC-sr` does not by itself infer true sub-spot cell positions. When the ST
H5AD contains segmentation-derived cell centers, REVISE uses them. Virtual
cells without a supplied center remain at their source spot center.

</details>

<details>
<summary><strong>Application YAML inputs and AnnData requirements</strong></summary>

Set `inputs.st.path` and `inputs.reference.path` explicitly. `paths.root_dir: .` means the launch current
working directory, not the application YAML directory; otherwise it must be an
existing absolute directory. Input and output values are relative children of
that root and cannot escape it with `..`. Package code in
`revise.config.authority` owns engine defaults and routing.

- Both inputs must have non-empty `X`, unique `obs_names`, unique `var_names`,
  and at least one shared gene.
- The ST input must contain finite two-dimensional coordinates in
  `obsm["spatial"]`.
- Every route requires `global_anchoring.broad_column` in reference `obs`. Only
  standard `sc-SVC` accepts and requires `local_refinement.subtype_column`.
  `sc-SVC-sr` composition and expression allocation use the broad assignment.
- A reference `filter_column` and `filter_value` must be provided together or
  both omitted. When present, filtering runs before the spatial and reference
  count/gene preprocessing steps.
- For sST inputs with segmentation-derived centers, the optional
  `uns["revise_cell_locations"]` table uses unique `cell_id` values as its
  index and contains `spot_name`, `x`, and `y`. Its cell IDs must agree with
  `uns["all_cells_in_spot"]`; `x/y` must use the same coordinate system and
  scale as `obsm["spatial"]`. Missing centers fall back to the spot center. The
  optional probability prior can be provided with the exact
  `inputs.pm_on_cell.path`. Its axes must exactly equal the active virtual-cell
  IDs and normalized cell types; values must be finite probabilities whose rows
  sum to one. REVISE reorders exact axes but does not clip or normalize PM.
  Without that field, coordinates are retained while cell types are assigned
  to the existing rows by a seeded random permutation of each spot's inferred
  quota.

</details>

<details>
<summary><strong>SpatialData and algorithm controls</strong></summary>

ST `format` accepts `h5ad`, `spatialdata`, or `auto`; the reference remains
H5AD. SpatialData table and spatial-element selections live under
`inputs.st.spatialdata.table` and `inputs.st.spatialdata.element`. REVISE does
not expose a coordinate-conversion promise through this entry point.

`execution.seed` and `algorithm.ot_method` are optional. The latter controls
both GA and LR; omitting it keeps the selected engine profile authoritative.
Standard `sc-SVC` defaults to TACCO and requires the `tacco` extra; REVISE never
changes solvers automatically. Unknown or route-inapplicable application
fields are errors; the engine YAML and generic key/value overrides are not
public application controls.

</details>

<details>
<summary><strong>Source installation, optional capabilities, and development setup</strong></summary>

To install the current repository source:

```bash
git clone https://github.com/wuys13/REVISE.git
cd REVISE
pip install .
```

The base package contains reconstruction, the POT implementation, clustering,
and the core scientific stack. Install additional capabilities only when
needed:

| Capability | Installation | Purpose |
| --- | --- | --- |
| Standard sc-SVC default solver | `pip install "revise-svc[tacco]"` | Installs TACCO 0.5.0, required by the default standard sc-SVC route |
| Pathway analysis | `pip install "revise-svc[pathway]"` | Dependencies used by pathway analysis notebooks |
| Cell-cell interaction analysis | `pip install "revise-svc[cci]"` | Dependencies used by CCI notebooks; databases and reference resources are prepared separately |
| Trajectory analysis | `pip install "revise-svc[trajectory]"` | Dependencies used by trajectory analysis notebooks |
| SpatialData input | `pip install "revise-svc[spatialdata]"` | SpatialData/Zarr input support |

For optional groups from a source checkout, replace `revise-svc` with `.`. For
development:

```bash
pip install -e ".[dev]"
```

Installation does not download research data or external analysis databases.

</details>

Detailed documentation: <https://revise-svc.readthedocs.io/en/latest/>

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

Application output names are declared in YAML. `sp-SVC` and `sc-SVC-sr` publish:

```text
<output-dir>/<output-name>.h5ad
```

If `output.name` is omitted, a single-output route uses `svc.h5ad`.
`sc-SVC` publishes `<output-name>_spatial.h5ad` and
`<output-name>_expr.h5ad`; without `output.name` it uses `spatial.h5ad` and
`expr.h5ad`. The run's
`provenance.json` records each result role together with the resolved route,
configuration, inputs, stages, and artifacts.

## Reproduce

Paper datasets and reproduced results are available from
<https://zenodo.org/records/17705737>. The repository keeps benchmark launchers
and curated notebooks under [`reproduce/`](reproduce/); see
[`reproduce/README.md`](reproduce/README.md) for the entry-point map.

<details>
<summary><strong>Run the Sim2Real-ST benchmarks</strong></summary>

From a source checkout, run one confounding family:

```bash
python reproduce/benchmark_main.py \
  --config configs/benchmark/segmentation.yaml \
  --data-root raw_data/Sim2Real-ST \
  --sample-name P2CRC/cut_part1 \
  --dataset-task segmentation \
  --evaluate true \
  --output-root output/benchmark
```

Supported values are `segmentation`, `bin2cell`, `batch_effect`, `spot_size`,
`gene_panel`, and `gene_dropout`. Run the bounded multi-family launcher with:

```bash
bash reproduce/benchmark_main.sh
```

The analysis notebooks are under [`reproduce/benchmark/`](reproduce/benchmark/).

</details>

<details>
<summary><strong>Open the application and downstream-analysis notebooks</strong></summary>

![SVC applications](png/SVC_applications.png)

<p align="center">Biological insights enabled by SVC reconstruction</p>

Application reconstruction and downstream-analysis notebooks are under
[`reproduce/case/`](reproduce/case/). Some require the optional `pathway`,
`cci`, or `trajectory` installation groups and their corresponding external
reference resources.

</details>

<details>
<summary><strong>Reproduction scope and validation status</strong></summary>

The notebooks preserve the paper workflows, but their presence does not mean
that the current source checkout has rerun every real-data analysis. Real-data
end-to-end validation remains a separate release step. Installation does not
download the paper datasets, reproduced results, or external analysis
databases.

</details>

## Python API

<details>
<summary><strong>Run the reconstruction pipeline from Python</strong></summary>

```python
from reconstruct import run_application

spatial_svc = run_application("configs/application/VisiumHD.yaml")
spatial, expression = run_application("configs/application/Xenium_T.yaml")

# sp-SVC and sc-SVC-sr return one AnnData; sc-SVC returns this fixed pair.
```

Application callers provide one YAML. `REVISEPipeline.run()` remains the
engine-level interface used by Benchmark and advanced integrations; it is not
the Application configuration API. Benchmark execution continues to use the
separate `cf` selector for its confounding family.

</details>

## Repository Layout

<details>
<summary><strong>Show the top-level repository structure</strong></summary>

- `revise/`: installable reconstruction and analysis package.
- `reproduce/`: benchmark launchers and benchmark/application notebooks.
- `docs/`: detailed user and method documentation.
- `logo/`, `png/`: public project and scientific figures.
- `tests/`: behavioral, scientific-contract, packaging, and CLI tests.
- `.github/`: continuous integration.
- `constraints/`: tested Python 3.10/3.11 dependency constraints.

</details>

## License

REVISE is released under the [MIT License](LICENSE). Third-party datasets and
notebook resources may have separate terms.
