# REVISE

<p align="center">
  <img src="./logo/REVISE.png" alt="REVISE logo" width="30%" />
  <img src="./logo/Sim2Real-ST.png" alt="Sim2Real-ST logo" width="30%" />
  <img src="./logo/SVC.png" alt="SVC logo" width="30%" />
</p>

[![PyPI](https://img.shields.io/pypi/v/revise-svc.svg)](https://pypi.org/project/revise-svc/)
[![Documentation Status](https://readthedocs.org/projects/revise-svc/badge/?version=latest)](https://revise-svc.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

[REVISE documentation](https://revise-svc.readthedocs.io/en/latest/) | [Dataset](https://zenodo.org/records/21921802) | Paper

REVISE reconstructs Spatially-inferred Virtual Cells (SVCs) from spatial
transcriptomics data and a matched single-cell reference. Choose the workflow
from what one row of the spatial data represents, edit a small YAML request,
and run one command. The same project also provides the Sim2Real-ST benchmark
workflows used to study reconstruction under controlled confounding factors.

## Install

REVISE supports Python 3.10 and 3.11.

```bash
python -m pip install revise-svc
python -m pip install "revise-svc[tacco]"  # alternative TACCO OT method
```

The base package includes POT, the default OT implementation. REVISE also
supports another OT method, TACCO; the maintained Xenium template selects it
explicitly, so install `revise-svc[tacco]` before running that template.

<details>
<summary>Optional data-reading or downstream bioinformatics dependencies</summary>

Install these only for the relevant data format or downstream analysis.

```bash
python -m pip install "revise-svc[spatialdata]"  # SpatialData/Zarr input
python -m pip install "revise-svc[pathway]"      # pathway analysis
python -m pip install "revise-svc[cci]"          # cell-cell interaction analysis
python -m pip install "revise-svc[trajectory]"   # trajectory analysis
```

</details>

See the [installation guide](https://revise-svc.readthedocs.io/en/latest/source/installation.html)
for source installation, optional dependencies, and documentation builds.

## Data downloads

Download the source material or reproduced H5AD that matches your workflow.
The records complement a YAML request; they do not configure a new run for you.

| Material | Download |
| --- | --- |
| Sim2Real-ST benchmark | [Zenodo](https://zenodo.org/records/21921802) |
| Reproduced benchmark results | [Zenodo](https://zenodo.org/records/21921802) |
| Real-world ST datasets | [Zenodo](https://zenodo.org/records/21921802) |
| Reproduced sp-SVC H5AD | [Zenodo](https://zenodo.org/records/18389835) |
| Reproduced sc-SVC H5AD | [Zenodo](https://zenodo.org/records/18389211) |

## Quick run

### Choose the application template that matches the ST data

Run from a local editable YAML copy. With a published package, first use the
[template-copy step](https://revise-svc.readthedocs.io/en/latest/source/installation.html#application-templates)
to copy one packaged template into your working directory, then use
`revise-reconstruct`. In a source checkout, the maintained copies are under
`configs/application/`; the source-checkout equivalent commands are shown with
each template.

| ST data | Start from | Public route |
| --- | --- | --- |
| Visium HD bins or pseudo-cells | [`configs/application/VisiumHD.yaml`](configs/application/VisiumHD.yaml) | `sp-SVC` |
| Xenium segmented cells | [`configs/application/Xenium.yaml`](configs/application/Xenium.yaml) | `sc-SVC`, cluster mode |
| Visium multi-cell spots | [`configs/application/Visium.yaml`](configs/application/Visium.yaml) | `sc-SVC`, sr mode |

`sc-SVC` always names its mode explicitly. Cluster mode is for segmented
imaging-ST cells; sr mode reconstructs virtual cells within multi-cell spots.
REVISE does not infer the choice from a filename or silently repair a
mismatched request.

Before running any YAML, update the input paths, reference filter, annotation
columns, preprocessing thresholds, and output root. The
[Application Reference](https://revise-svc.readthedocs.io/en/latest/source/application-reference.html)
defines the accepted fields and output contract.

For Visium HD bins or pseudo-cells, use `VisiumHD.yaml` and set the spatial and
reference paths plus the broad annotation column.

```bash
revise-reconstruct --config VisiumHD.yaml
# Source-checkout equivalent:
python reconstruct.py --config configs/application/VisiumHD.yaml
```

For Xenium segmented cells, use `Xenium.yaml`, set its input/filter/annotation
fields, and select one concrete broad cell type. The cluster route automatically
writes to `T/`, `Fibroblast/`, or `Mono_Macro/` below the configured output root.

```bash
revise-reconstruct --config Xenium.yaml --select-ct T
revise-reconstruct --config Xenium.yaml --select-ct Fibroblast
revise-reconstruct --config Xenium.yaml --select-ct "Mono/Macro"
# Source-checkout equivalent:
python reconstruct.py --config configs/application/Xenium.yaml --select-ct T
python reconstruct.py --config configs/application/Xenium.yaml --select-ct Fibroblast
python reconstruct.py --config configs/application/Xenium.yaml --select-ct "Mono/Macro"
```

For multi-cell Visium spots, use `Visium.yaml` and review the optional PM prior
when your run uses one.

```bash
revise-reconstruct --config Visium.yaml
# Source-checkout equivalent:
python reconstruct.py --config configs/application/Visium.yaml
```

`--config` is required. `--select-ct` is accepted only for `sc-SVC` cluster
mode and overrides the Xenium YAML's `local_refinement.select_cell_type` with
one concrete broad cell type.

## What is written

| Route | Files returned and published |
| --- | --- |
| `sp-SVC` | one `<output.dir>/<output.name>.h5ad`, or `svc.h5ad` when `output.name` is omitted |
| `sc-SVC`, cluster mode | `<final-dir>/<name>_spatial.h5ad` and `<final-dir>/<name>_expr.h5ad`, or `spatial.h5ad` and `expr.h5ad` without a name |
| `sc-SVC`, sr mode | one `<output.dir>/<output.name>.h5ad`, or `svc.h5ad` when `output.name` is omitted |

Each H5AD carries its normalized route and mode metadata and links to the
run's `provenance.json`. A run is successful only when its promised artifact
or artifacts and succeeded manifest exist.

For cluster mode, `output.dir` is the base directory and the final selected-cell-type
subdirectory is appended automatically after label normalization.

## Framework and supported data shapes

REVISE uses one reconstruction lifecycle for three Application choices and
for Sim2Real-ST Benchmark routes. Application handles user data and publishes
the promised artifact shape; Benchmark owns its own experimental-case
preparation and metrics.

![REVISE framework overview](png/REVISE_overview.png)

The concise scientific distinction is: Visium HD uses spatial refinement of
high-resolution units; Xenium uses cluster-mode cell-state refinement; Visium
uses sr-mode virtual-cell reconstruction. SR mode preserves spot-level
evidence and does not by itself establish true sub-spot cell locations. See
[Concepts](https://revise-svc.readthedocs.io/en/latest/source/concepts.html)
for the data-shape and evidence boundaries.

## Reproduce published material

### Sim2Real-ST Benchmark

Run the bounded benchmark batch launcher once from the repository root:

```bash
BENCHMARK_MAX_JOBS=1 bash reproduce/benchmark_main.sh
```

After the Benchmark run, use the notebooks in
[`reproduce/benchmark/`](reproduce/benchmark/) for the paper analyses. See
[reproduce/README.md](reproduce/README.md) for the data layout and detailed
single-family or batch commands, or the [Benchmark documentation](https://revise-svc.readthedocs.io/en/latest/)
for the current reference.

### Real-world ST Application

Run an Application template through the Quick run section, then use the
matching notebook in [`reproduce/case/`](reproduce/case/) to analyze its H5AD.
The [Application gallery](https://revise-svc.readthedocs.io/en/latest/source/gallery.html)
and [REVISE documentation](https://revise-svc.readthedocs.io/en/latest/) provide
the preserved analysis material and current reference. The notebooks are static
snapshots, not evidence of a current rerun or biological validation.

## Documentation index

The README is the first-run guide. Each detailed rule has one canonical owner:

| Question | Canonical documentation |
| --- | --- |
| How do I run a template quickly? | [Quick Start](https://revise-svc.readthedocs.io/en/latest/source/quickstart.html) |
| What does every Application YAML field mean? | [Application Reference](https://revise-svc.readthedocs.io/en/latest/source/application-reference.html) |
| Which SVC/mode fits my data and what does it prove? | [Concepts](https://revise-svc.readthedocs.io/en/latest/source/concepts.html) |
| How do I install dependencies? | [Installation](https://revise-svc.readthedocs.io/en/latest/source/installation.html) |
| How do I migrate a previous Application request? | [Application migration](https://revise-svc.readthedocs.io/en/latest/source/application-migration.html) |
| How are Application and Benchmark routed internally? | [Architecture](https://revise-svc.readthedocs.io/en/latest/source/architecture.html) |
| Which notebooks are preserved? | [Gallery](https://revise-svc.readthedocs.io/en/latest/source/gallery.html) |
| How do I reproduce the paper workflows? | [`reproduce/README.md`](reproduce/README.md) |
| How is the repository verified? | [`tests/README.md`](tests/README.md) |

## Repository layout

```text
configs/application/     maintained Application YAML templates
reconstruct.py           Application CLI and Python entry point
revise/                  package implementation and Sim2Real-ST benchmark API
reproduce/               benchmark launchers and preserved analysis notebooks
docs/                    ReadTheDocs source
tests/                   executable contracts
```

## Citation and license

This checkout does not currently contain citation metadata; consult the
[REVISE documentation](https://revise-svc.readthedocs.io/en/latest/) for the
current project citation. REVISE is released under the [MIT License](LICENSE).
