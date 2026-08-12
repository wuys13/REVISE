# REVISE paper reproduction snapshot

This branch preserves the benchmark and application workflows from
`wuys13/REVISE-legacy@reproduce` at commit
`3c546189d83abc1f9e19869bda27ef1c10d3e966`.

It is intentionally independent from the current REVISE runtime. The branch
contains one root commit and does not carry the legacy repository history.
Scientific algorithms and legacy APIs have not been rewritten to match the
current YAML-based application interface.

## Included material

- the legacy `revise` Python package;
- benchmark entry points for segmentation, bin2cell, batch effect, spot size,
  gene panel, and gene dropout experiments;
- sp-SVC and sc-SVC application entry points;
- benchmark and application notebooks with their code and Markdown retained;
- small pathway and CellPhoneDB resources used by the preserved analyses.

Notebook execution outputs and execution counters were removed to keep the
branch small. The large historical `reproduce/case/sp_SVC_case.ipynb` notebook,
project images, generated documentation, datasets, and result directories are
not included.

## Data

Download the paper data and generated results from:

<https://zenodo.org/records/17705737>

The command-line scripts expect the extracted data below the repository root:

```text
raw_data/
  Real_application/
  Sim2Real-ST/
```

The downloaded archives and generated outputs are not tracked by this branch.

## Environment

The original package metadata accepts Python 3.8 and later. Python 3.10 is the
recommended starting point for recreating the historical environment, but the
branch does not contain a fully locked dependency environment.

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

The default sc-SVC application path uses TACCO:

```bash
python -m pip install -e ".[annotation]"
```

Downstream notebooks may additionally require packages such as `omicverse`,
`scvelo`, `harmonypy`, `networkx`, or CellPhoneDB. These dependencies were not
locked in the legacy branch and should be installed only for the analysis being
run.

## Benchmark entry points

Run commands from the repository root. To launch all six benchmark families:

```bash
bash benchmark_main.sh
```

To launch one family, use its script under `reproduce/benchmark/`, for example:

```bash
bash reproduce/benchmark/benchmark_segmentation.sh
bash reproduce/benchmark/benchmark_gene_dropout.sh
```

The scripts use paper-specific filenames and directory conventions under
`raw_data/Sim2Real-ST`.

## Application entry points

Run the preserved application wrappers from the repository root:

```bash
bash application_sp_SVC_recon.sh
bash application_sc_SVC_recon.sh
```

Edit the shell variables or invoke the corresponding Python files directly
when using another patient, cell type, input filename, or output directory.

## Notebooks and historical assumptions

Benchmark notebooks are under `reproduce/benchmark/`; application notebooks
are at the repository root and under `reproduce/case/`.

Several notebooks preserve paths and working-directory assumptions from the
paper environment. Run a notebook with its own directory as the working
directory, inspect its first configuration cells, and update local data/output
paths before execution. Some cells also retain server-specific absolute paths;
these are historical workflow evidence and must be adapted locally when used.

The presence and syntax of these workflows do not establish that the current
machine has reproduced the paper results. Full reproduction requires the
external data, compatible dependencies, and a separate end-to-end run.
