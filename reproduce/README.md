# Reproducing REVISE analyses

This directory has two entry points: a Sim2Real-ST Benchmark launcher and the
Application reconstruction workflow. Data are not stored in this repository
and are not downloaded during package installation.

## Downloads

| Material | Download |
| --- | --- |
| Sim2Real-ST benchmark | [Zenodo](https://zenodo.org/records/21921802) |
| Reproduced benchmark results | [Zenodo](https://zenodo.org/records/21921802) |
| Real-world ST datasets | [Zenodo](https://zenodo.org/records/21921802) |
| Reproduced sp-SVC H5AD | [Zenodo](https://zenodo.org/records/18389835) |
| Reproduced sc-SVC H5AD | [Zenodo](https://zenodo.org/records/22046001) |

## Benchmark

The benchmark CLI resolves route inputs below
`<data-root>/<dataset-task>/<sample-name>`. The local checkout stores the
available inputs by route family rather than under the older
`P2CRC/cut_part1` example:

```text
raw_data/Sim2Real-ST/
├── imputation/part{1,2,3}/
│   ├── selected_xenium.h5ad
│   ├── real_sc_ref.h5ad
│   └── cut.png
├── seg/part{1,2,3}/
│   ├── selected_xenium.h5ad
│   ├── real_sc_ref.h5ad
│   ├── PM_on_cell.csv
│   └── spot_{1,2,3,4,8}/xenium_spot.h5ad
└── spot/part{1,2,3}/
    ├── selected_xenium.h5ad
    ├── simulated_xenium.h5ad
    ├── real_sc_ref_{part,all,others}.h5ad
    └── spot_{20,50,100,150,200}/
        ├── xenium_spot.h5ad
        ├── simulated_xenium_spot.h5ad
        └── cell_num_distribution.csv
```

The spot route also uses the shared `spot/PM_on_cell.csv` table; the three
`real_sc_ref_*` files live at each `spot/partN/` directory as shown above.

`gene_dropout` and `gene_panel` use the co-located `imputation/partN`
inputs. For example:

```bash
python reproduce/benchmark_main.py \
  --config configs/benchmark/gene_dropout.yaml \
  --data-root raw_data/Sim2Real-ST \
  --dataset-task imputation \
  --sample-name part1 \
  --output-root output/benchmark
```

The `seg` and `spot` trees keep shared reference files at `partN/` and put
route-specific spatial files under `spot_*`. The current benchmark input
resolver expects the reference and ground-truth roles under the selected
sample directory, so those routes need a prepared tree matching the YAML
contract before execution. The README's old `P2CRC/cut_part1` path and the
route names `segmentation`/`spot_size` do not describe the local storage tree.
This is an input-layout requirement; it does not change the benchmark
algorithm or rename the route API.

The bounded batch launcher expands all six benchmark route names. It assumes
that a matching input tree has already been prepared; it does not download,
rename, or rearrange data. Set its input/output environment variables only
after that preparation:

```bash
RAW_DATA_PATH=/data/Sim2Real-ST \
SAVE_PATH=output/benchmark \
SAMPLE_PATIENT=P2CRC \
SAMPLE_PARTS="part1" \
bash reproduce/benchmark_main.sh
```

The launcher expands the six Benchmark route families. Benchmark YAML selects
a confounding factor; it does not accept Application modes.

After a Benchmark run, use the preserved notebooks in [`benchmark/`](benchmark/)
for the associated analyses. The [REVISE documentation](https://revise-svc.readthedocs.io/en/latest/)
is the current website reference.

## Application

Choose the current Application template and run it from the repository root:

```bash
python reconstruct.py --config configs/application/VisiumHD.yaml
python reconstruct.py --config configs/application/Xenium.yaml --select-ct T
python reconstruct.py --config configs/application/Visium.yaml
```

The Xenium template is shared by T, Fibroblast, and Mono/Macro analyses. Use
the appropriate concrete ``--select-ct`` value. For the full configuration
and output contract, see the
[Application Reference](https://revise-svc.readthedocs.io/en/latest/source/application-reference.html).

After reconstruction, use the matching preserved notebook in [`case/`](case/)
to analyze the generated H5AD. The [Application gallery](https://revise-svc.readthedocs.io/en/latest/source/gallery.html)
is the website entry point.

The maintained Application profiles use `raw_data/` for reconstruction inputs
and set `output.dir` under `results/` for reconstructed H5ADs and provenance.
Case notebooks read those reconstruction artifacts from `results/` and write
their secondary analysis artifacts under `output/`; the three roots have
separate roles.

## Application notebooks

The preserved Application analysis notebooks are organized in website order:

1. [Visium HD sp-SVC](case/VisiumHD_sp_SVC.ipynb)
2. [Slide-seq mouse olfactory bulb sp-SVC](case/SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb)
3. [Slide-seq mouse colon sp-SVC](case/SlideSeq_mouse_colon_sp_SVC.ipynb)
4. [Stereo-seq zebrafish 5 hpf sp-SVC](case/StereoSeq_zebrafish_5hpf_sp_SVC.ipynb)
5. [CosMx SMI 267T_not sp-SVC](case/CosMx_SMI_267T_not_sp_SVC.ipynb)
6. [Xenium sc-SVC T cells](case/Xenium_sc_SVC_T.ipynb)
7. [Xenium sc-SVC Fibroblast](case/Xenium_sc_SVC_Fibroblast.ipynb)
8. [Xenium sc-SVC Mono/Macro](case/Xenium_sc_SVC_Monocyte.ipynb)
9. [osmFISH sc-SVC cluster](case/osmFISH_sc_SVC_cluster.ipynb)
10. [MERFISH Allen VISp sc-SVC cluster](case/MERFISH_Allen_VISp_sc_SVC_cluster.ipynb)
11. [Visium sc-SVC mouse brain](case/Visium_sc_SVC_mouse_brain.ipynb) — SR mode
12. [Visium sc-SVC P2CRC](case/Visium_sc_SVC_P2CRC.ipynb) — post-reconstruction tutorial using historical RNA-corrected results

Benchmark notebooks are under [`benchmark/`](benchmark/). The docs website
links these notebooks without executing them.

Some downstream notebooks require optional `pathway`, `cci`, or `trajectory`
installation groups and their associated external resources. Install only the
capabilities used by the analysis.

## Evidence boundary

Each notebook is a static historical snapshot. It can preserve a workflow and
displayed output, but it is not evidence that the current source has been
rerun, that all paper data were reprocessed, or that a biological conclusion
has been independently validated.
