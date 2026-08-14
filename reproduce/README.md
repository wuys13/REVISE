# Reproducing REVISE analyses

The reproduction material has two entry points:

| Workflow | Command entry | Notebook entry |
| --- | --- | --- |
| Benchmark | `python reproduce/benchmark_main.py ...` or `bash reproduce/benchmark_main.sh` from the repository root | [`benchmark/`](benchmark/) |
| Application | `python reconstruct.py ...` from the repository root, or `revise-reconstruct ...` after installation | [`case/`](case/) |

Download the paper datasets and reproduced results from
<https://zenodo.org/records/17705737>. Data are not stored in this Git
repository and are not downloaded during package installation.

Some application notebooks require optional installation groups such as
`pathway`, `cci`, or `trajectory`, as well as their corresponding external
reference resources.

## Canonical application-analysis notebooks

The maintained paper-facing analysis set is:

- [Xenium T sc-SVC](case/Xenium_sc_SVC_T.ipynb)
- [Xenium Fibroblast sc-SVC](case/Xenium_sc_SVC_Fibroblast.ipynb)
- [Xenium Monocyte sc-SVC](case/Xenium_sc_SVC_Monocyte.ipynb)
- [Visium HD sp-SVC](case/VisiumHD_sp_SVC.ipynb)

The three Xenium notebooks retain their validated execution snapshots and add
user-facing explanations for input acquisition, analysis methods, result
reading, direct observations, and paper context. The Visium HD notebook remains
an unexecuted reference workflow with an explicit bounded-audit boundary.

The separate
[Visium mouse-brain sc-SVC-sr notebook](case/sc_SVC_sr_case_Visium_mouse_brain.ipynb)
is preserved outside this four-notebook set. The three Xenium notebooks include
their analysis outputs for direct review. Visium HD preserves the analysis
workflow without claiming a full reconstruction or full-dataset execution.
