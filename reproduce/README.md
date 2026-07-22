# Reproducing REVISE analyses

The reproduction material has two entry points:

| Workflow | Command entry | Notebook entry |
| --- | --- | --- |
| Benchmark | `python benchmark_main.py ...` or `bash benchmark_main.sh` from the repository root | [`benchmark/`](benchmark/) |
| Application | `revise-reconstruct ...` or `python application_reconstruct.py ...` | [`case/`](case/) |

Download the paper datasets and reproduced results from
<https://zenodo.org/records/17705737>. Data are not stored in this Git
repository and are not downloaded during package installation.

Some application notebooks require optional installation groups such as
`pathway`, `cci`, or `trajectory`, as well as their corresponding external
reference resources.
