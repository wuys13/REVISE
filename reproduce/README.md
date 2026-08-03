# Reproducing REVISE analyses

Everything under `reproduce/` is 1.x historical reproduction material, not
current 2.0 output. Historical notebooks and carrier filenames remain unchanged;
use the current CLI contract for new 2.0 reconstructions.

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
