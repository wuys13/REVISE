from __future__ import annotations

import textwrap
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent


def md(source: str):
    return nbf.v4.new_markdown_cell(textwrap.dedent(source).strip() + "\n")


def code(source: str):
    return nbf.v4.new_code_cell(textwrap.dedent(source).strip() + "\n")


def notebook(
    cells,
    *,
    kernel_name="revise-tacco-py310",
    kernel_display_name="REVISE TACCO Python 3.10",
    python_version="3.10.14",
):
    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {
            "display_name": kernel_display_name,
            "language": "python",
            "name": kernel_name,
        },
        "language_info": {"name": "python", "version": python_version},
    }
    return nb


COMMON_SETUP = r'''
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
):
    os.environ[key] = "1"
os.environ["PYTHONHASHSEED"] = "42"
os.environ["MPLBACKEND"] = "Agg"
os.environ.pop("NUMBA_DISABLE_JIT", None)

import anndata as ad
from IPython.display import Image, Markdown, display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def locate_repo_root():
    current = Path.cwd().resolve()
    for candidate in (current, *current.parents):
        if (candidate / "reconstruct.py").is_file() and (candidate / "reproduce/case/tacco").is_dir():
            return candidate
    raise RuntimeError(f"Cannot locate REVISE repository from {current}")

REPO_ROOT = locate_repo_root()
TACCO_CASE_DIR = REPO_ROOT / "reproduce/case/tacco"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(TACCO_CASE_DIR))

from notebook_utils import (
    assigned_labels,
    audit_table,
    data_root,
    independent_umap,
    independent_umaps,
    load_json,
    manifest_from_output,
    plot_independent_panels,
    sha256_file,
)

DATA_ROOT = data_root()
PREPARE_SCRIPT = TACCO_CASE_DIR / "prepare_cases.py"
FORCE_RECONSTRUCTION = os.environ.get("REVISE_TACCO_FORCE_RUN", "0") == "1"

RUN_ENV = os.environ.copy()
RUN_ENV["PYTHONPATH"] = str(REPO_ROOT)
RUN_ENV["MPLBACKEND"] = "Agg"
RUN_ENV.pop("NUMBA_DISABLE_JIT", None)

environment = pd.DataFrame(
    [
        ("python", sys.version.split()[0]),
        ("executable", sys.executable),
        ("tacco", importlib.metadata.version("tacco")),
        ("scanpy", importlib.metadata.version("scanpy")),
        ("squidpy", importlib.metadata.version("squidpy")),
        ("setuptools", importlib.metadata.version("setuptools")),
        ("scipy", importlib.metadata.version("scipy")),
        ("pyamg", importlib.metadata.version("pyamg")),
        ("data root", str(DATA_ROOT)),
    ],
    columns=["item", "value"],
)
display(environment)
assert sys.version.split()[0] == "3.10.14"
assert importlib.metadata.version("tacco") == "0.5.0"
assert importlib.metadata.version("scanpy") == "1.11.4"
assert importlib.metadata.version("squidpy") == "1.6.3"
assert importlib.metadata.version("setuptools") == "80.9.0"
assert importlib.metadata.version("pyamg") == "5.2.1"

thread_environment = {
    key: os.environ[key]
    for key in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
    )
}
display(pd.DataFrame(thread_environment.items(), columns=["thread variable", "value"]))
assert set(thread_environment.values()) == {"1"}

disk = shutil.disk_usage(DATA_ROOT)
reserved_bytes = 2 * 1024**3
display(pd.DataFrame([{
    "free_bytes": disk.free,
    "reserved_bytes": reserved_bytes,
    "usable_above_reserve_bytes": disk.free - reserved_bytes,
}]))
if disk.free < reserved_bytes:
    raise RuntimeError(
        f"Insufficient disk headroom: free={disk.free:,}, "
        f"reserved={reserved_bytes:,}, deficit={reserved_bytes - disk.free:,} bytes"
    )
'''


ALLEN_RELAXED_SETUP = r'''
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

# Keep the interactive kernel usable across compatible local environments.  The
# formal child processes receive the reproducible launch settings below.
os.environ.pop("NUMBA_DISABLE_JIT", None)

import anndata as ad
from IPython.display import Image, Markdown, display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def locate_repo_root():
    current = Path.cwd().resolve()
    for candidate in (current, *current.parents):
        if (candidate / "reconstruct.py").is_file() and (candidate / "reproduce/case/tacco").is_dir():
            return candidate
    raise RuntimeError(f"Cannot locate REVISE repository from {current}")

def installed_version(distribution):
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"

REPO_ROOT = locate_repo_root()
TACCO_CASE_DIR = REPO_ROOT / "reproduce/case/tacco"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(TACCO_CASE_DIR))

from notebook_utils import (
    assigned_labels,
    audit_table,
    data_root,
    independent_umap,
    independent_umaps,
    load_json,
    manifest_from_output,
    plot_independent_panels,
    sha256_file,
)

DATA_ROOT = data_root()
PREPARE_SCRIPT = TACCO_CASE_DIR / "prepare_cases.py"
FORCE_RECONSTRUCTION = os.environ.get("REVISE_TACCO_FORCE_RUN", "0") == "1"

project_python = REPO_ROOT / ".venv/bin/python"
FORMAL_PYTHON = Path(
    os.environ.get("REVISE_TACCO_FORMAL_PYTHON", str(project_python))
).expanduser().absolute()

RUN_ENV = os.environ.copy()
for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
):
    RUN_ENV[key] = "1"
RUN_ENV["PYTHONHASHSEED"] = "42"
RUN_ENV["PYTHONPATH"] = str(REPO_ROOT)
RUN_ENV["MPLBACKEND"] = "Agg"
RUN_ENV.pop("NUMBA_DISABLE_JIT", None)

actual_versions = {
    "python": sys.version.split()[0],
    "tacco": installed_version("tacco"),
    "scanpy": installed_version("scanpy"),
    "squidpy": installed_version("squidpy"),
}
snapshot_versions = {
    "python": "3.10.14",
    "tacco": "0.5.0",
    "scanpy": "1.11.4",
    "squidpy": "1.6.3",
}
version_audit = []
for package, actual in actual_versions.items():
    if package == "python":
        compatible = sys.version_info[:2] in {(3, 10), (3, 11)}
        status = "compatible" if compatible else "outside project range (non-blocking audit)"
        recommended = ">=3.10,<3.12; snapshot 3.10.14"
    else:
        status = "snapshot match" if actual == snapshot_versions[package] else "different (non-blocking audit)"
        recommended = snapshot_versions[package]
    version_audit.append({
        "component": package,
        "interactive value": actual,
        "recommended / saved snapshot": recommended,
        "status": status,
    })

display(pd.DataFrame(version_audit))
display(pd.DataFrame([
    {
        "role": "interactive notebook",
        "executable": sys.executable,
        "purpose": "inspection, plotting, and transparent diagnostics",
    },
    {
        "role": "formal preparation/reconstruction",
        "executable": str(FORMAL_PYTHON),
        "purpose": "prepare_cases.py and reconstruct.py subprocesses",
    },
]))
if any(row["status"] not in {"compatible", "snapshot match"} for row in version_audit):
    display(Markdown(
        "**Compatibility notice:** the interactive kernel differs from the saved "
        "snapshot. This is recorded but does not stop the notebook. Formal preparation "
        "and reconstruction use `FORMAL_PYTHON`; genuine import or API incompatibilities "
        "will still surface where the affected operation is used."
    ))

disk = shutil.disk_usage(DATA_ROOT)
free_gib = disk.free / 1024**3
display(pd.DataFrame([{
    "data root": str(DATA_ROOT),
    "free GiB": round(free_gib, 2),
    "check": "informational only",
}]))
if free_gib < 2:
    display(Markdown(
        f"**Disk notice:** only {free_gib:.2f} GiB is free. A fresh reconstruction may "
        "run out of space, but existing artifacts can still be inspected."
    ))
'''


RUN_FUNCTION = r'''
def run_reconstruction(config_path, outputs, log_name, select_ct=None):
    outputs = [Path(path) for path in outputs]
    should_run = FORCE_RECONSTRUCTION or not all(path.is_file() for path in outputs)
    command = [
        sys.executable,
        str(REPO_ROOT / "reconstruct.py"),
        "--config",
        str(config_path),
    ]
    if select_ct is not None:
        command.extend(["--select-ct", select_ct])
    display(Markdown("**Formal REVISE command**\n\n```bash\n" + " ".join(command) + "\n```"))
    log_path = DATA_ROOT / "logs" / log_name
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if should_run:
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run(
                command,
                cwd=DATA_ROOT,
                env=RUN_ENV,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=True,
            )
    else:
        for output in outputs:
            existing = ad.read_h5ad(output, backed="r")
            try:
                manifest_from_output(existing)
            finally:
                existing.file.close()
        print("Reusing existing outputs after terminal provenance and artifact validation.")
    tail = log_path.read_text(encoding="utf-8").splitlines()[-80:] if log_path.exists() else []
    print("\n".join(tail))
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing promised REVISE outputs: {missing}")
    return command, log_path
'''


ALLEN_RUN_FUNCTION = RUN_FUNCTION.replace(
    "        sys.executable,\n",
    "        str(FORMAL_PYTHON),\n",
    1,
)


def build_human_liver():
    cells = [
        md(r'''
        # AM042 human liver MERFISH: REVISE sc-SVC cluster reconstruction

        This notebook is the executable AM042 human-liver MERFISH case.  The
        prepared segmented-cell matrix, snRNA-seq reference, and the independent
        MERFISH cell-type sidecar are audited before the formal runs.  REVISE is run
        serially for the three requested broad classes: `Hepatocyte`, `HSC`, and
        `Macrophage`.

        The sidecar is a post-hoc concordance source only.  It is never passed to
        preprocessing, global anchoring, local refinement, or any output-producing
        call.  The GA diagnostic below is the exact TACCO global-anchoring operation
        selected by the application route, applied once to the complete preprocessed
        AM042 target.  It is not a second published REVISE artifact.

        Every expression-space UMAP is fitted independently on its own native matrix.
        The reference, raw-ST, and SVC panels therefore have no shared coordinate
        system and are not observation-paired measurements of the same cells.
        '''),
        code(COMMON_SETUP),
        md(r'''
        ## 1. Audit AM042 source files, prepared artifacts, and manifest

        The source manifest is part of the prepared-case contract.  This cell checks
        source-record presence and the byte/hash identities of all three prepared
        artifacts without assuming a particular cell count, gene count, or label
        cardinality.  Such dimensions are displayed from the files and are carried
        into the later provenance table rather than copied into notebook assertions.
        '''),
        code(r'''
        CASE_ROOT = DATA_ROOT / "prepared/human_liver_merfish/AM042"
        ST_PATH = CASE_ROOT / "MERFISH.h5ad"
        REFERENCE_PATH = CASE_ROOT / "snRNAseq.h5ad"
        SIDECAR_PATH = CASE_ROOT / "MERFISH_Cell_Type.csv.gz"
        SOURCE_MANIFEST_PATH = DATA_ROOT / "manifests/human_liver_merfish_AM042_source_manifest.json"

        required_paths = (ST_PATH, REFERENCE_PATH, SIDECAR_PATH, SOURCE_MANIFEST_PATH)
        missing = [str(path) for path in required_paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "AM042 prepared inputs are incomplete; missing: " + ", ".join(missing)
            )

        raw_st = ad.read_h5ad(ST_PATH)
        reference = ad.read_h5ad(REFERENCE_PATH)
        source_manifest = load_json(SOURCE_MANIFEST_PATH)
        assert source_manifest["case"] == "human_liver_merfish_AM042"
        assert source_manifest["doi"] == "10.5061/dryad.37pvmcvsg"
        assert source_manifest["dryad"]["archive_version"] == "v5"
        assert source_manifest["dryad"]["api_version"] == 277671
        assert set(("merfish", "nucseq")) <= set(source_manifest["source_files"])
        prepared_files = source_manifest["prepared_files"]

        expected_sources = {
            "merfish": {
                "name": "adata_healthy_merfish.h5ad",
                "file_id": 2903637,
                "bytes": 350386069,
                "sha256": "6b7fc5879a78b03002d94500ceffb3dcebef4cea82c29b381d720bbd650eb40c",
            },
            "nucseq": {
                "name": "adata_healthy_nucseq.h5ad",
                "file_id": 2903639,
                "bytes": 1262101004,
                "sha256": "faa8897adfc86a436053c28e7d5012a422335a00f9faf306bc09761b76127ade",
            },
        }
        source_audit = []
        for role, expected in expected_sources.items():
            record = source_manifest["source_files"][role]
            for key, value in expected.items():
                assert record[key] == value
            source_path = Path(record["path"])
            assert source_path.is_file()
            assert source_path.stat().st_size == expected["bytes"]
            assert sha256_file(source_path) == expected["sha256"]
            source_audit.append({"role": role, "path": str(source_path), **expected})

        prepared_paths = {
            "st": ST_PATH,
            "reference": REFERENCE_PATH,
            "label_sidecar": SIDECAR_PATH,
        }
        prepared_audit = []
        for role, path in prepared_paths.items():
            record = prepared_files[role]
            assert path.stat().st_size == int(record["bytes"])
            assert sha256_file(path) == record["sha256"]
            prepared_audit.append({
                "role": role,
                "path": str(path),
                "bytes": int(path.stat().st_size),
                "sha256": record["sha256"],
            })

        display(audit_table({"raw AM042 MERFISH ST": raw_st, "snRNA-seq reference": reference}))
        display(pd.DataFrame(source_audit))
        display(pd.DataFrame([source_manifest["observed"]]))
        display(pd.DataFrame(prepared_audit))

        audit = audit_table({"raw AM042 MERFISH ST": raw_st, "snRNA-seq reference": reference})
        assert audit["unique_obs_names"].all()
        assert audit["unique_var_names"].all()
        assert audit["finite_X"].all()
        assert audit["nonnegative_X"].all()
        assert raw_st.n_obs > 1 and raw_st.n_vars > 1
        assert reference.n_obs > 1 and reference.n_vars > 1
        assert "spatial" in raw_st.obsm
        assert np.isfinite(np.asarray(raw_st.obsm["spatial"])).all()
        assert "Cell_Type" not in raw_st.obs
        assert set(raw_st.obs["sample_id"].astype(str)) == {"AM042"}
        assert set(reference.obs["sample_id"].astype(str)) == {"AM042"}
        assert "transcript_counts" in raw_st.obs and "transcript_counts" in reference.obs
        assert np.array_equal(
            raw_st.obs["transcript_counts"].to_numpy(),
            np.asarray(raw_st.X.sum(axis=1)).ravel(),
        )
        assert np.array_equal(
            reference.obs["transcript_counts"].to_numpy(),
            np.asarray(reference.X.sum(axis=1)).ravel(),
        )
        assert np.array_equal(
            np.asarray(raw_st.obsm["spatial"]),
            raw_st.obs[["x", "y"]].to_numpy(dtype=float),
        )
        assert reference.var_names.equals(
            raw_st.var_names[raw_st.var_names.isin(reference.var_names)]
        )
        assert "Level1" in reference.obs and "Level2" in reference.obs

        import revise
        revise_path = Path(revise.__file__).resolve()
        assert revise_path == (REPO_ROOT / "revise/__init__.py").resolve()
        display(pd.DataFrame([{
            "runtime executable": sys.executable,
            "revise.__file__": str(revise_path),
            "MPLBACKEND": RUN_ENV["MPLBACKEND"],
        }]))

        for key in tuple(RUN_ENV):
            if key == "PYTHONPATH" or key.startswith("CONDA_"):
                RUN_ENV.pop(key)
        RUN_ENV["PYTHONPATH"] = str(REPO_ROOT)
        assert not any(key.startswith("CONDA_") for key in RUN_ENV)
        assert RUN_ENV["PYTHONPATH"] == str(REPO_ROOT)

        selected_types = ("Hepatocyte", "HSC", "Macrophage")
        reference_broad = reference.obs["Level1"].astype(str)
        absent = [label for label in selected_types if label not in set(reference_broad)]
        if absent:
            raise ValueError(f"AM042 reference is missing selected broad labels: {absent}")
        selected_reference = reference[reference_broad.isin(selected_types)].copy()
        assert selected_reference.n_obs > 0
        display(reference_broad.value_counts().rename("reference cells").to_frame())
        display(
            selected_reference.obs.groupby("Level1", observed=True)["Level2"]
            .nunique()
            .rename("reference subtypes")
            .to_frame()
        )

        sidecar_table = pd.read_csv(
            SIDECAR_PATH,
            compression="gzip",
            index_col="cell_id",
        )
        if sidecar_table.empty:
            raise ValueError("AM042 MERFISH cell-type sidecar is empty")
        st_names = pd.Index(raw_st.obs_names.astype(str))
        sidecar_table.index = sidecar_table.index.astype(str)
        SIDECAR_LABEL_COLUMN = "Cell_Type"
        assert sidecar_table.columns.tolist() == [SIDECAR_LABEL_COLUMN]
        assert sidecar_table.index.is_unique
        assert sidecar_table.index.equals(st_names)
        display(
            sidecar_table[SIDECAR_LABEL_COLUMN]
            .value_counts()
            .rename("MERFISH sidecar cells")
            .to_frame()
        )
        display(pd.DataFrame([{
            "sidecar rows": int(sidecar_table.shape[0]),
            "sidecar columns": int(sidecar_table.shape[1]),
            "cell-id column": sidecar_table.index.name,
            "label column": SIDECAR_LABEL_COLUMN,
            "raw-ST ID overlap": int(len(st_names)),
        }]))
        '''),
        md(r'''
        ## 2. Run the formal sc-SVC applications serially

        Each selected broad class gets its own `spatial.h5ad` and `expr.h5ad` output.
        The output manifest is checked immediately after each run, including terminal
        success, publication artifact hashes, configuration identity, selected-cell
        type, and the declared output paths.  Running these applications serially keeps
        the three provenance records and their resource use unambiguous.
        '''),
        code(RUN_FUNCTION),
        code(r'''
        CONFIG_PATH = TACCO_CASE_DIR / "configs/MERFISH_human_liver_sc_SVC_cluster.yaml"
        SELECTED_TYPES = ("Hepatocyte", "HSC", "Macrophage")
        expected_config_sha256 = sha256_file(CONFIG_PATH)
        output_paths = {}
        run_records = []

        for broad in SELECTED_TYPES:
            output_dir = DATA_ROOT / "results/human_liver_merfish/AM042" / broad
            spatial_path = output_dir / "spatial.h5ad"
            expr_path = output_dir / "expr.h5ad"
            run_reconstruction(
                CONFIG_PATH,
                [spatial_path, expr_path],
                f"MERFISH_human_liver_sc_SVC_cluster_{broad}.log",
                select_ct=broad,
            )
            spatial_output = ad.read_h5ad(spatial_path)
            expr_output = ad.read_h5ad(expr_path)
            spatial_manifest_path, spatial_manifest = manifest_from_output(spatial_output)
            expr_manifest_path, expr_manifest = manifest_from_output(expr_output)
            assert spatial_manifest_path == expr_manifest_path
            assert spatial_manifest == expr_manifest
            spatial_metadata = spatial_output.uns["revise_reconstruction"]
            expr_metadata = expr_output.uns["revise_reconstruction"]
            assert spatial_metadata["output_role"] == "spatial"
            assert expr_metadata["output_role"] == "expression"
            assert spatial_manifest["application_config"]["source_sha256"] == expected_config_sha256
            effective_request = spatial_manifest["application_config"]["effective_request"]
            assert effective_request["application_route"] == "sc-SVC"
            assert effective_request["application_mode"] == "cluster"
            assert effective_request["selected_cell_type"] == broad
            declared_paths = spatial_manifest["application_config"]["output_paths"]
            assert Path(declared_paths["spatial"]).resolve() == spatial_path.resolve()
            assert Path(declared_paths["expression"]).resolve() == expr_path.resolve()
            assert spatial_manifest["input_identities"]
            assert spatial_output.n_obs > 1 and expr_output.n_obs > 1
            assert spatial_output.var_names.equals(expr_output.var_names)
            assert spatial_output.obs_names.is_unique and expr_output.obs_names.is_unique
            assert {"Level1", "Level2", "SVC_cluster"} <= set(spatial_output.obs)
            assert {"Level1", "Level2", "SVC_cluster"} <= set(expr_output.obs)
            assert "spatial" in spatial_output.obsm
            assert np.isfinite(np.asarray(spatial_output.obsm["spatial"])).all()
            assert "SVC_cluster" in expr_output.obsm
            assert "Cell_Type" not in spatial_output.obs
            assert "Cell_Type" not in expr_output.obs
            output_paths[broad] = {
                "spatial_path": spatial_path,
                "expr_path": expr_path,
                "spatial": spatial_output,
                "expr": expr_output,
                "manifest_path": spatial_manifest_path,
                "manifest": spatial_manifest,
            }
            run_records.append({
                "broad class": broad,
                "run status": spatial_manifest["run"]["status"],
                "spatial shape": str(spatial_output.shape),
                "expression shape": str(expr_output.shape),
                "manifest": str(spatial_manifest_path),
            })

        display(pd.DataFrame(run_records))
        '''),
        md(r'''
        ## 3. Exact global anchoring, GA UMAP, and complete spatial map

        The following call mirrors the `sc-SVC:cluster` route itself: its GA solver and
        TACCO parameters are read from the route registry, rather than re-entered as a
        notebook-only approximation.  The complete preprocessed target is retained so
        that all broad labels are visible in the diagnostic UMAP and measured-space
        map.  The selected outputs are checked against this exact GA by cell identity
        and spatial coordinates.
        '''),
        code(r'''
        import scanpy as sc
        from reconstruct import preprocess_data
        from revise.application.config import compile_application_config, load_application_yaml
        from revise.backend.kernels.ot import OTKernel
        from revise.config import ENGINE_DEFAULTS, ROUTES
        from revise.utils.deterministic import set_global_seed

        config_source, config_document = load_application_yaml(CONFIG_PATH)
        app_config = compile_application_config(
            config_document,
            source=config_source,
            cwd=DATA_ROOT,
        )
        assert app_config.svc_type == "sc-SVC"
        assert app_config.mode == "cluster"
        assert app_config.ot_method == "tacco"
        assert app_config.broad_column == "Level1"
        assert app_config.subtype_column == "Level2"
        assert app_config.seed == 42
        assert app_config.spatial_min_transcript_counts == 15
        assert app_config.spatial_min_cell_counts == 1
        assert app_config.reference_min_transcript_counts is None
        assert app_config.reference_min_cell_counts == 1
        assert app_config.local_refinement_alpha == 0.2
        assert app_config.local_refinement_resolutions == (0.6, 0.7, 0.8)

        full_spatial, full_reference = preprocess_data(
            raw_st.copy(),
            reference.copy(),
            app_config,
        )
        assert full_spatial.n_obs > 1 and full_spatial.n_vars > 1
        assert full_reference.n_obs > 1 and full_reference.n_vars > 1
        assert full_spatial.obs_names.is_unique and full_reference.obs_names.is_unique
        assert "spatial" in full_spatial.obsm
        assert np.isfinite(np.asarray(full_spatial.obsm["spatial"])).all()
        assert set(SELECTED_TYPES) <= set(full_reference.obs[app_config.broad_column].astype(str))

        route_overrides = ROUTES["application"]["sc-SVC:cluster"].overrides
        ga_method = route_overrides["ot"]["ga"]["solver"]
        tacco_parameters = route_overrides["sc"]["tacco_annotate"]
        confidence_key = ENGINE_DEFAULTS["columns"]["confidence_col"]
        unknown_key = ENGINE_DEFAULTS["columns"]["unknown_key"]
        assert ga_method == "tacco"
        assert tacco_parameters == {"multi_center": 1, "lamb": 0.001}

        set_global_seed(app_config.seed)
        ga_spatial = OTKernel.annotate(
            full_spatial,
            full_reference,
            method=ga_method,
            annotation_key=app_config.broad_column,
            confidence_key=confidence_key,
            unknown_key=unknown_key,
            **tacco_parameters,
        )
        broad_order = list(pd.unique(
            full_reference.obs[app_config.broad_column].astype(str)
        ))
        ga_labels = ga_spatial.obs[app_config.broad_column].astype(str)
        posterior = ga_spatial.obsm[app_config.broad_column]
        assert ga_spatial.shape == full_spatial.shape
        assert ga_spatial.obs_names.equals(full_spatial.obs_names)
        assert posterior.shape == (ga_spatial.n_obs, len(broad_order))
        assert list(map(str, posterior.columns)) == broad_order
        assert set(ga_labels).issubset(set(broad_order))
        assert set(SELECTED_TYPES) <= set(ga_labels)
        observed_broad_order = [
            label for label in broad_order if label in set(ga_labels)
        ]

        selected_spatial = {}
        for broad in SELECTED_TYPES:
            spatial_output = output_paths[broad]["spatial"].copy()
            observed = spatial_output.obs[app_config.broad_column].astype(str)
            assert set(observed.unique()) == {broad}
            assert spatial_output.obs_names.isin(full_spatial.obs_names).all()
            expected_names = ga_spatial.obs_names[ga_labels == broad]
            assert spatial_output.obs_names.equals(expected_names)
            assert np.array_equal(
                np.asarray(spatial_output.obsm["spatial"]),
                np.asarray(ga_spatial[expected_names].obsm["spatial"]),
            )
            selected_spatial[broad] = spatial_output

        selected_names = pd.Index([])
        for broad in SELECTED_TYPES:
            selected_names = selected_names.append(selected_spatial[broad].obs_names)
        assert selected_names.is_unique

        ga_umap = independent_umap(
            ga_spatial,
            categorical_labels=ga_labels,
            source_name="complete_ga",
        )
        broad_palette = plt.get_cmap("tab20").colors
        broad_colors = {
            label: broad_palette[index % len(broad_palette)]
            for index, label in enumerate(observed_broad_order)
        }
        ga_umap.obs[app_config.broad_column] = pd.Categorical(
            ga_umap.obs[app_config.broad_column].astype(str),
            categories=observed_broad_order,
            ordered=True,
        )
        fig, ax = plt.subplots(figsize=(8.2, 5.8))
        sc.pl.umap(
            ga_umap,
            color=app_config.broad_column,
            palette=[broad_colors[label] for label in observed_broad_order],
            size=10,
            title="AM042 MERFISH complete GA broad labels",
            legend_fontsize=7,
            frameon=False,
            ax=ax,
            show=False,
        )
        GA_UMAP_PATH = DATA_ROOT / "figures/MERFISH_human_liver_AM042_GA_umap.png"
        GA_UMAP_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(GA_UMAP_PATH, dpi=220, bbox_inches="tight")
        display(Image(filename=str(GA_UMAP_PATH)))
        plt.close(fig)
        display(pd.DataFrame([ga_umap.uns["independent_umap_contract"]]))

        background = np.asarray(ga_spatial.obsm["spatial"])
        ga_counts = ga_labels.value_counts().reindex(observed_broad_order)
        fig, ax = plt.subplots(figsize=(8.2, 8.2), constrained_layout=True)
        point_colors = np.asarray([broad_colors[label] for label in ga_labels])
        ax.scatter(
            background[:, 0], background[:, 1], s=2.0, c=point_colors,
            alpha=0.85, linewidths=0, rasterized=True,
        )
        from matplotlib.lines import Line2D
        ax.legend(
            handles=[
                Line2D(
                    [0], [0], marker="o", linestyle="", color=broad_colors[label],
                    markersize=5, label=f"{label} ({ga_counts.loc[label]:,})",
                )
                for label in observed_broad_order
            ],
            bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
            fontsize=7, title="GA Level1", title_fontsize=8,
        )
        ax.set_title("AM042 MERFISH spatial localization of complete GA")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        GA_SPATIAL_PATH = DATA_ROOT / "figures/MERFISH_human_liver_AM042_GA_spatial.png"
        fig.savefig(GA_SPATIAL_PATH, dpi=220, bbox_inches="tight")
        display(Image(filename=str(GA_SPATIAL_PATH)))
        plt.close(fig)
        display(pd.DataFrame({
            "figure": [GA_UMAP_PATH, GA_SPATIAL_PATH],
            "role": ["complete GA expression UMAP", "complete GA measured-space map"],
        }))
        '''),
        md(r'''
        ## 4. Post-hoc sidecar concordance and subtype-versus-SVC spatial panels

        Concordance is reported only where sidecar cell IDs intersect a formal output.
        It is a descriptive comparison of independently produced labels; it is not a
        supervised input, a validation split, or an observation-paired biological
        estimate.  The three spatial rows retain the same measured-space limits while
        comparing each reference `Level2` subtype view with the formal `SVC_cluster`
        view for Hepatocyte, HSC, and Macrophage.
        '''),
        code(r'''
        from sklearn.metrics import adjusted_rand_score

        def _posthoc_concordance(left, right, *, pair_name):
            left = pd.Series(left, dtype="string")
            right = pd.Series(right, dtype="string")
            common = left.index.intersection(right.index)
            if len(common) == 0:
                return {
                    "comparison": pair_name,
                    "n_common": 0,
                    "adjusted_rand_index": np.nan,
                    "table": pd.DataFrame(),
                }
            left_common = left.reindex(common).astype(str)
            right_common = right.reindex(common).astype(str)
            ari = (
                adjusted_rand_score(left_common, right_common)
                if len(pd.unique(left_common)) > 1 and len(pd.unique(right_common)) > 1
                else np.nan
            )
            return {
                "comparison": pair_name,
                "n_common": int(len(common)),
                "adjusted_rand_index": float(ari) if np.isfinite(ari) else np.nan,
                "table": pd.crosstab(left_common, right_common),
            }

        sidecar_concordance = []
        sidecar_labels = sidecar_table[SIDECAR_LABEL_COLUMN]
        sidecar_broad_labels = sidecar_labels.map(
            lambda label: (
                "Hepatocyte" if str(label).startswith("Hep") else
                "HSC" if str(label).startswith("HSC") else
                "Macrophage" if str(label).startswith("Mac") else
                str(label)
            )
        )
        sidecar_concordance.append(
            _posthoc_concordance(
                sidecar_broad_labels,
                ga_labels,
                pair_name="broad-mapped sidecar versus exact GA Level1",
            )
        )
        for broad in SELECTED_TYPES:
            spatial_output = selected_spatial[broad]
            sidecar_concordance.append(
                _posthoc_concordance(
                    sidecar_labels,
                    spatial_output.obs["SVC_cluster"],
                    pair_name=f"sidecar versus {broad} SVC_cluster",
                )
            )
        display(pd.DataFrame([
            {key: value for key, value in record.items() if key != "table"}
            for record in sidecar_concordance
        ]))
        for record in sidecar_concordance:
            if not record["table"].empty:
                display(Markdown(f"**{record['comparison']}**"))
                display(record["table"])

        spatial_x = (float(background[:, 0].min()), float(background[:, 0].max()))
        spatial_y = (float(background[:, 1].min()), float(background[:, 1].max()))
        fig, axes = plt.subplots(
            len(SELECTED_TYPES), 2, figsize=(13, 4.8 * len(SELECTED_TYPES)),
            squeeze=False, constrained_layout=True,
        )
        for row, broad in enumerate(SELECTED_TYPES):
            spatial_output = selected_spatial[broad]
            coords = np.asarray(spatial_output.obsm["spatial"])
            for column, key in enumerate(("Level2", "SVC_cluster")):
                assert key in spatial_output.obs
                ax = axes[row, column]
                labels = spatial_output.obs[key].astype(str)
                categories = labels.value_counts().index.tolist()
                palette = plt.get_cmap("tab20").colors
                for index, label in enumerate(categories):
                    selected = labels.to_numpy() == label
                    ax.scatter(
                        coords[selected, 0], coords[selected, 1], s=1.2,
                        color=palette[index % len(palette)], label=label,
                        linewidths=0, rasterized=True,
                    )
                title = "reference Level2 subtype" if key == "Level2" else "REVISE SVC_cluster"
                ax.set_title(f"{broad} | {title}")
                ax.set_xlim(spatial_x)
                ax.set_ylim(spatial_y)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.legend(
                    bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False,
                    fontsize=7, markerscale=3,
                )
        SVC_SPATIAL_PATH = DATA_ROOT / "figures/MERFISH_human_liver_AM042_subtype_vs_SVC_spatial.png"
        fig.savefig(SVC_SPATIAL_PATH, dpi=220, bbox_inches="tight")
        display(Image(filename=str(SVC_SPATIAL_PATH)))
        plt.close(fig)
        '''),
        md(r'''
        ## 5. Independent reference, raw-ST, and SVC UMAPs

        The raw-ST panel is colored with the exact GA `Level1` assignment, while the
        SVC panel is colored by the reconstructed `SVC_cluster` labels.  These are
        intentionally three independently fitted UMAPs, not a joint embedding and not
        a shared-coordinate comparison.  The expression and spatial artifacts also
        have different observation roles; the SVC expression-side rows are not the
        same cells as the raw MERFISH spatial observations.
        '''),
        code(r'''
        raw_for_umap = raw_st[ga_spatial.obs_names, :].copy()
        raw_for_umap_labels = ga_labels.reindex(raw_for_umap.obs_names)
        assert raw_for_umap_labels.notna().all()
        svc_expression = ad.concat(
            [output_paths[broad]["expr"] for broad in SELECTED_TYPES],
            join="inner",
            label="_selected_type",
            index_unique="::",
            merge="same",
        )
        assert "SVC_cluster" in svc_expression.obs
        svc_label_column = "SVC_cluster"
        svc_expression_labels = svc_expression.obs[svc_label_column].astype(str)
        reference_for_umap = full_reference.copy()
        reference_for_umap_labels = reference_for_umap.obs[app_config.broad_column].astype(str)
        embeddings = independent_umaps(
            {
                "reference": reference_for_umap,
                "raw_st": raw_for_umap,
                "svc": svc_expression,
            },
            categorical_labels={
                "reference": reference_for_umap_labels,
                "raw_st": raw_for_umap_labels,
                "svc": svc_expression_labels,
            },
        )
        display(pd.DataFrame([
            embeddings[source].uns["independent_umap_contract"]
            for source in ("reference", "raw_st", "svc")
        ]))
        fig = plot_independent_panels(
            embeddings,
            ["reference", "raw_st", "svc"],
            titles={
                "reference": "AM042 snRNA-seq reference",
                "raw_st": "Raw AM042 MERFISH ST | exact GA Level1",
                "svc": "REVISE AM042 sc-SVC | SVC_cluster",
            },
        )
        INDEPENDENT_UMAP_PATH = DATA_ROOT / "figures/MERFISH_human_liver_AM042_independent_umap.png"
        INDEPENDENT_UMAP_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(INDEPENDENT_UMAP_PATH, dpi=220, bbox_inches="tight")
        display(Image(filename=str(INDEPENDENT_UMAP_PATH)))
        plt.close(fig)

        paired_overlap = raw_for_umap.obs_names.intersection(svc_expression.obs_names)
        display(pd.DataFrame([{
            "raw ST observations": int(raw_for_umap.n_obs),
            "SVC expression observations": int(svc_expression.n_obs),
            "raw-ST/SVC ID overlap": int(len(paired_overlap)),
            "interpretation": "not observation-paired; compare panels descriptively only",
        }]))
        '''),
        md(r'''
        ## Interpretation boundary

        The sidecar concordance and all figures are downstream diagnostics.  They do
        not turn one AM042 specimen into independent biological replication.  In
        particular, the reference, raw-ST, and SVC expression objects are not matched
        cell-by-cell; independent UMAP coordinates cannot be compared across panels,
        and spatial subtype/SVC panels describe the formal reconstruction rather than
        an independent validation set.
        '''),
    ]
    return notebook(cells)


# ---------------------------------------------------------------------------
# Root TACCO gallery cases
# ---------------------------------------------------------------------------
#
# The original AM042 builder above is retained byte-for-byte because its
# nested notebook is a historical, unexecuted source snapshot.  The four root
# gallery cases use the smaller contract below: an explicit data-root
# environment variable, serial preparation/reconstruction, and no output
# reuse.  Keeping these definitions separate also prevents a gallery cleanup
# from changing the human-liver source snapshot.

TARGET_COMMON_SETUP = r'''
import os
from pathlib import Path
import subprocess
import sys

for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS",
):
    os.environ[key] = "1"
os.environ["PYTHONHASHSEED"] = "42"
os.environ["MPLBACKEND"] = "Agg"
SEED = 42

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(SEED)

REPO_ROOT = Path.cwd().resolve()
TACCO_CASE_DIR = REPO_ROOT / "reproduce/case/tacco"

DATA_ROOT = Path(os.environ["REVISE_TACCO_DATA_ROOT"]).expanduser().resolve()
os.environ["REVISE_TACCO_DATA_ROOT"] = str(DATA_ROOT)
PREPARE_SCRIPT = TACCO_CASE_DIR / "prepare_cases.py"
RUN_ENV = os.environ.copy()
RUN_ENV["REVISE_TACCO_DATA_ROOT"] = str(DATA_ROOT)
RUN_ENV["PYTHONPATH"] = str(REPO_ROOT)
RUN_ENV["MPLBACKEND"] = "Agg"
'''

def build_olfactory():
    cells = [
        md(r'''
        # Case study: Slide-seq mouse olfactory bulb sp-SVC

        ## Notebook Guide

        Prepare TACCO's OB1 Slide5 puck, run the formal REVISE sp-SVC application,
        and inspect independently fitted reference, raw-ST, and SVC UMAPs.  Raw-ST
        colors are the formal `svc.obs["Level1"]` labels aligned by observation ID.
        Set `REVISE_TACCO_DATA_ROOT` first; preparation also requires
        `REVISE_READSEURAT_PYTHON` for the RDS converter.
        The PNG in the figure cell is a historical snapshot; it was not re-run here.
        Phase 1 snapshot: this streamlined code has not been rerun.

        **Question.** Does the formal sp-SVC route yield a readable same-case
        reconstruction?  **Method.** Prepare, reconstruct, and fit source-local UMAPs.
        **Direct observation.** The retained image shows those source-local panels.
        **Interpretation boundary.** It is not independent biological validation.
        '''),
        code(TARGET_COMMON_SETUP),
        md("## 1. Prepare the original TACCO inputs"),
        code(r'''
        subprocess.run(
            [sys.executable, str(PREPARE_SCRIPT), "olfactory"],
            cwd=REPO_ROOT,
            env=RUN_ENV,
            check=True,
        )
        CASE_ROOT = DATA_ROOT / "prepared/slideseq_olfactory_bulb"
        ST_PATH = CASE_ROOT / "OB1_Slide5.h5ad"
        REFERENCE_PATH = CASE_ROOT / "GSE121891_reference.h5ad"
        raw_st = ad.read_h5ad(ST_PATH)
        reference = ad.read_h5ad(REFERENCE_PATH)
        '''),
        md("## 2. Run the formal TACCO-backed REVISE sp-SVC route"),
        code(r'''
        CONFIG_PATH = TACCO_CASE_DIR / (
            "configs/SlideSeq_mouse_olfactory_bulb_sp_SVC.yaml"
        )
        command = [
            sys.executable,
            str(REPO_ROOT / "reconstruct.py"),
            "--config",
            str(CONFIG_PATH),
        ]
        subprocess.run(command, cwd=DATA_ROOT, env=RUN_ENV, check=True)
        '''),
        code(r'''
        SVC_PATH = DATA_ROOT / (
            "results/slideseq_olfactory_bulb/"
            "SlideSeq_mouse_olfactory_bulb_sp_SVC.h5ad"
        )
        svc = ad.read_h5ad(SVC_PATH)
        '''),
        md(r'''
        ## 3. Independent expression-space UMAPs

        Each source is fitted independently with direct Scanpy preprocessing and
        UMAP.  UMAP coordinates are descriptive and are not comparable across
        panels.  The raw-ST panel is colored from the formal SVC `Level1` labels,
        not from a second notebook-side global anchoring run.
        '''),
        code(r'''
        import scanpy as sc

        raw_for_umap = raw_st[
            np.asarray(raw_st.X.sum(axis=1)).ravel() >= 50
        ].copy()
        umap_objects = {
            "reference": reference.copy(),
            "raw_st": raw_for_umap.copy(),
            "svc": svc.copy(),
        }
        svc_level1 = svc.obs["Level1"].astype(str)
        umap_objects["raw_st"].obs["Level1"] = (
            svc_level1.loc[raw_for_umap.obs_names].to_numpy()
        )
        for adata in umap_objects.values():
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=30, random_state=SEED)
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
            sc.tl.umap(adata, random_state=SEED)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), squeeze=False)
        for ax, (name, title) in zip(
            axes.ravel(),
            (
                ("reference", "GSE121891 single-cell reference"),
                ("raw_st", "Raw OB1 Slide5 ST | formal SVC Level1"),
                ("svc", "REVISE reconstructed sp-SVC"),
            ),
        ):
            sc.pl.umap(
                umap_objects[name],
                color="Level1",
                title=title,
                size=4,
                frameon=False,
                legend_loc="none",
                ax=ax,
                show=False,
            )
        plt.tight_layout()
        plt.show()
        '''),
        md(r'''
        ## Interpretation boundary

        The three panels have independent geometry and therefore do not support
        cross-panel distance or orientation claims.  Historical figure outputs are
        retained for the gallery while this source remains explicitly unexecuted.
        '''),
    ]
    return notebook(cells)


def build_colon():
    cells = [
        md(r'''
        # Case study: Slide-seq mouse colon sp-SVC

        ## Notebook Guide

        Download the five authenticated SCP2038 files, prepare the normal puck
        `2020-09-14_Puck_200701_21`, run the formal REVISE sp-SVC route, and inspect
        independent UMAPs.  Raw-ST colors come from formal `svc.obs["Level1"]`
        labels aligned by observation ID.  The PNG is a historical snapshot and was
        not re-run here; missing portal authentication fails naturally in preparation.
        Set `REVISE_TACCO_DATA_ROOT` before running.
        Phase 1 snapshot: this streamlined code has not been rerun.

        **Question.** What does the authenticated normal puck produce under sp-SVC?
        **Method.** Prepare the exact portal files and run the configured route.
        **Direct observation.** The retained image is the historical UMAP snapshot.
        **Interpretation boundary.** Missing authentication stops preparation; no
        substitute data are inferred.
        '''),
        code(TARGET_COMMON_SETUP),
        md(r'''
        ## 1. Prepare the authenticated SCP2038 source

        The preparation command requires the exact authenticated Single Cell Portal
        artifacts; no substitute data are accepted.
        '''),
        code(r'''
        subprocess.run(
            [sys.executable, str(PREPARE_SCRIPT), "colon"],
            cwd=REPO_ROOT,
            env=RUN_ENV,
            check=True,
        )
        '''),
        code(r'''
        CASE_ROOT = DATA_ROOT / "prepared/slideseq_mouse_colon"
        ST_PATH = CASE_ROOT / "normal_Puck_200701_21.h5ad"
        REFERENCE_PATH = CASE_ROOT / "normal_scRNAseq.h5ad"
        raw_st = ad.read_h5ad(ST_PATH)
        reference = ad.read_h5ad(REFERENCE_PATH)
        '''),
        md("## 2. Run the formal TACCO-backed REVISE sp-SVC route"),
        code(r'''
        CONFIG_PATH = TACCO_CASE_DIR / "configs/SlideSeq_mouse_colon_sp_SVC.yaml"
        command = [
            sys.executable,
            str(REPO_ROOT / "reconstruct.py"),
            "--config",
            str(CONFIG_PATH),
        ]
        subprocess.run(command, cwd=DATA_ROOT, env=RUN_ENV, check=True)
        '''),
        code(r'''
        SVC_PATH = DATA_ROOT / (
            "results/slideseq_mouse_colon/SlideSeq_mouse_colon_sp_SVC.h5ad"
        )
        svc = ad.read_h5ad(SVC_PATH)
        '''),
        md(r'''
        ## 3. Independent expression-space UMAPs

        Each source is fitted independently with direct Scanpy preprocessing and
        UMAP.  The raw-ST panel uses the formal
        SVC `Level1` labels, with index alignment, and does not repeat complete GA.
        The selected configuration uses `strength = 0.0` and joint-graph `alpha = 0.8`.
        '''),
        code(r'''
        import scanpy as sc

        raw_for_umap = raw_st[
            np.asarray(raw_st.X.sum(axis=1)).ravel() >= 50
        ].copy()
        umap_objects = {
            "reference": reference.copy(),
            "raw_st": raw_for_umap.copy(),
            "svc": svc.copy(),
        }
        svc_level1 = svc.obs["Level1"].astype(str)
        umap_objects["raw_st"].obs["Level1"] = (
            svc_level1.loc[raw_for_umap.obs_names].to_numpy()
        )
        for adata in umap_objects.values():
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=30, random_state=SEED)
            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
            sc.tl.umap(adata, random_state=SEED)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), squeeze=False)
        for ax, (name, title) in zip(
            axes.ravel(),
            (
                ("reference", "Normal colon single-cell reference"),
                ("raw_st", "Raw normal Slide-seq ST | formal SVC Level1"),
                ("svc", "REVISE sp-SVC | strength 0, graph alpha 0.8"),
            ),
        ):
            sc.pl.umap(
                umap_objects[name],
                color="Level1",
                title=title,
                size=4,
                frameon=False,
                legend_loc="none",
                ax=ax,
                show=False,
            )
        plt.tight_layout()
        plt.show()
        '''),
        md(r'''
        ## Interpretation boundary

        UMAP orientation, scale, and distance are source-local.  This gallery snapshot
        records a formal reconstruction workflow; it is not an independent biological
        validation of the SCP2038 case.
        '''),
    ]
    return notebook(cells)


def build_osmfish():
    cells = [
        md(r'''
        # Case study: osmFISH TonT sc-SVC cluster reconstruction

        ## Notebook Guide

        Prepare the TACCO osmFISH TonT segmented-cell matrix, run the formal REVISE
        sc-SVC cluster route for Pyramidal and Inhibitory cells, and inspect a complete
        12-class GA overview plus selected spatial panels.  `TonT_labels.csv.gz` is a
        sidecar for context only and is not passed into preparation, GA, or either
        reconstruction.  The three PNGs are historical snapshots; they were not
        re-run here.
        Set `REVISE_TACCO_DATA_ROOT` before running.
        Phase 1 snapshot: this streamlined code has not been rerun.

        **Question.** How do selected cluster reconstructions sit within a complete GA
        overview?  **Method.** Run both formal routes and one in-memory complete GA.
        **Direct observation.** The retained images show broad and selected spatial
        views.  **Interpretation boundary.** The sidecar is not reconstruction truth.
        '''),
        code(TARGET_COMMON_SETUP),
        md("## 1. Prepare the TonT segmented-cell input"),
        code(r'''
        subprocess.run(
            [sys.executable, str(PREPARE_SCRIPT), "osmfish"],
            cwd=REPO_ROOT,
            env=RUN_ENV,
            check=True,
        )
        CASE_ROOT = DATA_ROOT / "prepared/osmfish"
        ST_PATH = CASE_ROOT / "TonT_segmented_cells.h5ad"
        REFERENCE_PATH = CASE_ROOT / "osmFISH_reference.h5ad"
        raw_st = ad.read_h5ad(ST_PATH)
        reference = ad.read_h5ad(REFERENCE_PATH)
        '''),
        md(r'''
        The molecule-derived TonT sidecar is deliberately excluded from the REVISE
        input contract.  It is not used for preprocessing, global anchoring, local
        refinement, output labels, or the complete-GA overview.
        '''),
        md("## 2. Run the Pyramidal and Inhibitory sc-SVC cluster routes"),
        code(r'''
        CONFIG_PATH = TACCO_CASE_DIR / "configs/osmFISH_sc_SVC_cluster.yaml"
        '''),
        code(r'''
        output_paths = {}
        for broad in ("Pyramidal", "Inhibitory"):
            output_dir = DATA_ROOT / "results/osmfish" / broad
            spatial_path = output_dir / "spatial.h5ad"
            expr_path = output_dir / "expr.h5ad"
            command = [
                sys.executable,
                str(REPO_ROOT / "reconstruct.py"),
                "--config",
                str(CONFIG_PATH),
                "--select-ct",
                broad,
            ]
            subprocess.run(command, cwd=DATA_ROOT, env=RUN_ENV, check=True)
            spatial_output = ad.read_h5ad(spatial_path)
            expr_output = ad.read_h5ad(expr_path)
            output_paths[broad] = {
                "spatial": spatial_output,
                "expr": expr_output,
            }
        '''),
        md(r'''
        The two published artifacts share reconstructed phenotype semantics but have
        different observation roles.  The overview below repeats the exact complete
        TACCO GA stage in memory solely to show all broad classes; it does not create
        another application artifact.
        '''),
        md("## 3. Complete GA overview and selected spatial maps"),
        code(r'''
        import scanpy as sc
        import pandas as pd
        from reconstruct import preprocess_data
        from revise.application.config import compile_application_config, load_application_yaml
        from revise.backend.kernels.ot import OTKernel
        from revise.config import ENGINE_DEFAULTS, ROUTES
        from revise.utils.deterministic import set_global_seed

        config_source, config_document = load_application_yaml(CONFIG_PATH)
        app_config = compile_application_config(
            config_document,
            source=config_source,
            cwd=DATA_ROOT,
        )
        full_spatial, full_reference = preprocess_data(
            raw_st.copy(), reference.copy(), app_config
        )
        route_overrides = ROUTES["application"]["sc-SVC:cluster"].overrides
        ga_method = route_overrides["ot"]["ga"]["solver"]
        tacco_parameters = route_overrides["sc"]["tacco_annotate"]
        confidence_key = ENGINE_DEFAULTS["columns"]["confidence_col"]
        unknown_key = ENGINE_DEFAULTS["columns"]["unknown_key"]
        set_global_seed(SEED)
        ga_spatial = OTKernel.annotate(
            full_spatial,
            full_reference,
            method=ga_method,
            annotation_key=app_config.broad_column,
            confidence_key=confidence_key,
            unknown_key=unknown_key,
            **tacco_parameters,
        )
        broad_order = list(pd.unique(
            full_reference.obs[app_config.broad_column].astype(str)
        ))
        ga_labels = ga_spatial.obs[app_config.broad_column].astype(str)
        selected_spatial = {
            broad: output_paths[broad]["spatial"]
            for broad in ("Pyramidal", "Inhibitory")
        }
        ga_umap = ga_spatial.copy()
        sc.pp.normalize_total(ga_umap, target_sum=1e4)
        sc.pp.log1p(ga_umap)
        sc.pp.pca(ga_umap, n_comps=30, random_state=SEED)
        sc.pp.neighbors(ga_umap, n_neighbors=15, n_pcs=30, random_state=SEED)
        sc.tl.umap(ga_umap, random_state=SEED)
        broad_palette = plt.get_cmap("tab20").colors
        broad_colors = {
            label: broad_palette[index % len(broad_palette)]
            for index, label in enumerate(broad_order)
        }
        fig, ax = plt.subplots(figsize=(8.2, 5.8))
        sc.pl.umap(
            ga_umap,
            color=app_config.broad_column,
            size=10,
            title="All broad cell types after complete GA",
            frameon=False,
            legend_loc="none",
            ax=ax,
            show=False,
        )
        coordinates = np.asarray(ga_spatial.obsm["spatial"])
        fig, ax = plt.subplots(figsize=(7.5, 8.5), constrained_layout=True)
        point_colors = np.asarray([broad_colors[label] for label in ga_labels])
        ax.scatter(
            coordinates[:, 0], coordinates[:, 1], s=2.0, c=point_colors,
            alpha=0.85, linewidths=0, rasterized=True,
        )
        ax.set_title("Spatial localization of complete GA broad classes")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        x_limits = (float(coordinates[:, 0].min()), float(coordinates[:, 0].max()))
        y_limits = (float(coordinates[:, 1].min()), float(coordinates[:, 1].max()))
        fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
        for row, broad in enumerate(("Pyramidal", "Inhibitory")):
            spatial_output = selected_spatial[broad]
            coords = np.asarray(spatial_output.obsm["spatial"])
            for column, key in enumerate(("Level2", "SVC_cluster")):
                ax = axes[row, column]
                labels = spatial_output.obs[key].astype(str)
                categories = labels.value_counts().index.tolist()
                palette = plt.get_cmap("tab20").colors
                for index, label in enumerate(categories):
                    selected = labels.to_numpy() == label
                    ax.scatter(
                        coords[selected, 0], coords[selected, 1], s=1.2,
                        color=palette[index % len(palette)],
                        linewidths=0, rasterized=True,
                    )
                ax.set_title(f"{broad} | {key}")
                ax.set_xlim(x_limits)
                ax.set_ylim(y_limits)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xticks([])
                ax.set_yticks([])
        plt.show()
        '''),
        md(r'''
        ## Interpretation boundary

        The complete GA UMAP and measured-space map contain all broad classes.  The
        selected SVC panels are separate application outputs and are not an independent
        benchmark; all spatial and expression UMAP geometry is source-local.
        '''),
    ]
    return notebook(cells)


def build_allen_visp_merfish():
    cells = [
        md(r'''
        # Case study: Allen VISp MERFISH sc-SVC cluster reconstruction

        ## Notebook Guide

        Prepare the Allen VISp MERFISH table and same-species Tasic cortex reference,
        run formal REVISE sc-SVC cluster applications for three broad classes, and
        inspect the complete GA overview plus reconstructed clusters.  The three PNG
        outputs are historical snapshots; this notebook source was not re-run.
        Set `REVISE_TACCO_DATA_ROOT` before running.
        Phase 1 snapshot: this streamlined code has not been rerun.

        **Question.** What broad and cluster structure is visible in the Allen VISp
        case?  **Method.** Prepare, run three serial applications, and repeat complete
        GA in memory.  **Direct observation.** The retained PNGs show the overview and
        reconstructed panels.  **Interpretation boundary.** The reference is a
        same-species cortex match, not cell-, animal-, or experiment-paired truth.
        '''),
        code(TARGET_COMMON_SETUP),
        md(r'''
        ## 1. Prepare source identities

        Preparation runs the explicit Allen converter through `prepare_cases.py` and
        then reads the prepared H5AD objects.  The archive is paired with an external
        Tasic cortex reference by species and region only; no cell, animal, specimen,
        or experiment pairing is claimed.
        '''),
        code(r'''
        subprocess.run(
            [sys.executable, str(PREPARE_SCRIPT), "allen_visp"],
            cwd=REPO_ROOT,
            env=RUN_ENV,
            check=True,
        )
        CASE_ROOT = DATA_ROOT / "prepared/allen_visp_merfish"
        ST_PATH = CASE_ROOT / "MERFISH.h5ad"
        REFERENCE_PATH = CASE_ROOT / "Tasic2018_sc_mouse_cortex.h5ad"
        raw_st = ad.read_h5ad(ST_PATH)
        reference = ad.read_h5ad(REFERENCE_PATH)
        '''),
        md("## 2. Run three formal TACCO-backed sc-SVC applications serially"),
        code(r'''
        CONFIG_PATH = TACCO_CASE_DIR / (
            "configs/MERFISH_Allen_VISp_sc_SVC_cluster.yaml"
        )
        SELECTED_TYPES = ("Glutamatergic", "GABAergic", "Non-Neuronal")
        '''),
        code(r'''
        output_paths = {}
        for broad in SELECTED_TYPES:
            output_dir = DATA_ROOT / "results/allen_visp_merfish" / broad
            spatial_path = output_dir / "spatial.h5ad"
            expr_path = output_dir / "expr.h5ad"
            command = [
                sys.executable,
                str(REPO_ROOT / "reconstruct.py"),
                "--config",
                str(CONFIG_PATH),
                "--select-ct",
                broad,
            ]
            subprocess.run(command, cwd=DATA_ROOT, env=RUN_ENV, check=True)
            spatial_output = ad.read_h5ad(spatial_path)
            expr_output = ad.read_h5ad(expr_path)
            output_paths[broad] = {
                "spatial": spatial_output,
                "expr": expr_output,
            }
        '''),
        md(r'''
        ## 3. Complete GA overview: expression UMAP and measured-space map

        The exact route-level TACCO GA operation is repeated once on the complete
        preprocessed target for an overview of all broad labels.  This in-memory
        diagnostic is not an additional published artifact.
        '''),
        code(r'''
        import scanpy as sc
        import pandas as pd
        from reconstruct import preprocess_data
        from revise.application.config import compile_application_config, load_application_yaml
        from revise.backend.kernels.ot import OTKernel
        from revise.config import ENGINE_DEFAULTS, ROUTES
        from revise.utils.deterministic import set_global_seed

        config_source, config_document = load_application_yaml(CONFIG_PATH)
        app_config = compile_application_config(
            config_document, source=config_source, cwd=DATA_ROOT
        )
        full_spatial, full_reference = preprocess_data(
            raw_st.copy(), reference.copy(), app_config
        )
        route_overrides = ROUTES["application"]["sc-SVC:cluster"].overrides
        ga_method = route_overrides["ot"]["ga"]["solver"]
        tacco_parameters = route_overrides["sc"]["tacco_annotate"]
        confidence_key = ENGINE_DEFAULTS["columns"]["confidence_col"]
        unknown_key = ENGINE_DEFAULTS["columns"]["unknown_key"]
        set_global_seed(SEED)
        ga_spatial = OTKernel.annotate(
            full_spatial,
            full_reference,
            method=ga_method,
            annotation_key=app_config.broad_column,
            confidence_key=confidence_key,
            unknown_key=unknown_key,
            **tacco_parameters,
        )
        broad_order = list(pd.unique(
            full_reference.obs[app_config.broad_column].astype(str)
        ))
        ga_labels = ga_spatial.obs[app_config.broad_column].astype(str)
        selected_spatial = {
            broad: output_paths[broad]["spatial"] for broad in SELECTED_TYPES
        }
        observed_broad = [
            label for label in broad_order if label in set(ga_labels)
        ]
        broad_palette = plt.get_cmap("tab10").colors
        broad_colors = {
            label: broad_palette[index % len(broad_palette)]
            for index, label in enumerate(observed_broad)
        }
        ga_umap = ga_spatial.copy()
        sc.pp.normalize_total(ga_umap, target_sum=1e4)
        sc.pp.log1p(ga_umap)
        sc.pp.pca(ga_umap, n_comps=30, random_state=SEED)
        sc.pp.neighbors(ga_umap, n_neighbors=15, n_pcs=30, random_state=SEED)
        sc.tl.umap(ga_umap, random_state=SEED)
        fig, ax = plt.subplots(figsize=(7.8, 5.6))
        sc.pl.umap(
            ga_umap,
            color=app_config.broad_column,
            size=12,
            title="Allen VISp MERFISH | complete GA Level1",
            frameon=False,
            legend_loc="none",
            ax=ax,
            show=False,
        )
        coordinates = np.asarray(ga_spatial.obsm["spatial"])
        fig, ax = plt.subplots(figsize=(7.8, 7.8), constrained_layout=True)
        point_colors = np.asarray([broad_colors[label] for label in ga_labels])
        ax.scatter(
            coordinates[:, 0], coordinates[:, 1], s=4, c=point_colors,
            alpha=0.9, linewidths=0, rasterized=True,
        )
        ax.set_title("Allen VISp MERFISH | complete GA spatial map")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.show()
        '''),
        md(r'''
        ## 4. Transferred Tasic subtypes versus reconstructed SVC clusters

        `Level2` is the transferred Tasic subtype and `SVC_cluster` is the formal
        reconstructed cluster.  The six panels use shared MERFISH coordinates for
        the three selected broad classes and are a within-reconstruction comparison.
        '''),
        code(r'''
        coordinates = np.asarray(ga_spatial.obsm["spatial"])
        x_limits = (float(coordinates[:, 0].min()), float(coordinates[:, 0].max()))
        y_limits = (float(coordinates[:, 1].min()), float(coordinates[:, 1].max()))
        fig, axes = plt.subplots(
            len(SELECTED_TYPES),
            2,
            figsize=(13, 4.2 * len(SELECTED_TYPES)),
            squeeze=False,
            constrained_layout=True,
        )
        for row, broad in enumerate(SELECTED_TYPES):
            spatial_output = selected_spatial[broad]
            coords = np.asarray(spatial_output.obsm["spatial"])
            for column, key in enumerate(("Level2", "SVC_cluster")):
                ax = axes[row, column]
                labels = spatial_output.obs[key].astype(str)
                categories = labels.value_counts().index.tolist()
                palette = plt.get_cmap("tab10").colors
                for index, label in enumerate(categories):
                    selected = labels.to_numpy() == label
                    ax.scatter(
                        coords[selected, 0],
                        coords[selected, 1],
                        s=4,
                        color=palette[index % len(palette)],
                        linewidths=0,
                        rasterized=True,
                    )
                ax.set_title(f"{broad} | {key}")
                ax.set_xlim(x_limits)
                ax.set_ylim(y_limits)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xticks([])
                ax.set_yticks([])
        plt.show()
        '''),
        md(r'''
        ## Interpretation boundary

        These measured-space panels compare `Level2` and `SVC_cluster` within one
        reconstruction.  They are not author-label validation, independent biological
        validation, or evidence of cell-level ground truth.
        '''),
    ]
    return notebook(cells, kernel_name="python3", kernel_display_name="Python 3 (REVISE compatible)", python_version="3.10")



def main():
    case_root = HERE.parent
    outputs = {
        case_root / "SlideSeq_mouse_olfactory_bulb_sp_SVC.ipynb": build_olfactory(),
        case_root / "SlideSeq_mouse_colon_sp_SVC.ipynb": build_colon(),
        case_root / "osmFISH_sc_SVC_cluster.ipynb": build_osmfish(),
        HERE / "MERFISH_human_liver_sc_SVC_cluster.ipynb": build_human_liver(),
        case_root / "MERFISH_Allen_VISp_sc_SVC_cluster.ipynb": build_allen_visp_merfish(),
    }
    for path, nb in outputs.items():
        nbf.write(nb, path)
        print(path)


if __name__ == "__main__":
    main()
