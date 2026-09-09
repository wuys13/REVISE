from __future__ import annotations

import argparse
import gc
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import tacco as tc
from scipy import io as scipy_io
from scipy import sparse

TACCO_EXAMPLES_COMMIT = "ed61ddc584be72217fe83d33b8995589264efd50"
RESERVED_FREE_BYTES = 2 * 1024**3

OLFACTORY_URLS = {
    "rds_gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5173nnn/"
        "GSM5173929/suppl/GSM5173929_OB1_Slide5.rds.gz"
    ),
    "expression_gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121891/suppl/"
        "GSE121891_OB_6_runs.raw.dge.csv.gz"
    ),
    "metadata_gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121891/suppl/"
        "GSE121891_OB_metaData_seurat.csv.gz"
    ),
}
OLFACTORY_EXPECTED_BYTES = {
    "rds_gz": 445_764_676,
    "expression_gz": 62_652_087,
    "metadata_gz": 901_809,
}

COLON_FILENAMES = (
    "Slideseq_raw.mtx.gz",
    "Slideseq_raw_cells.tsv.gz",
    "Slideseq_raw_genes.tsv.gz",
    "scRNAseq.h5ad",
    "singlecellportal_metadata.tsv.gz",
)

OSMFISH_ORIGINAL_URLS = {
    "loom": "http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom",
    "molecules": (
        "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/"
        "data/mRNA_coords_raw_counting.hdf5"
    ),
    "segmentation": (
        "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/osmFISH/"
        "data/polyT_seg.pkl"
    ),
}
OSMFISH_EFFECTIVE_URLS = {
    "loom": (
        "https://d24h2xsgaj29mf.cloudfront.net/raw/"
        "osmfish_codeluppi_2018_nat-methods_somatosensory-cortex/"
        "osmFISH_SScortex_mouse_all_cells.loom"
    ),
    "molecules": (
        "https://d24h2xsgaj29mf.cloudfront.net/raw/"
        "osmfish_codeluppi_2018_nat-methods_somatosensory-cortex/"
        "mRNA_coords_raw_counting.hdf5"
    ),
    "segmentation": (
        "https://d24h2xsgaj29mf.cloudfront.net/raw/"
        "osmfish_codeluppi_2018_nat-methods_somatosensory-cortex/"
        "polyT_seg.pkl"
    ),
}
OSMFISH_EXPECTED_BYTES = {
    "loom": 1_120_947,
    "molecules": 31_643_024,
    "segmentation": 1_461_350_426,
}

HUMAN_LIVER_DOI = "10.5061/dryad.37pvmcvsg"
HUMAN_LIVER_DRYAD_LANDING_URL = (
    f"https://datadryad.org/stash/dataset/doi:{HUMAN_LIVER_DOI}"
)
HUMAN_LIVER_DRYAD_VERSION = "v5"
HUMAN_LIVER_DRYAD_API_VERSION = 277671
HUMAN_LIVER_SOURCE_FILES = {
    "merfish": {
        "name": "adata_healthy_merfish.h5ad",
        "file_id": 2_903_637,
        "bytes": 350_386_069,
        "size": 350_386_069,
        "sha256": "6b7fc5879a78b03002d94500ceffb3dcebef4cea82c29b381d720bbd650eb40c",
        "url": "https://datadryad.org/downloads/file_stream/2903637",
    },
    "nucseq": {
        "name": "adata_healthy_nucseq.h5ad",
        "file_id": 2_903_639,
        "bytes": 1_262_101_004,
        "size": 1_262_101_004,
        "sha256": "faa8897adfc86a436053c28e7d5012a422335a00f9faf306bc09761b76127ade",
        "url": "https://datadryad.org/downloads/file_stream/2903639",
    },
}

ALLEN_VISP_ARCHIVE = {
    "name": "merfish_spatialdata_0.7.3a2.dev27+gde3ed360b.zip",
    "bytes": 53_359_643,
    "sha256": "978d837b8c9ca25ace028c164c4d8d9db1dc82b30850d6dcb28f2d4ae803bfe3",
    "spatialdata_version": "0.7.3a2.dev27+gde3ed360b",
}
ALLEN_VISP_TASIC_REFERENCE = {
    "name": "sc_mouse_cortex.h5ad",
    "bytes": 3_254_624_804,
    "sha256": "3e0a26e1af06c1ea8f53a808ee683bf950de8cc03ee48bd291f95eeca6056aac",
    "url": "https://exampledata.scverse.org/squidpy/sc_mouse_cortex.h5ad",
    "doi": "10.1038/s41586-018-0654-5",
    "squidpy_registry_commit": "be17fcf6afddbf06011429a6cf43b70cb9fac9b4",
}
ALLEN_VISP_SOURCE_EVIDENCE = {
    "converter_repository": "https://github.com/giovp/spatialdata-sandbox",
    "converter_commit_audited": "565b712307fd52c85b9f6c93c7990e3ca4a4cd84",
    "converter_script": "merfish/download.py",
    "spacejam_repository": "https://github.com/spacetx-spacejam/data",
    "spacejam_commit_audited": "45df7f4929db424cc003c906c7943b0545442c32",
    "spacejam_dataset": "MERFISH (Allen: VISp)",
}


def data_root() -> Path:
    return Path(os.environ["REVISE_TACCO_DATA_ROOT"]).expanduser().resolve()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_headroom(root: Path, incoming_bytes: int, label: str) -> None:
    root.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(root).free
    required = int(incoming_bytes) + RESERVED_FREE_BYTES
    if free < required:
        deficit = required - free
        raise RuntimeError(
            f"Insufficient disk for {label}: free={free:,}, required={required:,}, "
            f"deficit={deficit:,} bytes"
        )


def _probe_url(url: str) -> dict:
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request) as response:
            status = int(response.status)
            final_url = response.geturl()
            headers = response.headers
    except urllib.error.HTTPError as error:
        status = int(error.code)
        final_url = error.geturl()
        headers = error.headers
    content_length = headers.get("Content-Length")
    return {
        "http_status": status,
        "http_final_url": final_url,
        "http_content_length": int(content_length) if content_length else None,
        "http_etag": headers.get("ETag"),
        "http_last_modified": headers.get("Last-Modified"),
    }


def _download(url: str, target: Path, *, expected_bytes: int) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    probe = _probe_url(url)
    if probe["http_status"] != 200:
        raise RuntimeError(f"Source URL is not HTTP 200: {url} -> {probe}")
    content_length = probe["http_content_length"]
    if content_length is not None and content_length != expected_bytes:
        raise RuntimeError(
            f"Source Content-Length mismatch for {url}: "
            f"{content_length:,} != {expected_bytes:,}"
        )
    if target.exists():
        observed = target.stat().st_size
        if observed != expected_bytes:
            raise RuntimeError(
                f"Existing download has wrong size: {target} ({observed:,} != "
                f"{expected_bytes:,})"
            )
        return {
            "url": url,
            "path": str(target),
            "bytes": observed,
            "sha256": sha256_file(target),
            "downloaded_now": False,
            **probe,
        }

    _require_headroom(target.parent, expected_bytes, target.name)
    temporary = target.with_suffix(target.suffix + ".part")
    if temporary.exists():
        raise RuntimeError(f"Incomplete download must be inspected first: {temporary}")
    digest = hashlib.sha256()
    written = 0
    with urllib.request.urlopen(url) as response, temporary.open("wb") as handle:
        while True:
            chunk = response.read(8 * 1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
            digest.update(chunk)
            written += len(chunk)
    if written != expected_bytes:
        raise RuntimeError(
            f"Downloaded byte count mismatch for {url}: {written:,} != {expected_bytes:,}"
        )
    os.replace(temporary, target)
    return {
        "url": url,
        "path": str(target),
        "bytes": written,
        "sha256": digest.hexdigest(),
        "downloaded_now": True,
        **probe,
    }


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _json_compatible(value):
    if isinstance(value, np.ndarray):
        return [_json_compatible(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item) for item in value]
    return value


def _atomic_h5ad(path: Path, adata: ad.AnnData) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.h5ad")
    if temporary.exists():
        raise RuntimeError(f"Temporary H5AD must be inspected first: {temporary}")
    adata.write_h5ad(temporary, compression="gzip")
    check = ad.read_h5ad(temporary, backed="r")
    observed_shape = tuple(check.shape)
    check.file.close()
    if observed_shape != tuple(adata.shape):
        raise RuntimeError(f"H5AD shape changed during write: {observed_shape} != {adata.shape}")
    os.replace(temporary, path)


def _matrix_min(matrix) -> float:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return float(values.min()) if values.size else 0.0


def _validate_counts(adata: ad.AnnData, name: str) -> None:
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError(f"{name} has non-unique observation or variable names")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or _matrix_min(adata.X) < 0:
        raise ValueError(f"{name} contains non-finite or negative expression values")
    if "spatial" in adata.obsm and not np.isfinite(
        np.asarray(adata.obsm["spatial"])
    ).all():
        raise ValueError(f"{name} contains non-finite spatial coordinates")


def _adata_identity(path: Path, adata: ad.AnnData) -> dict:
    _validate_counts(adata, path.name)
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "unique_obs_names": bool(adata.obs_names.is_unique),
        "unique_var_names": bool(adata.var_names.is_unique),
        "has_spatial": "spatial" in adata.obsm,
    }


def _runtime_identity() -> dict:
    versions = {}
    for package in (
        "tacco",
        "scanpy",
        "squidpy",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "pyamg",
        "setuptools",
    ):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    tracked_diff = subprocess.run(
        ["git", "diff", "--binary"],
        cwd=repo_root(),
        check=True,
        capture_output=True,
    ).stdout
    case_files = subprocess.run(
        [
            "git",
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "--",
            "reproduce/case/tacco",
        ],
        cwd=repo_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    dirty_state = hashlib.sha256(tracked_diff)
    for relative in sorted(case_files):
        path = repo_root() / relative
        if not path.is_file():
            continue
        dirty_state.update(relative.encode("utf-8") + b"\0")
        dirty_state.update(path.read_bytes())
    return {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "packages": versions,
        "repo_head": git_head,
        "repo_tracked_diff_sha256": hashlib.sha256(tracked_diff).hexdigest(),
        "repo_dirty_diff_sha256": dirty_state.hexdigest(),
        "repo_case_files": sorted(case_files),
        "tacco_examples_commit": TACCO_EXAMPLES_COMMIT,
        "thread_environment": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMBA_NUM_THREADS",
                "PYTHONHASHSEED",
            )
        },
    }


def _broad_osmfish(label: str) -> str:
    text = str(label)
    lowered = text.lower()
    for prefix, broad in (
        ("pyramidal", "Pyramidal"),
        ("inhibitory", "Inhibitory"),
        ("astrocyte", "Astrocyte"),
        ("oligodendrocyte", "Oligodendrocyte"),
        ("endothelial", "Endothelial"),
    ):
        if lowered.startswith(prefix):
            return broad
    return text


def _write_manifest(case: str, payload: dict) -> Path:
    path = data_root() / "manifests" / f"{case}_source_manifest.json"
    _atomic_json(
        path,
        {
            "case": case,
            "runtime": _runtime_identity(),
            **payload,
        },
    )
    return path


def _prepare_olfactory_reference(expression_path: Path, metadata_path: Path) -> ad.AnnData:
    metadata = pd.read_csv(metadata_path, compression="gzip", index_col=0)
    metadata.index = metadata.index.astype(str)

    cluster2type = {
        "OEC": [f"OEC{i + 1}" for i in range(5)],
        "N": [f"N{i + 1}" for i in range(16)],
        "Astro": [f"Astro{i + 1}" for i in range(3)],
        "EC": [f"EC{i + 1}" for i in range(2)],
        "MicroG": [f"MicroG{i + 1}" for i in range(3)],
        "Mural": [f"Mural{i + 1}" for i in range(2)],
        "Mes": [f"Mes{i + 1}" for i in range(2)],
    }
    type2long = {
        "olfactory ensheathing cell-based (Sox10+)": "OEC",
        "neuronal (Syt1+/Tubb3+)": "N",
        "astrocytic (Gfap+)": "Astro",
        "endothelial (Slco1c1+)": "EC",
        "microglia (Aif1+/Siglech+)": "MicroG",
        "myelinating-oligodendrocyte-based (Mag+)": "MyOligo",
        "mural (Pdgfrb+)": "Mural",
        "mesenchymal": "Mes",
        "monocyte (Aif1+/Cd74+)": "Mono",
        "macrophage (Aif1+/Cd52+)": "Mφ",
        "oligodendrocyte-precursor-based (Olig2+)": "OPC",
        "red blood cell (Aif1+/Hba-a1+)": "RBCs",
    }
    tc.utils.merge_annotation(metadata, "ClusterName", cluster2type, "type")
    tc.utils.merge_annotation(metadata, "type", type2long, "long")

    expression_cells = pd.Index(
        pd.read_csv(expression_path, compression="gzip", nrows=0).columns.astype(str)
    )
    cells = metadata.index.tolist()
    missing = sorted(set(cells) - set(expression_cells))
    if missing:
        raise ValueError(f"GSE121891 expression is missing metadata cells: {missing[:5]}")
    extra = expression_cells.difference(pd.Index(cells))
    if len(expression_cells) != 52_549 or len(extra) != 1_123:
        raise ValueError(
            "Unexpected GSE121891 cell axes: "
            f"expression={len(expression_cells):,}, metadata={len(cells):,}, "
            f"expression_only={len(extra):,}"
        )

    blocks = []
    genes: list[str] = []
    reader = pd.read_csv(
        expression_path,
        compression="gzip",
        index_col=0,
        usecols=cells,
        chunksize=128,
    )
    for chunk in reader:
        if chunk.columns.astype(str).tolist() != cells:
            chunk = chunk.loc[:, cells]
        genes.extend(chunk.index.astype(str).tolist())
        values = chunk.to_numpy(dtype=np.int16, copy=False)
        blocks.append(sparse.csr_matrix(values))
    expression = sparse.vstack(blocks, format="csr").T.tocsr()
    reference = ad.AnnData(
        X=expression,
        obs=metadata.loc[cells].copy(),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    reference.obs["Level1"] = reference.obs["type"].astype(str)
    reference.obs["Level2"] = reference.obs["ClusterName"].astype(str)
    reference.obs["transcript_counts"] = np.asarray(reference.X.sum(axis=1)).ravel()
    return reference


def prepare_olfactory() -> Path:
    root = data_root()
    raw = root / "raw" / "slideseq_olfactory_bulb"
    prepared = root / "prepared" / "slideseq_olfactory_bulb"
    scratch = root / "scratch" / "slideseq_olfactory_bulb"
    raw.mkdir(parents=True, exist_ok=True)
    prepared.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)

    raw_paths = {
        "rds_gz": raw / "GSM5173929_OB1_Slide5.rds.gz",
        "expression_gz": raw / "GSE121891_OB_6_runs.raw.dge.csv.gz",
        "metadata_gz": raw / "GSE121891_OB_metaData_seurat.csv.gz",
    }
    source_records = {
        key: _download(
            OLFACTORY_URLS[key],
            path,
            expected_bytes=OLFACTORY_EXPECTED_BYTES[key],
        )
        for key, path in raw_paths.items()
    }

    st_path = prepared / "OB1_Slide5.h5ad"
    if not st_path.exists():
        _require_headroom(scratch, 3 * 1024**3, "decompressed olfactory RDS")
        rds_path = scratch / "GSM5173929_OB1_Slide5.rds"
        converted_path = scratch / "GSM5173929_OB1_Slide5_readseurat.h5ad"
        with gzip.open(raw_paths["rds_gz"], "rb") as source, rds_path.open("wb") as target:
            shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
        converter_python = Path(
            os.environ["REVISE_READSEURAT_PYTHON"]
        ).expanduser().resolve()
        converter_code = """
import importlib.metadata
import platform
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from readseurat.rdata import read_rds

if tuple(map(int, platform.python_version_tuple())) < (3, 11, 0):
    raise RuntimeError(f"The isolated RDS converter requires Python >=3.11, got {platform.python_version()}")
if importlib.metadata.version("readseurat") != "0.1.0":
    raise RuntimeError("The isolated RDS converter requires readseurat==0.1.0")

seurat = read_rds(sys.argv[1])
version = tuple(np.asarray(seurat.version[0], dtype=int).tolist())
if version != (3, 2, 0):
    raise RuntimeError(f"Expected the audited Seurat 3.2.0 object, got {version}")
if set(seurat.assays) != {"RNA", "SCT"}:
    raise RuntimeError(f"Unexpected Seurat assays: {list(seurat.assays)}")
if "Spatial" not in seurat.reductions:
    raise RuntimeError("The audited Seurat object lacks its Spatial reduction")

rna = seurat.assays["RNA"]
counts = rna.counts
meta = getattr(seurat, "meta.data")
features = getattr(rna, "meta.features")
spatial = getattr(seurat.reductions["Spatial"], "cell.embeddings")
cells = pd.Index(meta.index.astype(str), name="cell")
genes = pd.Index(features.index.astype(str), name="gene")
spatial_cells = pd.Index(spatial.coords[spatial.dims[0]].values.astype(str), name="cell")
spatial_columns = spatial.coords[spatial.dims[1]].values.astype(str).tolist()

if not sparse.isspmatrix_csc(counts):
    raise RuntimeError(f"Expected CSC RNA counts, got {type(counts)!r}")
if counts.shape != (22_170, 55_086) or counts.shape != (len(genes), len(cells)):
    raise RuntimeError(f"Unexpected RNA counts shape: {counts.shape}")
if meta.shape != (55_086, 11) or features.shape != (22_170, 0):
    raise RuntimeError(f"Unexpected Seurat metadata shapes: meta={meta.shape}, features={features.shape}")
if tuple(spatial.shape) != (55_086, 2) or spatial_columns != ["Spatial_1", "Spatial_2"]:
    raise RuntimeError(f"Unexpected Spatial reduction: shape={spatial.shape}, columns={spatial_columns}")
if not cells.is_unique or not genes.is_unique or not spatial_cells.is_unique:
    raise RuntimeError("The audited Seurat axes must be unique")
if not cells.equals(spatial_cells):
    raise RuntimeError("Seurat meta.data and Spatial cell axes are not identically ordered")
if not np.isfinite(counts.data).all() or (counts.data < 0).any():
    raise RuntimeError("RNA counts contain non-finite or negative values")
coordinates = np.asarray(spatial.values, dtype=np.float64)
if not np.isfinite(coordinates).all():
    raise RuntimeError("Spatial coordinates contain non-finite values")

converted = ad.AnnData(
    X=counts.T.tocsr().astype(np.float32),
    obs=pd.DataFrame(index=cells),
    var=pd.DataFrame(index=genes),
)
converted.obsm["X_Spatial"] = coordinates
converted.uns["source_conversion"] = {
    "seurat_version": ".".join(map(str, version)),
    "assay": "RNA",
    "slot": "counts",
    "gene_axis": "RNA@meta.features.index",
    "cell_axis": "meta.data.index",
    "spatial_axis": "reductions$Spatial@cell.embeddings",
    "readseurat_version": importlib.metadata.version("readseurat"),
    "python": platform.python_version(),
}
converted.write_h5ad(sys.argv[2], compression="gzip")
print(converted.shape, list(converted.obsm), converted.uns["source_conversion"])
"""
        subprocess.run(
            [str(converter_python), "-c", converter_code, str(rds_path), str(converted_path)],
            check=True,
        )
        st = ad.read_h5ad(converted_path)
        if tuple(st.shape) != (55_086, 22_170):
            raise ValueError(f"Unexpected olfactory RDS shape: {st.shape}")
        if "X_Spatial" not in st.obsm:
            raise KeyError("readseurat output lacks the expected X_Spatial reduction")
        coords = np.asarray(st.obsm["X_Spatial"], dtype=float)[:, :2]
        st.obs = pd.DataFrame(index=st.obs_names.copy())
        st.obs["x"] = coords[:, 0]
        st.obs["y"] = coords[:, 1]
        st.obs["transcript_counts"] = np.asarray(st.X.sum(axis=1)).ravel()
        st.obsm.clear()
        st.obsm["spatial"] = coords
        _atomic_h5ad(st_path, st)
        rds_path.unlink()
        converted_path.unlink()
        del st
        gc.collect()

    reference_path = prepared / "GSE121891_reference.h5ad"
    if not reference_path.exists():
        reference = _prepare_olfactory_reference(
            raw_paths["expression_gz"], raw_paths["metadata_gz"]
        )
        if tuple(reference.shape) != (51_426, 18_560):
            raise ValueError(f"Unexpected olfactory reference shape: {reference.shape}")
        if reference.obs["Level2"].nunique() != 38:
            raise ValueError("Expected 38 olfactory ClusterName labels")
        if reference.obs["Level1"].nunique() != 12:
            raise ValueError("Expected 12 olfactory broad labels")
        _atomic_h5ad(reference_path, reference)
        del reference
        gc.collect()

    st = ad.read_h5ad(st_path)
    reference = ad.read_h5ad(reference_path)
    retained_beads = np.asarray(st.X.sum(axis=1)).ravel() >= 50
    filtered_beads = int(retained_beads.sum())
    raw_overlap = st.var_names.intersection(reference.var_names, sort=False)
    spatial_expressed = (
        np.asarray((st[retained_beads, raw_overlap].X != 0).sum(axis=0)).ravel() >= 1
    )
    reference_expressed = (
        np.asarray((reference[:, raw_overlap].X != 0).sum(axis=0)).ravel() >= 1
    )
    effective_overlap = raw_overlap[spatial_expressed & reference_expressed]
    if filtered_beads != 44_311:
        raise ValueError(f"Expected 44,311 olfactory beads after counts>=50, got {filtered_beads}")
    if len(raw_overlap) != 17_410:
        raise ValueError(f"Expected 17,410 raw shared olfactory genes, got {len(raw_overlap)}")
    if len(effective_overlap) != 17_402:
        raise ValueError(
            "Expected 17,402 effective olfactory genes after the configured filters, "
            f"got {len(effective_overlap)}"
        )
    manifest = _write_manifest(
        "slideseq_olfactory_bulb",
        {
            "source_files": source_records,
            "transformations": [
                "readseurat 0.1.0 conversion of the RNA counts and Spatial reduction",
                "Level1=type; Level2=ClusterName",
                "REVISE config applies counts>=50 to spatial beads",
            ],
            "observed": {
                "raw_st_shape": list(st.shape),
                "reference_shape": list(reference.shape),
                "beads_after_counts_50": filtered_beads,
                "raw_shared_genes": len(raw_overlap),
                "shared_genes": len(effective_overlap),
                "level1_labels": int(reference.obs["Level1"].nunique()),
                "level2_labels": int(reference.obs["Level2"].nunique()),
                "rds_conversion": dict(st.uns.get("source_conversion", {})),
            },
            "prepared_files": {
                "st": _adata_identity(st_path, st),
                "reference": _adata_identity(reference_path, reference),
            },
        },
    )
    print(manifest)
    return manifest


def prepare_colon() -> Path:
    root = data_root()
    raw = root / "raw" / "slideseq_mouse_colon"
    prepared = root / "prepared" / "slideseq_mouse_colon"
    raw.mkdir(parents=True, exist_ok=True)
    missing = [name for name in COLON_FILENAMES if not (raw / name).is_file()]
    if missing:
        raise RuntimeError(
            "SCP2038_AUTH_REQUIRED: sign in at "
            "https://singlecell.broadinstitute.org/single_cell/study/SCP2038, "
            f"download {list(COLON_FILENAMES)!r}, and place them in {str(raw)!r}. "
            f"Missing now: {missing!r}"
        )

    source_records = {
        name: {
            "portal": "https://singlecell.broadinstitute.org/single_cell/study/SCP2038",
            "path": str(raw / name),
            "bytes": (raw / name).stat().st_size,
            "sha256": sha256_file(raw / name),
        }
        for name in COLON_FILENAMES
    }
    with gzip.open(raw / "Slideseq_raw.mtx.gz", "rb") as handle:
        matrix = scipy_io.mmread(handle).T.tocsr().astype(np.float32)
    cells = pd.read_csv(
        raw / "Slideseq_raw_cells.tsv.gz", sep="\t", header=None
    ).to_numpy().ravel().astype(str)
    genes = pd.read_csv(
        raw / "Slideseq_raw_genes.tsv.gz", sep="\t", header=None
    ).to_numpy().ravel().astype(str)
    spatial = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
    )
    reference = ad.read_h5ad(raw / "scRNAseq.h5ad")
    metadata = pd.read_csv(
        raw / "singlecellportal_metadata.tsv.gz",
        sep="\t",
        skiprows=[1],
        index_col="NAME",
    )
    spatial.obs["x"] = metadata.loc[spatial.obs_names, "x_spatial"].to_numpy(dtype=float)
    spatial.obs["y"] = metadata.loc[spatial.obs_names, "y_spatial"].to_numpy(dtype=float)
    spatial.obs["sample"] = metadata.loc[spatial.obs_names, "puck"].astype(str).to_numpy()
    spatial.obs["State"] = metadata.loc[
        spatial.obs_names, "disease__ontology_label"
    ].astype(str).to_numpy()
    spatial.obs["SampleID"] = metadata.loc[
        spatial.obs_names, "biosample_id"
    ].astype(str).to_numpy()
    spatial = spatial[
        (spatial.obs["State"] == "normal")
        & (spatial.obs["sample"] == "2020-09-14_Puck_200701_21")
    ].copy()
    reference = reference[reference.obs["State"].astype(str) == "normal"].copy()
    spatial.obs[["x", "y"]] = spatial.obs[["x", "y"]] / 0.65
    expressed = np.asarray(spatial.X.sum(axis=0)).ravel() != 0
    spatial = spatial[:, expressed].copy()
    spatial.obsm["spatial"] = spatial.obs[["x", "y"]].to_numpy(dtype=float)
    spatial.obs["transcript_counts"] = np.asarray(spatial.X.sum(axis=1)).ravel()
    reference.obs["Level1"] = reference.obs["labels"].astype(str)
    reference.obs["Level2"] = reference.obs["labels"].astype(str)
    reference.obs["transcript_counts"] = np.asarray(reference.X.sum(axis=1)).ravel()
    spatial.var_names_make_unique()
    reference.var_names_make_unique()

    if tuple(spatial.shape) != (33_673, 20_388):
        raise ValueError(f"Unexpected colon puck shape: {spatial.shape}")
    if tuple(reference.shape) != (17_512, 31_053):
        raise ValueError(f"Unexpected colon reference shape: {reference.shape}")
    if reference.obs["Level1"].nunique() != 9:
        raise ValueError("Expected nine colon reference labels")
    filtered_beads = int((spatial.obs["transcript_counts"] >= 50).sum())
    if filtered_beads != 22_089:
        raise ValueError(f"Expected 22,089 colon beads after counts>=50, got {filtered_beads}")
    raw_overlap = spatial.var_names.intersection(reference.var_names)
    if len(raw_overlap) != 18_214:
        raise ValueError(f"Expected 18,214 raw colon shared genes, got {len(raw_overlap)}")
    retained_beads = spatial.obs["transcript_counts"].to_numpy(dtype=float) >= 50
    spatial_expressed = (
        np.asarray((spatial[retained_beads, raw_overlap].X != 0).sum(axis=0)).ravel() >= 1
    )
    reference_expressed = (
        np.asarray((reference[:, raw_overlap].X != 0).sum(axis=0)).ravel() >= 1
    )
    effective_overlap = raw_overlap[spatial_expressed & reference_expressed]
    if len(effective_overlap) != 16_841:
        raise ValueError(
            "Expected 16,841 effective colon genes after the configured filters, "
            f"got {len(effective_overlap)}"
        )

    st_path = prepared / "normal_Puck_200701_21.h5ad"
    reference_path = prepared / "normal_scRNAseq.h5ad"
    _atomic_h5ad(st_path, spatial)
    _atomic_h5ad(reference_path, reference)
    manifest = _write_manifest(
        "slideseq_mouse_colon",
        {
            "source_files": source_records,
            "transformations": [
                "State=normal and puck=2020-09-14_Puck_200701_21",
                "spatial x/y divided by 0.65",
                "all-zero spatial genes removed",
                "Level1=Level2=labels",
                "REVISE config applies counts>=50 to spatial beads",
            ],
            "observed": {
                "raw_st_shape": list(spatial.shape),
                "reference_shape": list(reference.shape),
                "beads_after_counts_50": filtered_beads,
                "raw_shared_genes": len(raw_overlap),
                "shared_genes": len(effective_overlap),
                "reference_labels": int(reference.obs["Level1"].nunique()),
            },
            "prepared_files": {
                "st": _adata_identity(st_path, spatial),
                "reference": _adata_identity(reference_path, reference),
            },
        },
    )
    print(manifest)
    return manifest


def _human_liver_source_records(raw: Path) -> dict:
    missing = [
        specification["name"]
        for specification in HUMAN_LIVER_SOURCE_FILES.values()
        if not (raw / specification["name"]).is_file()
    ]
    if missing:
        raise RuntimeError(
            "HUMAN_LIVER_SOURCE_REQUIRED: download the two pinned Dryad v5 files "
            "from the official file pages and place them in "
            f"{str(raw)!r}; missing now: {missing!r}"
        )

    records = {}
    for role, specification in HUMAN_LIVER_SOURCE_FILES.items():
        path = raw / specification["name"]
        expected_bytes = int(specification["bytes"])
        observed_bytes = path.stat().st_size
        if observed_bytes != expected_bytes:
            raise RuntimeError(
                f"Human-liver {role} source size mismatch for {path}: "
                f"{observed_bytes:,} != {expected_bytes:,}"
            )
        observed_sha256 = sha256_file(path)
        if observed_sha256 != specification["sha256"]:
            raise RuntimeError(
                f"Human-liver {role} source SHA-256 mismatch for {path}: "
                f"{observed_sha256} != {specification['sha256']}"
            )
        records[role] = {
            **specification,
            "path": str(path),
            "bytes": observed_bytes,
            "sha256": observed_sha256,
        }
    return records


def _human_liver_backed_matrix(
    matrix,
    rows: np.ndarray,
    columns: np.ndarray | None,
):
    if columns is None:
        result = matrix[rows, :]
    elif isinstance(matrix, np.ndarray):
        result = matrix[np.ix_(rows, columns)]
    else:
        try:
            result = matrix[rows, columns]
        except (IndexError, TypeError, ValueError):
            # h5py datasets allow only one vector-valued index at a time. Read
            # contiguous source-gene runs so the AM042 subset is never widened
            # to the full snRNA-seq gene axis before materialization.
            order = np.argsort(columns)
            sorted_columns = columns[order]
            blocks = []
            start = 0
            while start < len(sorted_columns):
                stop = start
                while (
                    stop + 1 < len(sorted_columns)
                    and sorted_columns[stop + 1] == sorted_columns[stop] + 1
                ):
                    stop += 1
                left = int(sorted_columns[start])
                right = int(sorted_columns[stop]) + 1
                blocks.append(np.asarray(matrix[rows, left:right]))
                start = stop + 1
            sorted_result = np.concatenate(blocks, axis=1)
            result = sorted_result[:, np.argsort(order)]
    if sparse.issparse(result):
        return result.tocsr()
    return np.asarray(result)


def _human_liver_raw_subset(
    source: ad.AnnData,
    *,
    role: str,
    sample_id: str,
    panel: pd.Index | None = None,
) -> tuple[sparse.csr_matrix | np.ndarray, pd.DataFrame, pd.DataFrame]:
    if source.raw is None:
        raise KeyError(f"Human-liver {role} source lacks raw.X")
    if "sample_id" not in source.obs:
        raise KeyError(f"Human-liver {role} source lacks obs['sample_id']")
    if source.raw.n_obs != source.n_obs:
        raise ValueError(
            f"Human-liver {role} raw.X and obs row counts differ: "
            f"{source.raw.n_obs} != {source.n_obs}"
        )

    sample_values = source.obs["sample_id"].astype(str).to_numpy()
    rows = np.flatnonzero(sample_values == sample_id)
    if len(rows) == 0:
        raise ValueError(
            f"Human-liver {role} source has no sample_id={sample_id!r} cells"
        )

    raw_var_names = pd.Index(source.raw.var_names.astype(str), name="gene")
    if not raw_var_names.is_unique:
        raise ValueError(f"Human-liver {role} raw gene names are not unique")
    if panel is None:
        selected_var_names = raw_var_names
        columns = None
    else:
        selected_var_names = pd.Index(panel.astype(str), name="gene").intersection(
            raw_var_names, sort=False
        )
        if len(selected_var_names) == 0:
            raise ValueError(
                f"Human-liver {role} source has no genes in the requested panel"
            )
        columns = raw_var_names.get_indexer(selected_var_names)
    matrix = _human_liver_backed_matrix(source.raw.X, rows, columns)
    obs = source.obs.iloc[rows].copy()
    if columns is None:
        var = source.raw.var.copy()
        var.index = raw_var_names
    else:
        var = source.raw.var.iloc[columns].copy()
        var.index = selected_var_names
    return matrix, obs, var


def _human_liver_level1(label: str) -> str:
    text = str(label)
    if text.startswith("Hep"):
        return "Hepatocyte"
    if text.startswith("HSC"):
        return "HSC"
    if text.startswith("Mac"):
        return "Macrophage"
    return text


def prepare_human_liver() -> Path:
    root = data_root()
    raw = root / "raw" / "human_liver_merfish"
    prepared = root / "prepared" / "human_liver_merfish" / "AM042"
    raw.mkdir(parents=True, exist_ok=True)
    prepared.mkdir(parents=True, exist_ok=True)

    source_records = _human_liver_source_records(raw)
    merfish_source_path = raw / HUMAN_LIVER_SOURCE_FILES["merfish"]["name"]
    nucseq_source_path = raw / HUMAN_LIVER_SOURCE_FILES["nucseq"]["name"]

    merfish_source = ad.read_h5ad(merfish_source_path, backed="r")
    try:
        merfish_matrix, merfish_obs, merfish_var = _human_liver_raw_subset(
            merfish_source,
            role="MERFISH",
            sample_id="AM042",
        )
    finally:
        merfish_source.file.close()

    required_merfish_obs = {"x", "y", "Cell_Type"}
    missing_merfish_obs = sorted(required_merfish_obs - set(merfish_obs.columns))
    if missing_merfish_obs:
        raise KeyError(
            f"Human-liver MERFISH source is missing obs fields: {missing_merfish_obs}"
        )
    merfish_coordinates = merfish_obs[["x", "y"]].to_numpy(dtype=float)
    if not np.isfinite(merfish_coordinates).all():
        raise ValueError("Human-liver MERFISH x/y coordinates are non-finite")
    if merfish_obs["Cell_Type"].isna().any():
        raise ValueError("Human-liver MERFISH Cell_Type contains missing values")

    label_sidecar = merfish_obs[["Cell_Type"]].copy()
    label_sidecar.index = pd.Index(label_sidecar.index.astype(str), name="cell_id")
    merfish_obs = merfish_obs.drop(columns=["Cell_Type"])
    spatial = ad.AnnData(X=merfish_matrix, obs=merfish_obs, var=merfish_var)
    spatial.obsm["spatial"] = merfish_coordinates
    spatial.obs["transcript_counts"] = np.asarray(
        spatial.X.sum(axis=1)
    ).ravel()
    spatial.uns["human_liver_preparation"] = {
        "dataset": "Watson et al. healthy human liver MERFISH and snRNA-seq",
        "doi": HUMAN_LIVER_DOI,
        "dryad_version": HUMAN_LIVER_DRYAD_VERSION,
        "dryad_api_version": HUMAN_LIVER_DRYAD_API_VERSION,
        "role": "MERFISH",
        "sample_id": "AM042",
        "expression_source": "raw.X",
        "annotation_source": "sidecar/MERFISH_Cell_Type.csv.gz",
        "source_file": source_records["merfish"]["name"],
        "source_sha256": source_records["merfish"]["sha256"],
    }

    nucseq_source = ad.read_h5ad(nucseq_source_path, backed="r")
    try:
        merfish_panel = pd.Index(spatial.var_names.astype(str), name="gene")
        nucseq_matrix, nucseq_obs, nucseq_var = _human_liver_raw_subset(
            nucseq_source,
            role="snRNA-seq",
            sample_id="AM042",
            panel=merfish_panel,
        )
    finally:
        nucseq_source.file.close()

    required_nucseq_obs = {"cell_type_final", "seurat_clusters"}
    missing_nucseq_obs = sorted(required_nucseq_obs - set(nucseq_obs.columns))
    if missing_nucseq_obs:
        raise KeyError(
            "Human-liver snRNA-seq source is missing obs fields: "
            f"{missing_nucseq_obs}"
        )
    shared_genes = pd.Index(nucseq_var.index.astype(str), name="gene")

    reference = ad.AnnData(X=nucseq_matrix, obs=nucseq_obs, var=nucseq_var)
    reference.obs["Level1"] = reference.obs["cell_type_final"].map(
        _human_liver_level1
    )
    reference.obs["Level2"] = (
        reference.obs["Level1"].astype(str)
        + "__seurat_"
        + reference.obs["seurat_clusters"].astype(str)
    )
    reference.obs["transcript_counts"] = np.asarray(
        reference.X.sum(axis=1)
    ).ravel()
    reference.uns["human_liver_preparation"] = {
        "dataset": "Watson et al. healthy human liver MERFISH and snRNA-seq",
        "doi": HUMAN_LIVER_DOI,
        "dryad_version": HUMAN_LIVER_DRYAD_VERSION,
        "dryad_api_version": HUMAN_LIVER_DRYAD_API_VERSION,
        "role": "snRNA-seq reference",
        "sample_id": "AM042",
        "expression_source": "raw.X",
        "annotation_source": "obs/cell_type_final and obs/seurat_clusters",
        "gene_panel_source": "MERFISH raw.var_names",
        "source_file": source_records["nucseq"]["name"],
        "source_sha256": source_records["nucseq"]["sha256"],
    }

    st_path = prepared / "MERFISH.h5ad"
    reference_path = prepared / "snRNAseq.h5ad"
    sidecar_path = prepared / "MERFISH_Cell_Type.csv.gz"
    _atomic_h5ad(st_path, spatial)
    _atomic_h5ad(reference_path, reference)
    temporary_sidecar = sidecar_path.with_suffix(sidecar_path.suffix + ".tmp")
    label_sidecar.to_csv(temporary_sidecar, compression="gzip")
    os.replace(temporary_sidecar, sidecar_path)

    manifest = _write_manifest(
        "human_liver_merfish_AM042",
        {
            "dataset": "Watson et al. healthy human liver MERFISH and snRNA-seq",
            "doi": HUMAN_LIVER_DOI,
            "dryad": {
                "landing_url": HUMAN_LIVER_DRYAD_LANDING_URL,
                "archive_version": HUMAN_LIVER_DRYAD_VERSION,
                "api_version": HUMAN_LIVER_DRYAD_API_VERSION,
            },
            "source_files": source_records,
            "transformations": [
                "backed-read the pinned Dryad v5 H5AD files",
                "select sample_id=AM042 before materializing raw.X",
                "MERFISH expression=raw.X; spatial=[obs/x,obs/y]",
                "MERFISH Cell_Type is stored only in MERFISH_Cell_Type.csv.gz",
                "snRNA-seq expression=raw.X restricted to the MERFISH gene panel",
                "Level1 maps Hep*=Hepatocyte, HSC*=HSC, Mac*=Macrophage; other labels are retained",
                "Level2=Level1 + '__seurat_' + seurat_clusters",
            ],
            "observed": {
                "sample_id": "AM042",
                "st_shape": list(spatial.shape),
                "reference_shape": list(reference.shape),
                "shared_panel_genes": len(shared_genes),
                "reference_level1": reference.obs["Level1"]
                .value_counts()
                .astype(int)
                .to_dict(),
                "reference_level2": int(reference.obs["Level2"].nunique()),
                "sidecar_rows": len(label_sidecar),
            },
            "prepared_files": {
                "st": _adata_identity(st_path, spatial),
                "reference": _adata_identity(reference_path, reference),
                "label_sidecar": {
                    "path": str(sidecar_path),
                    "bytes": sidecar_path.stat().st_size,
                    "sha256": sha256_file(sidecar_path),
                    "rows": len(label_sidecar),
                    "column": "Cell_Type",
                },
            },
        },
    )
    print(manifest)
    return manifest


def _verified_local_source(root: Path, specification: dict, role: str) -> dict:
    path = root / specification["name"]
    if not path.is_file():
        raise RuntimeError(f"Missing {role} source file: {path}")
    observed_bytes = path.stat().st_size
    expected_bytes = int(specification["bytes"])
    if observed_bytes != expected_bytes:
        raise RuntimeError(
            f"{role} source size mismatch for {path}: "
            f"{observed_bytes:,} != {expected_bytes:,}"
        )
    observed_sha256 = sha256_file(path)
    if observed_sha256 != specification["sha256"]:
        raise RuntimeError(
            f"{role} source SHA-256 mismatch for {path}: "
            f"{observed_sha256} != {specification['sha256']}"
        )
    return {
        **specification,
        "path": str(path),
        "bytes": observed_bytes,
        "sha256": observed_sha256,
    }


def _extract_allen_visp_archive(archive_path: Path, zarr_path: Path) -> None:
    if zarr_path.is_dir():
        return
    import zipfile

    _require_headroom(zarr_path.parent, 3 * 1024**3, "Allen VISp Zarr3 extraction")
    extraction_root = zarr_path.parent / "allen_visp_extracting"
    if extraction_root.exists():
        raise RuntimeError(
            f"Incomplete Allen VISp extraction must be inspected first: {extraction_root}"
        )
    extraction_root.mkdir(parents=True)
    with zipfile.ZipFile(archive_path) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise RuntimeError(f"Allen VISp ZIP CRC failure: {bad_member}")
        roots = {Path(name).parts[0] for name in archive.namelist() if name}
        if roots != {"data.zarr"}:
            raise ValueError(f"Unexpected Allen VISp ZIP roots: {sorted(roots)!r}")
        archive.extractall(extraction_root)
    extracted = extraction_root / "data.zarr"
    if not (extracted / "zarr.json").is_file():
        raise RuntimeError("Allen VISp archive did not produce data.zarr/zarr.json")
    os.replace(extracted, zarr_path)
    extraction_root.rmdir()


def _convert_allen_visp_zarr3(zarr_path: Path, output_path: Path) -> dict:
    converter = Path(__file__).with_name("convert_allen_visp_zarr3.py")
    if output_path.is_file():
        converted = ad.read_h5ad(output_path, backed="r")
        try:
            if tuple(converted.shape) != (2_389, 268):
                raise ValueError(
                    f"Unexpected existing Zarr3 conversion shape: {converted.shape}"
                )
            conversion = _json_compatible(
                dict(converted.uns.get("zarr3_conversion", {}))
            )
        finally:
            converted.file.close()
        if not conversion:
            raise KeyError(f"Existing conversion lacks provenance: {output_path}")
        return {
            "path": str(output_path),
            "bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
            "conversion": conversion,
            "converter_path": str(converter),
            "converter_sha256": sha256_file(converter),
            "created_now": False,
        }

    temporary = output_path.with_suffix(".tmp.h5ad")
    if temporary.exists():
        raise RuntimeError(
            f"Incomplete Allen VISp conversion must be inspected first: {temporary}"
        )
    command = [
        "uv",
        "run",
        "--isolated",
        "--python",
        "3.11",
        "--with",
        "zarr>=3,<4",
        "--with",
        "anndata>=0.12,<0.13",
        "python",
        str(converter),
        str(zarr_path),
        str(temporary),
    ]
    converter_environment = os.environ.copy()
    converter_environment.pop("PYTHONPATH", None)
    for key in tuple(converter_environment):
        if key.startswith("CONDA_"):
            converter_environment.pop(key)
    subprocess.run(
        command,
        cwd=repo_root(),
        env=converter_environment,
        check=True,
    )
    converted = ad.read_h5ad(temporary, backed="r")
    try:
        if tuple(converted.shape) != (2_389, 268):
            raise ValueError(f"Unexpected Zarr3 conversion shape: {converted.shape}")
        conversion = _json_compatible(dict(converted.uns["zarr3_conversion"]))
    finally:
        converted.file.close()
    os.replace(temporary, output_path)
    return {
        "path": str(output_path),
        "bytes": output_path.stat().st_size,
        "sha256": sha256_file(output_path),
        "conversion": conversion,
        "converter_path": str(converter),
        "converter_sha256": sha256_file(converter),
        "created_now": True,
    }


def _backed_matrix_columns(matrix, columns: np.ndarray, n_obs: int):
    blocks = []
    for start in range(0, n_obs, 2_048):
        block = matrix[start : min(start + 2_048, n_obs), columns]
        blocks.append(
            block.tocsr() if sparse.issparse(block) else sparse.csr_matrix(block)
        )
    return sparse.vstack(blocks, format="csr")


def prepare_allen_visp() -> Path:
    import geopandas as gpd

    root = data_root()
    raw = root / "raw" / "allen_visp_merfish"
    prepared = root / "prepared" / "allen_visp_merfish"
    raw.mkdir(parents=True, exist_ok=True)
    prepared.mkdir(parents=True, exist_ok=True)

    archive_record = _verified_local_source(raw, ALLEN_VISP_ARCHIVE, "MERFISH ZIP")
    reference_record = _verified_local_source(
        raw, ALLEN_VISP_TASIC_REFERENCE, "Tasic cortex reference"
    )
    zarr_path = raw / "data.zarr"
    _extract_allen_visp_archive(Path(archive_record["path"]), zarr_path)
    converted_path = prepared / "source_table_zarr3.h5ad"
    conversion_record = _convert_allen_visp_zarr3(zarr_path, converted_path)

    source_table = ad.read_h5ad(converted_path)
    if tuple(source_table.shape) != (2_389, 268):
        raise ValueError(f"Unexpected converted MERFISH table shape: {source_table.shape}")
    required_obs = {"cell_id", "region"}
    if set(source_table.obs.columns) != required_obs:
        raise ValueError(
            f"Unexpected converted MERFISH obs fields: {list(source_table.obs.columns)!r}"
        )
    blank_mask = source_table.var_names.str.startswith("Blank-")
    blank_genes = source_table.var_names[blank_mask].astype(str).tolist()
    if set(blank_genes) != {f"Blank-{index}" for index in range(1, 11)}:
        raise ValueError(f"Unexpected MERFISH blank controls: {blank_genes}")
    panel = pd.Index(source_table.var_names[~blank_mask].astype(str), name="gene")
    if len(panel) != 258 or not panel.is_unique:
        raise ValueError(f"Unexpected non-blank MERFISH panel: {len(panel)} genes")

    cell_ids = source_table.obs["cell_id"].to_numpy(dtype=int)
    if not np.array_equal(np.sort(cell_ids), np.arange(2_389)):
        raise ValueError("MERFISH cell_id is not a permutation of 0..2388")
    shapes_path = zarr_path / "shapes" / "cells" / "shapes.parquet"
    cells = gpd.read_parquet(shapes_path)
    cells.index = pd.Index(cells.index.to_numpy(dtype=int), name="cell_id")
    if len(cells) != 2_389 or not cells.index.is_unique:
        raise ValueError(f"Unexpected Allen VISp cell shapes: {cells.shape}")
    if set(cells.geometry.geom_type.astype(str)) != {"Point"}:
        raise ValueError("Allen VISp linked cell shapes must all be Points")
    cells = cells.loc[cell_ids]
    coordinates = np.column_stack(
        [cells.geometry.x.to_numpy(dtype=float), cells.geometry.y.to_numpy(dtype=float)]
    )
    if not np.isfinite(coordinates).all():
        raise ValueError("Allen VISp linked cell coordinates are non-finite")

    spatial_obs = source_table.obs.drop(columns=["cell_id"]).copy()
    spatial_obs["source_obs_name"] = source_table.obs_names.astype(str)
    spatial_obs.index = pd.Index(cell_ids.astype(str), name="cell_id")
    spatial_obs["x"] = coordinates[:, 0]
    spatial_obs["y"] = coordinates[:, 1]
    spatial_obs["radius"] = cells["radius"].to_numpy(dtype=float)
    spatial_matrix = sparse.csr_matrix(
        np.asarray(source_table[:, panel].X, dtype=np.int32)
    )
    spatial = ad.AnnData(
        X=spatial_matrix,
        obs=spatial_obs,
        var=source_table.var.loc[panel].copy(),
    )
    spatial.obsm["spatial"] = coordinates
    spatial.obs["transcript_counts"] = np.asarray(
        spatial.X.sum(axis=1)
    ).ravel()
    spatial.uns["allen_visp_preparation"] = {
        "dataset": "Allen Institute prototype MERFISH, mouse VISp",
        "role": "segmented-cell MERFISH target",
        "expression_source": "SpatialData tables/table/X",
        "coordinate_source": "SpatialData shapes/cells Point geometry",
        "coordinate_units": "not specified by the archive",
        "removed_controls": blank_genes,
        "archive_sha256": archive_record["sha256"],
    }

    reference_source = ad.read_h5ad(reference_record["path"], backed="r")
    try:
        if tuple(reference_source.shape) != (21_697, 36_826):
            raise ValueError(
                f"Unexpected Squidpy mouse-cortex shape: {reference_source.shape}"
            )
        if reference_source.raw is None:
            raise KeyError("Squidpy mouse-cortex reference lacks raw.X counts")
        required_reference_obs = {"cell_class", "cell_subclass", "cell_cluster"}
        missing_reference_obs = sorted(
            required_reference_obs - set(reference_source.obs.columns)
        )
        if missing_reference_obs:
            raise KeyError(
                f"Squidpy mouse-cortex reference lacks labels: {missing_reference_obs}"
            )
        raw_var_names = pd.Index(reference_source.raw.var_names.astype(str), name="gene")
        columns = raw_var_names.get_indexer(panel)
        if (columns < 0).any():
            raise ValueError(
                "Tasic mouse-cortex reference is missing MERFISH panel genes: "
                f"{panel[columns < 0].tolist()}"
            )
        reference_matrix = _backed_matrix_columns(
            reference_source.raw.X, columns, reference_source.n_obs
        )
        reference_obs = reference_source.obs.copy()
        reference_var = reference_source.raw.var.iloc[columns].copy()
        reference_var.index = panel
    finally:
        reference_source.file.close()
    if not np.isfinite(reference_matrix.data).all() or (
        reference_matrix.data < 0
    ).any():
        raise ValueError("Tasic raw.X panel contains invalid counts")
    if not np.allclose(reference_matrix.data, np.rint(reference_matrix.data)):
        raise ValueError("Tasic raw.X panel is not integer-valued count data")

    reference = ad.AnnData(
        X=reference_matrix.astype(np.float32),
        obs=reference_obs,
        var=reference_var,
    )
    reference.obs["Level1"] = reference.obs["cell_class"].astype(str)
    reference.obs["Level2"] = reference.obs["cell_cluster"].astype(str)
    reference.obs["transcript_counts"] = np.asarray(
        reference.X.sum(axis=1)
    ).ravel()
    reference.uns["allen_visp_preparation"] = {
        "dataset": "Tasic et al. 2018 mouse cortex scRNA-seq",
        "doi": ALLEN_VISP_TASIC_REFERENCE["doi"],
        "role": "external same-species cortex reference",
        "expression_source": "raw.X restricted to the MERFISH non-blank panel",
        "Level1": "cell_class",
        "Level2": "cell_cluster",
        "pairing_boundary": (
            "regional/species match only; not the same cells, animal, or experiment"
        ),
        "source_sha256": reference_record["sha256"],
    }

    if not spatial.var_names.equals(reference.var_names):
        raise ValueError("Prepared MERFISH and Tasic gene axes are not identical")
    spatial_expressed = np.asarray((spatial.X != 0).sum(axis=0)).ravel() >= 1
    reference_expressed = np.asarray((reference.X != 0).sum(axis=0)).ravel() >= 1
    effective_genes = spatial.var_names[spatial_expressed & reference_expressed]
    if len(effective_genes) < 2:
        raise ValueError("Allen VISp and Tasic reference have fewer than two usable genes")
    selected_types = ("Glutamatergic", "GABAergic", "Non-Neuronal")
    broad_counts = reference.obs["Level1"].value_counts().astype(int)
    subtype_counts = (
        reference.obs.groupby("Level1", observed=True)["Level2"]
        .nunique()
        .astype(int)
    )
    for label in selected_types:
        if broad_counts.get(label, 0) < 2 or subtype_counts.get(label, 0) < 2:
            raise ValueError(
                f"Tasic reference lacks a usable {label} subset: "
                f"cells={broad_counts.get(label, 0)}, "
                f"subtypes={subtype_counts.get(label, 0)}"
            )

    st_path = prepared / "MERFISH.h5ad"
    reference_path = prepared / "Tasic2018_sc_mouse_cortex.h5ad"
    _atomic_h5ad(st_path, spatial)
    _atomic_h5ad(reference_path, reference)
    manifest = _write_manifest(
        "allen_visp_merfish_tasic2018",
        {
            "dataset": "Allen Institute prototype MERFISH mouse VISp",
            "source_identity": ALLEN_VISP_SOURCE_EVIDENCE,
            "identity_boundary": (
                "The archive generator and SpaceJam source identify Allen VISp; "
                "the archive does not substantiate a Moffitt 2018 hypothalamic identity."
            ),
            "reference_pairing": {
                "dataset": "Tasic et al. 2018 mouse cortex scRNA-seq",
                "doi": ALLEN_VISP_TASIC_REFERENCE["doi"],
                "relationship": "same species and cortex-region match",
                "not_claimed": ["same cell", "same animal", "same experiment"],
            },
            "source_files": {
                "merfish_zip": archive_record,
                "tasic_reference": reference_record,
                "zarr3_conversion": conversion_record,
            },
            "transformations": [
                "convert only tables/table from Zarr3 to an intermediate H5AD in an isolated Python 3.11 environment",
                "replace source obs_names with table cell_id and retain source_obs_name as metadata",
                "join shapes/cells Point centroids to table rows by cell_id",
                "remove Blank-1 through Blank-10 controls",
                "recompute MERFISH transcript_counts from non-blank counts",
                "use Tasic raw.X counts restricted to the 258-gene MERFISH panel",
                "Level1=cell_class; Level2=cell_cluster",
            ],
            "observed": {
                "source_table_shape": [2_389, 268],
                "st_shape": list(spatial.shape),
                "reference_shape": list(reference.shape),
                "nonblank_panel_genes": len(panel),
                "effective_genes_after_min_cell_1": len(effective_genes),
                "st_cells_after_counts_15": int(
                    (spatial.obs["transcript_counts"] >= 15).sum()
                ),
                "reference_level1_cells": broad_counts.to_dict(),
                "reference_level1_subtypes": subtype_counts.to_dict(),
                "selected_types": list(selected_types),
                "coordinate_bounds": {
                    "x_min": float(coordinates[:, 0].min()),
                    "x_max": float(coordinates[:, 0].max()),
                    "y_min": float(coordinates[:, 1].min()),
                    "y_max": float(coordinates[:, 1].max()),
                },
            },
            "prepared_files": {
                "st": _adata_identity(st_path, spatial),
                "reference": _adata_identity(reference_path, reference),
            },
        },
    )
    print(manifest)
    return manifest


def _decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O"}:
        return np.asarray(
            [value.decode() if isinstance(value, bytes) else value for value in values]
        )
    return values


def _read_osmfish_loom(path: Path) -> ad.AnnData:
    with h5py.File(path, "r") as handle:
        matrix = np.asarray(handle["matrix"]).T
        row_attrs = {
            key: _decode(np.asarray(dataset))
            for key, dataset in handle["row_attrs"].items()
        }
        col_attrs = {
            key: _decode(np.asarray(dataset))
            for key, dataset in handle["col_attrs"].items()
        }
    if "Gene" not in row_attrs or "CellID" not in col_attrs:
        raise KeyError("osmFISH loom lacks Gene or CellID attributes")
    obs = pd.DataFrame(col_attrs)
    obs.index = pd.Index(obs["CellID"].astype(str), name="CellID")
    var = pd.DataFrame(row_attrs)
    var.index = pd.Index(var["Gene"].astype(str), name="Gene")
    return ad.AnnData(X=sparse.csr_matrix(matrix), obs=obs, var=var)


def _hybridization_annotation(label: str) -> str:
    if label.endswith("_Hybridization4"):
        return "Partly_out_of_focus"
    if label.startswith(("Klk6_", "Lum_")):
        return "Low quality"
    if label == "Tbr1_Hybridization11":
        return "Internal control"
    if label.startswith(("Cnr1", "Plp1_", "Vtn_")):
        return "Repeat Round 4"
    return ""


def _molecule_to_original_cells(
    coordinates: np.ndarray,
    segmentation_path: Path,
) -> np.ndarray:
    import pickle

    codes = (coordinates[:, 0].astype(np.uint64) << np.uint64(32)) | coordinates[
        :, 1
    ].astype(np.uint64)
    unique_codes, inverse = np.unique(codes, return_inverse=True)
    assignment = np.full(len(unique_codes), "", dtype=object)
    with segmentation_path.open("rb") as handle:
        segments = pickle.load(handle)
    for cell_id, pixels in segments.items():
        array = np.asarray(pixels)
        if array.size == 0:
            continue
        pixel_codes = (
            array[:, 0].astype(np.uint64) << np.uint64(32)
        ) | array[:, 1].astype(np.uint64)
        positions = np.searchsorted(unique_codes, pixel_codes)
        valid = positions < len(unique_codes)
        positions = positions[valid]
        pixel_codes = pixel_codes[valid]
        matched = unique_codes[positions] == pixel_codes
        positions = positions[matched]
        positions = positions[assignment[positions] == ""]
        assignment[positions] = str(cell_id)
    del segments
    gc.collect()
    return assignment[inverse]


def _dominant_label(frame: pd.DataFrame, label: str) -> pd.Series:
    counts = (
        frame.loc[frame[label].notna(), ["segment", label]]
        .groupby(["segment", label], observed=True)
        .size()
        .rename("n")
        .reset_index()
        .sort_values(["segment", "n", label], ascending=[True, False, True])
        .drop_duplicates("segment")
        .set_index("segment")[label]
    )
    return counts.astype(str)


def prepare_osmfish() -> Path:
    root = data_root()
    raw = root / "raw" / "osmfish"
    prepared = root / "prepared" / "osmfish"
    raw.mkdir(parents=True, exist_ok=True)
    prepared.mkdir(parents=True, exist_ok=True)
    raw_paths = {
        "loom": raw / "osmFISH_SScortex_mouse_all_cells.loom",
        "molecules": raw / "mRNA_coords_raw_counting.hdf5",
        "segmentation": raw / "polyT_seg.pkl",
    }
    source_records = {}
    for key, path in raw_paths.items():
        record = _download(
            OSMFISH_EFFECTIVE_URLS[key],
            path,
            expected_bytes=OSMFISH_EXPECTED_BYTES[key],
        )
        record["original_url"] = OSMFISH_ORIGINAL_URLS[key]
        record["original_url_probe"] = _probe_url(OSMFISH_ORIGINAL_URLS[key])
        record["mirror_url"] = OSMFISH_EFFECTIVE_URLS[key]
        source_records[key] = record

    reference_all = _read_osmfish_loom(raw_paths["loom"])
    if tuple(reference_all.shape) != (6_471, 33):
        raise ValueError(f"Unexpected osmFISH loom shape: {reference_all.shape}")
    um_per_pixel = float(
        np.sqrt(
            reference_all.obs["size_um2"].astype(float)
            / reference_all.obs["size_pix"].astype(float)
        ).mean()
    )
    reference = reference_all[
        reference_all.obs["ClusterName"].astype(str) != "Excluded"
    ].copy()
    if tuple(reference.shape) != (4_839, 33):
        raise ValueError(f"Unexpected filtered osmFISH reference shape: {reference.shape}")
    reference.obs["Level2"] = reference.obs["ClusterName"].astype(str)
    reference.obs["Level1"] = reference.obs["Level2"].map(_broad_osmfish)
    reference.obs["transcript_counts"] = np.asarray(reference.X.sum(axis=1)).ravel()
    reference.obs["ClusterName"] = pd.Categorical(reference.obs["ClusterName"].astype(str))
    reference_path = prepared / "osmFISH_reference.h5ad"
    if not reference_path.exists():
        _atomic_h5ad(reference_path, reference)

    st_path = prepared / "TonT_segmented_cells.h5ad"
    sidecar_path = prepared / "TonT_labels.csv.gz"
    if not st_path.exists() or not sidecar_path.exists():
        genes_by_channel = {}
        coordinates = []
        genes = []
        with h5py.File(raw_paths["molecules"], "r") as handle:
            channels = list(handle.keys())
            good_channels = [
                channel
                for channel in channels
                if _hybridization_annotation(channel)
                not in {"Partly_out_of_focus", "Low quality", "Internal control"}
            ]
            channel_genes = {channel: channel.split("_")[0] for channel in good_channels}
            sm2seg = {
                source: target
                for source, target in zip(
                    sorted(channel_genes.values()), sorted(reference.var_names.astype(str))
                )
            }
            for channel in good_channels:
                gene = sm2seg[channel_genes[channel]]
                values = np.asarray(handle[channel], dtype=np.uint32)
                coordinates.append(values)
                genes.append(np.full(len(values), gene, dtype=object))
                genes_by_channel[channel] = gene
        coordinates_array = np.concatenate(coordinates, axis=0)
        gene_array = np.concatenate(genes)
        if len(coordinates_array) != 1_802_589:
            raise ValueError(
                f"Expected 1,802,589 filtered osmFISH molecules, got {len(coordinates_array)}"
            )

        original_cells = _molecule_to_original_cells(
            coordinates_array, raw_paths["segmentation"]
        )
        cluster_lookup = reference_all.obs["ClusterName"].astype(str).to_dict()
        original_cluster = pd.Series(original_cells).map(cluster_lookup).fillna("").to_numpy()
        valid_reference = (original_cluster != "") & (original_cluster != "Excluded")
        max_valid_x = int(coordinates_array[valid_reference, 0].max())
        valid_region = coordinates_array[:, 0] <= max_valid_x
        coordinates_array = coordinates_array[valid_region]
        gene_array = gene_array[valid_region]
        original_cluster = original_cluster[valid_region]
        del original_cells
        gc.collect()

        molecules = pd.DataFrame(
            {
                "gene": pd.Categorical(
                    gene_array,
                    categories=reference.var_names.astype(str),
                    ordered=True,
                ),
                "x": coordinates_array[:, 0].astype(np.float32) * um_per_pixel,
                "y": coordinates_array[:, 1].astype(np.float32) * um_per_pixel,
                "original_cluster": pd.Categorical(
                    pd.Series(original_cluster).replace("Excluded", "").replace("", np.nan)
                ),
            }
        )
        del coordinates_array, gene_array, original_cluster
        gc.collect()

        tc.tl.annotate_single_molecules(
            molecules,
            reference=reference,
            method="OT",
            annotation_key="ClusterName",
            result_key="tacco",
            bin_size=10,
            n_shifts=3,
            bisections=4,
            bisection_divisor=3,
            platform_iterations=-1,
        )
        tc.tl.segment(
            molecules,
            distance_scale=2.0,
            max_size=1600,
            result_key="tacco_seg",
            position_scale=10.0,
            position_range=2,
            annotation_key="tacco",
            annotation_distance=None,
        )
        segmented = molecules.loc[molecules["tacco_seg"].notna()].copy()
        segmented["segment"] = segmented["tacco_seg"].astype(str)
        row_codes, segment_names = pd.factorize(segmented["segment"], sort=False)
        gene_codes = segmented["gene"].cat.codes.to_numpy()
        if (gene_codes < 0).any():
            raise ValueError("Segmented osmFISH molecules contain genes outside reference")
        counts = sparse.coo_matrix(
            (
                np.ones(len(segmented), dtype=np.int16),
                (row_codes, gene_codes),
            ),
            shape=(len(segment_names), reference.n_vars),
        ).tocsr()
        n_per_segment = np.bincount(row_codes).astype(float)
        mean_x = np.bincount(row_codes, weights=segmented["x"].to_numpy()) / n_per_segment
        mean_y = np.bincount(row_codes, weights=segmented["y"].to_numpy()) / n_per_segment
        obs = pd.DataFrame(index=pd.Index(segment_names.astype(str), name="tacco_seg"))
        obs["x"] = mean_x
        obs["y"] = mean_y
        obs["transcript_counts"] = np.asarray(counts.sum(axis=1)).ravel()
        st = ad.AnnData(X=counts, obs=obs, var=reference.var.copy())
        st.obsm["spatial"] = obs[["x", "y"]].to_numpy(dtype=float)

        sidecar = pd.DataFrame(index=obs.index)
        sidecar["tacco_label"] = _dominant_label(segmented, "tacco").reindex(obs.index)
        sidecar["original_cluster"] = _dominant_label(
            segmented, "original_cluster"
        ).reindex(obs.index)
        sidecar["tacco_broad"] = sidecar["tacco_label"].map(_broad_osmfish)
        sidecar["transcript_counts"] = obs["transcript_counts"]
        _atomic_h5ad(st_path, st)
        temporary_sidecar = sidecar_path.with_suffix(sidecar_path.suffix + ".tmp")
        sidecar.to_csv(temporary_sidecar, compression="gzip")
        os.replace(temporary_sidecar, sidecar_path)
        del molecules, segmented, st, sidecar, counts
        gc.collect()

    st = ad.read_h5ad(st_path)
    sidecar = pd.read_csv(sidecar_path, index_col=0)
    sidecar.index = sidecar.index.astype(str)
    valid_counts = sidecar["transcript_counts"] >= 20
    class_counts = (
        sidecar.loc[valid_counts, "tacco_broad"].value_counts().astype(int).to_dict()
    )
    expected_reference = {
        "Pyramidal": {"cells": 1_944, "subtypes": 8},
        "Inhibitory": {"cells": 779, "subtypes": 7},
    }
    reference_class_counts = reference.obs["Level1"].value_counts().astype(int).to_dict()
    reference_subtype_counts = (
        reference.obs.groupby("Level1", observed=True)["Level2"]
        .nunique()
        .astype(int)
        .to_dict()
    )
    for label, expected in expected_reference.items():
        if (
            reference_class_counts.get(label) != expected["cells"]
            or reference_subtype_counts.get(label) != expected["subtypes"]
        ):
            raise ValueError(
                f"Unexpected osmFISH {label} reference subset: "
                f"cells={reference_class_counts.get(label)}, "
                f"subtypes={reference_subtype_counts.get(label)}"
            )
    manifest = _write_manifest(
        "osmfish",
        {
            "source_files": source_records,
            "transformations": [
                "remove Partly_out_of_focus, Low quality, and Internal control channels",
                "map genes by the TACCO sorted-name rule",
                "molecule coordinates converted by loom-derived um_per_pixel",
                "TACCO base annotate_single_molecules parameters only",
                "TonT=tacco_seg on TACCO molecule labels; aggregate by tacco_seg",
                "Level2=ClusterName; Level1 from explicit prefix mapping",
                "TonT labels are stored in a sidecar and excluded from REVISE ST input",
            ],
            "observed": {
                "raw_reference_shape": list(reference_all.shape),
                "filtered_reference_shape": list(reference.shape),
                "reference_cluster_names": int(reference.obs["Level2"].nunique()),
                "um_per_pixel": um_per_pixel,
                "tont_shape": list(st.shape),
                "tont_counts_ge_20_by_broad_label": class_counts,
                "reference_cells_by_broad_label": reference_class_counts,
                "reference_subtypes_by_broad_label": reference_subtype_counts,
            },
            "prepared_files": {
                "st": _adata_identity(st_path, st),
                "reference": _adata_identity(reference_path, reference),
                "label_sidecar": {
                    "path": str(sidecar_path),
                    "bytes": sidecar_path.stat().st_size,
                    "sha256": sha256_file(sidecar_path),
                    "rows": len(sidecar),
                },
            },
        },
    )
    print(manifest)
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Prepare TACCO-derived and matched-reference datasets for REVISE"
    )
    parser.add_argument(
        "case",
        choices=(
            "olfactory",
            "colon",
            "osmfish",
            "human_liver",
            "allen_visp",
        ),
    )
    args = parser.parse_args(argv)
    if args.case == "olfactory":
        prepare_olfactory()
    elif args.case == "colon":
        prepare_colon()
    elif args.case == "human_liver":
        prepare_human_liver()
    elif args.case == "allen_visp":
        prepare_allen_visp()
    elif args.case == "osmfish":
        prepare_osmfish()


if __name__ == "__main__":
    main()
