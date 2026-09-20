"""Prepare immutable, spatially contiguous real inputs and reconstruction requests.

Run from the REVISE checkout. This only selects input observations and writes
requests; reconstruction remains owned by reconstruct.py. Existing mini inputs
are refused. Full requests use the original sources, never the mini inputs.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SAMPLES = {
    "P2CRC_Xenium": ("Xenium.yaml", 0.2125, "confirmed by user; original counts"),
    "P1CRC_HD": ("VisiumHD.yaml", 0.27380817798463214, "source expression history unconfirmed"),
    "P2CRC_Visium": ("Visium.yaml", 0.73, "counts confirmed by user; coordinate scale provisional"),
}
COUNTS = {"identity": "raw_counts", "scale": "untransformed_nonnegative"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def square_roi(coordinates: np.ndarray, fraction: float) -> tuple[np.ndarray, dict]:
    """Keep every original observation in a median-centered square, including ties."""
    xy = np.asarray(coordinates, dtype=float)
    if xy.ndim != 2 or xy.shape[1] < 2 or len(xy) == 0 or not np.isfinite(xy).all():
        raise ValueError("ROI requires nonempty finite spatial coordinates")
    if not np.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    xy = xy[:, :2]
    center = np.median(xy, axis=0)
    distances = np.max(np.abs(xy - center), axis=1)
    target = int(np.ceil(len(xy) * fraction))
    radius = float(np.partition(distances, target - 1)[target - 1])
    selected = np.flatnonzero(distances <= radius)
    return selected, {
        "rule": "coordinate-wise median; Chebyshev order statistic; include boundary ties",
        "fraction": fraction,
        "source_n_obs": len(xy),
        "target_n_obs": target,
        "selected_n_obs": len(selected),
        "center": center.tolist(),
        "half_side": radius,
        "lower": (center - radius).tolist(),
        "upper": (center + radius).tolist(),
    }


def hst_calibration(source) -> dict:
    """Check stored pixel calibration against the observed Visium HD bin grid."""
    spatial = source.uns["spatial"]["Visium_HD_Human_Colon_Cancer_P1"]
    factors = spatial["scalefactors"]
    scale = float(factors["microns_per_pixel"])
    expected = float(factors["bin_size_um"])
    index = np.linspace(0, source.n_obs - 1, min(5000, source.n_obs), dtype=int)
    xy = np.asarray(source.obsm["spatial"])[index, :2]
    grid = source.obs.iloc[index][["array_row", "array_col"]].to_numpy()
    design = np.column_stack([grid, np.ones(len(grid))])
    fit = np.linalg.lstsq(design, xy, rcond=None)[0]
    steps = np.linalg.norm(fit[:2], axis=1) * scale
    if not np.allclose(steps, expected, rtol=0.01, atol=0):
        raise ValueError("Stored HD calibration does not match the observed bin grid")
    return {"source": "uns.spatial.scalefactors.microns_per_pixel",
            "microns_per_coordinate": scale, "bin_size_um": expected,
            "observed_grid_steps_um": steps.tolist(),
            "affine_max_residual_pixels": float(np.max(np.abs(design @ fit - xy)))}


def request(root: Path, sample: str, source_path: str, output: str, sample_id: str, mapping: str) -> dict:
    template, scale, _ = SAMPLES[sample]
    doc = yaml.safe_load((root / "configs/application" / template).read_text())
    doc["paths"]["root_dir"] = "."
    doc["inputs"]["st"]["path"] = source_path
    doc["inputs"]["reference"]["path"] = "raw_data/Real_application/adata_sc_all_reanno.h5ad"
    doc["inputs"]["reference"]["expression"] = dict(COUNTS)
    if sample.startswith("P2CRC"):
        doc["inputs"]["st"]["expression"] = dict(COUNTS)
        doc["inputs"]["reference"].update(filter_column="Patient", filter_value="P2CRC")
    else:
        doc["inputs"]["st"]["expression"] = {"identity": "unknown", "scale": "unknown"}
    doc["inputs"].pop("pm_on_cell", None)
    doc["output"] = {"dir": output}
    if sample == "P2CRC_Xenium":
        doc["output"]["ist_mapping"] = mapping
    doc["delivery"] = {"sample_id": sample_id,
                       "coordinates": {"key": "spatial", "unit": "pixel", "microns_per_coordinate": scale}}
    return doc


def prepare(root: Path, output: Path, configs: Path, fraction: float) -> dict:
    for sample in SAMPLES:
        mappings = ("random", "mean", "within_cluster", "outside_cluster") if sample == "P2CRC_Xenium" else ("default",)
        for size in ("mini", "full"):
            for mapping in mappings:
                path = configs / size / f"{sample}-{mapping}.yaml"
                if path.exists():
                    raise FileExistsError(f"Existing request is immutable; select a new --configs directory: {path}")
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    inputs.mkdir()
    evidence = {"kind": "real-data ROI preparation; not reconstruction or scientific acceptance",
                "fraction": fraction, "samples": {}}
    for sample, (_, scale, note) in SAMPLES.items():
        original = root / "raw_data/Real_application" / f"{sample}.h5ad"
        before = sha256(original)
        source = ad.read_h5ad(original, backed="r")
        try:
            positions, selection = square_roi(source.obsm["spatial"], fraction)
            calibration = hst_calibration(source) if sample == "P1CRC_HD" else {
                "source": "user confirmation 2026-09-20", "microns_per_coordinate": scale,
                "provisional": sample == "P2CRC_Visium"}
            mini = source[positions, :].to_memory()
            mini.uns["revise_acceptance_input"] = {
                "source_path": str(original.relative_to(root)), "source_sha256": before,
                "original_sample_id": sample, "mini_sample_id": f"{sample}_mini",
                "selection": deepcopy(selection), "coordinate_calibration": calibration,
                "source_note": note,
            }
            path = inputs / f"{sample}_mini.h5ad"
            mini.write_h5ad(path)
            selection.update(obs_ids=mini.obs_names.astype(str).tolist(), n_vars=mini.n_vars)
        finally:
            source.file.close()
        if sha256(original) != before:
            raise RuntimeError(f"Source changed during preparation: {original}")
        evidence["samples"][sample] = {
            "source": str(original.relative_to(root)), "source_sha256": before,
            "mini_input": str(path.relative_to(root)), "mini_sha256": sha256(path),
            "selection": selection, "coordinate_calibration": calibration, "source_note": note,
        }
        mappings = ("random", "mean", "within_cluster", "outside_cluster") if sample == "P2CRC_Xenium" else ("default",)
        for size in ("mini", "full"):
            for mapping in mappings:
                sample_id = f"{sample}_mini" if size == "mini" else sample
                st_path = path.relative_to(root) if size == "mini" else original.relative_to(root)
                destination = (output / "delivery" / sample_id / mapping).relative_to(root) if size == "mini" else Path("results/acceptance-full") / sample_id / mapping
                doc = request(root, sample, str(st_path), str(destination), sample_id, mapping)
                config = configs / size / f"{sample}-{mapping}.yaml"
                config.parent.mkdir(parents=True, exist_ok=True)
                header = f"# {size} real-data acceptance; launch from the REVISE root.\n# {note}\n"
                config.write_text(header + yaml.safe_dump(doc, sort_keys=False))
        del mini
    (output / "input_manifest.json").write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n")
    return evidence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=Path("output/mini-acceptance/20260920"))
    parser.add_argument("--configs", type=Path, default=Path("configs/acceptance"))
    parser.add_argument("--fraction", type=float, default=0.01)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (root / args.output).resolve()
    configs = (root / args.configs).resolve()
    evidence = prepare(root, output, configs, args.fraction)
    print(json.dumps({k: {"selected": v["selection"]["selected_n_obs"], "path": v["mini_input"]}
                      for k, v in evidence["samples"].items()}, indent=2))


if __name__ == "__main__":
    main()
