"""Two-stage Sim2Real-ST pseudo-spot workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml
from anndata import AnnData, read_h5ad

from revise.application.preprocess import (
    filter_reference,
    preprocess_reference,
    preprocess_spatial,
)
from revise.backend.kernels.ot import OTKernel
from revise.preprocess.sim2real_pseudospot.aggregation import aggregate_real_cells
from revise.preprocess.sim2real_pseudospot.regions import (
    Candidate,
    candidate_table,
    filter_region_cells,
    plot_candidates,
    plot_selected_region,
    propose_candidates,
)


REGION_PARTS = {
    "leading_edge": "part1",
    "normal_core": "part2",
    "tumor_core": "part3",
}
REQUIRED_SPOT_SIZES = (50, 100, 150, 200)


@dataclass(frozen=True)
class SampleConfig:
    xenium_path: Path
    output_dir: Path


@dataclass(frozen=True)
class PreprocessingConfig:
    transcript_counts_min: int
    gene_min_cells: int
    label_key: str
    unknown_label: str
    min_cells_per_type: int


@dataclass(frozen=True)
class ProposalConfig:
    base_width: float
    base_height: float
    scales: tuple[float, ...]
    step: float
    min_cells: int
    max_iou: float


@dataclass(frozen=True)
class Sim2RealPseudospotConfig:
    path: Path
    reference_path: Path
    patient_column: str
    samples: dict[str, SampleConfig]
    preprocessing: PreprocessingConfig
    template_root: Path
    proposal: ProposalConfig
    spot_sizes: tuple[int, ...]
    seed: int


@dataclass(frozen=True)
class ProposalResult:
    sample: str
    region_selection_dir: Path
    annotated_path: Path
    image_path: Path
    table_path: Path
    proposal_path: Path
    candidates: dict[str, list[Candidate]]


@dataclass(frozen=True)
class BuildResult:
    sample: str
    output_dir: Path
    parts: dict[str, Path]


def _required(mapping: dict[str, Any], key: str) -> Any:
    try:
        return mapping[key]
    except KeyError as exc:
        raise ValueError(f"Missing required configuration key: {key}") from exc


def _resolve(base: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else (base / path).resolve()


def load_config(path: str | Path) -> Sim2RealPseudospotConfig:
    """Load the versioned, path-stable Sim2Real-ST workflow configuration."""
    config_path = Path(path).resolve()
    with config_path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError("Sim2Real pseudo-spot configuration must be a mapping.")
    if raw.get("schema_version") != 1:
        raise ValueError("schema_version must be 1.")

    base = config_path.parent
    reference = _required(raw, "reference")
    preprocessing = _required(raw, "preprocessing")
    proposal = _required(raw, "proposal")
    raw_samples = _required(raw, "samples")
    if not isinstance(reference, dict) or not isinstance(raw_samples, dict):
        raise ValueError("reference and samples must be mappings.")
    if not isinstance(preprocessing, dict) or not isinstance(proposal, dict):
        raise ValueError("preprocessing and proposal must be mappings.")

    samples: dict[str, SampleConfig] = {}
    for sample, values in raw_samples.items():
        if not isinstance(values, dict):
            raise ValueError(f"samples.{sample} must be a mapping.")
        samples[str(sample)] = SampleConfig(
            xenium_path=_resolve(base, _required(values, "xenium_path")),
            output_dir=_resolve(base, _required(values, "output_dir")),
        )
    if not samples:
        raise ValueError("samples must not be empty.")

    spot_sizes = tuple(int(value) for value in _required(raw, "spot_sizes"))
    if spot_sizes != REQUIRED_SPOT_SIZES:
        raise ValueError(f"spot_sizes must be {list(REQUIRED_SPOT_SIZES)}.")

    return Sim2RealPseudospotConfig(
        path=config_path,
        reference_path=_resolve(base, _required(reference, "path")),
        patient_column=str(_required(reference, "patient_column")),
        samples=samples,
        preprocessing=PreprocessingConfig(
            transcript_counts_min=int(_required(preprocessing, "transcript_counts_min")),
            gene_min_cells=int(_required(preprocessing, "gene_min_cells")),
            label_key=str(_required(preprocessing, "label_key")),
            unknown_label=str(_required(preprocessing, "unknown_label")),
            min_cells_per_type=int(_required(preprocessing, "min_cells_per_type")),
        ),
        template_root=_resolve(base, _required(raw, "template_root")),
        proposal=ProposalConfig(
            base_width=float(_required(proposal, "base_width")),
            base_height=float(_required(proposal, "base_height")),
            scales=tuple(float(value) for value in _required(proposal, "scales")),
            step=float(_required(proposal, "step")),
            min_cells=int(_required(proposal, "min_cells")),
            max_iou=float(_required(proposal, "max_iou")),
        ),
        spot_sizes=spot_sizes,
        seed=int(_required(raw, "seed")),
    )


def _as_config(config: str | Path | Sim2RealPseudospotConfig) -> Sim2RealPseudospotConfig:
    return config if isinstance(config, Sim2RealPseudospotConfig) else load_config(config)


def _sample_config(config: Sim2RealPseudospotConfig, sample: str) -> SampleConfig:
    try:
        return config.samples[sample]
    except KeyError as exc:
        raise ValueError(f"Unknown sample {sample!r}; configured samples are {sorted(config.samples)}.") from exc


def _prepare_annotation(
    config: Sim2RealPseudospotConfig, sample: str, sample_config: SampleConfig
) -> AnnData:
    spatial = read_h5ad(sample_config.xenium_path)
    reference = filter_reference(
        read_h5ad(config.reference_path),
        filter_column=config.patient_column,
        filter_value=sample,
    )
    spatial = preprocess_spatial(
        spatial,
        min_transcript_counts=config.preprocessing.transcript_counts_min,
        min_cell_counts=config.preprocessing.gene_min_cells,
    )
    reference = preprocess_reference(
        reference,
        min_cell_counts=config.preprocessing.gene_min_cells,
    )
    overlap = spatial.var_names.intersection(reference.var_names)
    if overlap.empty:
        raise ValueError(f"{sample} has no overlapping genes after preprocessing.")
    spatial = spatial[:, overlap].copy()
    reference = reference[:, overlap].copy()
    if "spatial" not in spatial.obsm:
        if {"x", "y"}.issubset(spatial.obs.columns):
            spatial.obsm["spatial"] = spatial.obs[["x", "y"]].to_numpy()
        else:
            raise KeyError("Xenium input must provide obsm['spatial'] or obs[['x', 'y']].")
    if config.preprocessing.label_key not in reference.obs:
        raise KeyError(
            f"Reference is missing obs[{config.preprocessing.label_key!r}] for Level 1 annotation."
        )
    np.random.seed(config.seed)
    return OTKernel.annotate(
        spatial,
        reference,
        method="tacco",
        annotation_key=config.preprocessing.label_key,
        confidence_key=f"{config.preprocessing.label_key}_confidence",
        unknown_key=config.preprocessing.unknown_label,
    )


def _template_composition(
    config: Sim2RealPseudospotConfig, role: str
) -> dict[str, float]:
    part = REGION_PARTS[role]
    template_path = config.template_root / part / "selected_xenium.h5ad"
    if not template_path.is_file():
        raise FileNotFoundError(f"P2 template is missing: {template_path}")
    template = read_h5ad(template_path)
    retained = filter_region_cells(
        template,
        bounds=(
            float(np.asarray(template.obsm["spatial"])[:, 0].min()),
            float(np.asarray(template.obsm["spatial"])[:, 0].max()),
            float(np.asarray(template.obsm["spatial"])[:, 1].min()),
            float(np.asarray(template.obsm["spatial"])[:, 1].max()),
        ),
        label_key=config.preprocessing.label_key,
        unknown_label=config.preprocessing.unknown_label,
        min_cells_per_type=config.preprocessing.min_cells_per_type,
    )
    counts = retained.obs[config.preprocessing.label_key].astype(str).value_counts()
    if counts.empty:
        raise ValueError(f"P2 template {part} has no retained Level 1 cells.")
    return (counts / counts.sum()).to_dict()


def _candidate_record(candidate: Candidate) -> dict[str, Any]:
    return {
        "candidate_id": candidate.candidate_id,
        "bounds": [float(value) for value in candidate.bounds],
        "cell_count": candidate.cell_count,
        "scores": {
            "composition": candidate.composition_score,
            "simplicity": candidate.simplicity_score,
            "continuity": candidate.continuity_score,
            "coverage": candidate.coverage_score,
            "total": candidate.total_score,
        },
        "composition": dict(candidate.composition),
    }


def propose_regions(
    config: str | Path | Sim2RealPseudospotConfig, sample: str
) -> ProposalResult:
    """Annotate one sample and write three human-reviewable region proposals."""
    loaded = _as_config(config)
    sample_config = _sample_config(loaded, sample)
    region_dir = sample_config.output_dir / "region_selection"
    proposal_path = region_dir / "proposal.yaml"
    annotated_path = region_dir / "annotated_xenium.h5ad"
    if proposal_path.exists() or annotated_path.exists():
        raise FileExistsError(
            f"A proposal already exists in {region_dir}; do not overwrite a review artifact."
        )

    annotated = _prepare_annotation(loaded, sample, sample_config)
    candidates = {
        role: propose_candidates(
            annotated,
            role=role,
            template_composition=_template_composition(loaded, role),
            config=loaded.proposal,
            label_key=loaded.preprocessing.label_key,
            unknown_label=loaded.preprocessing.unknown_label,
            min_cells_per_type=loaded.preprocessing.min_cells_per_type,
        )
        for role in REGION_PARTS
    }
    short = [role for role, values in candidates.items() if len(values) != 3]
    if short:
        raise ValueError(f"Could not produce three valid candidates for: {', '.join(short)}")

    region_dir.mkdir(parents=True, exist_ok=False)
    image_path = region_dir / "candidate_regions.png"
    table_path = region_dir / "candidate_composition.csv"
    annotated.uns["sim2real_pseudospot"] = {
        "sample": sample,
        "config": str(loaded.path),
        "seed": loaded.seed,
    }
    annotated.write_h5ad(annotated_path)
    plot_candidates(
        annotated,
        label_key=loaded.preprocessing.label_key,
        candidates_by_role=candidates,
        output_path=image_path,
    )
    candidate_table(candidates).to_csv(table_path, index=False)
    proposal = {
        "schema_version": 1,
        "sample": sample,
        "config": str(loaded.path),
        "annotated_xenium": str(annotated_path),
        "regions": {
            role: [_candidate_record(candidate) for candidate in values]
            for role, values in candidates.items()
        },
    }
    proposal_path.write_text(yaml.safe_dump(proposal, sort_keys=False), encoding="utf-8")
    return ProposalResult(
        sample=sample,
        region_selection_dir=region_dir,
        annotated_path=annotated_path,
        image_path=image_path,
        table_path=table_path,
        proposal_path=proposal_path,
        candidates=candidates,
    )


def _bounds_from_confirmation(
    confirmation: Mapping[str, Any], proposal: Mapping[str, Any], sample: str
) -> dict[str, tuple[float, float, float, float]]:
    if confirmation.get("schema_version") != 1:
        raise ValueError("confirmation schema_version must be 1.")
    if confirmation.get("sample") != sample:
        raise ValueError("confirmation sample does not match the requested sample.")
    selections = confirmation.get("regions")
    if not isinstance(selections, Mapping):
        raise ValueError("confirmation regions must be a mapping.")
    bounds_by_role: dict[str, tuple[float, float, float, float]] = {}
    for role in REGION_PARTS:
        selected = selections.get(role)
        if not isinstance(selected, Mapping):
            raise ValueError(f"confirmation is missing region {role!r}.")
        if "bounds" in selected:
            raw_bounds = selected["bounds"]
        else:
            candidate_id = selected.get("candidate_id")
            matches = [
                candidate
                for candidate in proposal["regions"][role]
                if candidate["candidate_id"] == candidate_id
            ]
            if len(matches) != 1:
                raise ValueError(f"Unknown confirmed candidate {candidate_id!r} for {role}.")
            raw_bounds = matches[0]["bounds"]
        if not isinstance(raw_bounds, list) or len(raw_bounds) != 4:
            raise ValueError(f"confirmation bounds for {role} must have four values.")
        xmin, xmax, ymin, ymax = (float(value) for value in raw_bounds)
        if xmin >= xmax or ymin >= ymax:
            raise ValueError(f"confirmation bounds for {role} are invalid.")
        bounds_by_role[role] = (xmin, xmax, ymin, ymax)
    return bounds_by_role


def build_real_pseudospots(
    config: str | Path | Sim2RealPseudospotConfig,
    sample: str,
    confirmation: str | Path,
) -> BuildResult:
    """Build final real-only Sim2Real-ST outputs from confirmed region choices."""
    loaded = _as_config(config)
    sample_config = _sample_config(loaded, sample)
    confirmation_path = Path(confirmation).resolve()
    if not confirmation_path.is_file():
        raise FileNotFoundError(f"confirmation manifest does not exist: {confirmation_path}")
    with confirmation_path.open(encoding="utf-8") as handle:
        confirmed = yaml.safe_load(handle)
    if not isinstance(confirmed, Mapping):
        raise ValueError("confirmation manifest must be a mapping.")

    region_dir = sample_config.output_dir / "region_selection"
    proposal_path = region_dir / "proposal.yaml"
    if not proposal_path.is_file():
        raise FileNotFoundError(f"proposal manifest does not exist: {proposal_path}")
    with proposal_path.open(encoding="utf-8") as handle:
        proposal = yaml.safe_load(handle)
    if not isinstance(proposal, Mapping) or proposal.get("sample") != sample:
        raise ValueError("proposal manifest does not match the requested sample.")
    bounds_by_role = _bounds_from_confirmation(confirmed, proposal, sample)

    final_dir = sample_config.output_dir / "spot"
    temporary_dir = sample_config.output_dir / ".spot-building"
    if final_dir.exists() or temporary_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {final_dir}")
    annotated_path = Path(str(proposal.get("annotated_xenium", "")))
    if not annotated_path.is_file():
        raise FileNotFoundError(f"proposal annotated Xenium file is missing: {annotated_path}")
    annotated = read_h5ad(annotated_path)
    reference = filter_reference(
        read_h5ad(loaded.reference_path),
        filter_column=loaded.patient_column,
        filter_value=sample,
    )

    temporary_dir.mkdir(parents=True)
    parts: dict[str, Path] = {}
    for role, part in REGION_PARTS.items():
        selected = filter_region_cells(
            annotated,
            bounds=bounds_by_role[role],
            label_key=loaded.preprocessing.label_key,
            unknown_label=loaded.preprocessing.unknown_label,
            min_cells_per_type=loaded.preprocessing.min_cells_per_type,
        )
        if selected.n_obs == 0:
            raise ValueError(f"Confirmed {role} contains no retained cells.")
        selected.obs_names = selected.obs["cell_id"].astype(str)
        selected.obs_names.name = "cell_id"

        part_dir = temporary_dir / part
        part_dir.mkdir()
        selected.write_h5ad(part_dir / "selected_xenium.h5ad")
        plot_selected_region(
            selected,
            label_key=loaded.preprocessing.label_key,
            output_path=part_dir / "cut.png",
        )

        cell_types = selected.obs[loaded.preprocessing.label_key].astype(str).unique()
        ref_part = reference[
            reference.obs[loaded.preprocessing.label_key].astype(str).isin(cell_types)
        ].copy()
        ref_part.obs["clusters"] = ref_part.obs[loaded.preprocessing.label_key].astype(str)
        ref_part.write_h5ad(part_dir / "real_sc_ref_part.h5ad")

        for spot_size in loaded.spot_sizes:
            spots, distribution = aggregate_real_cells(selected, spot_size=spot_size)
            spots.uns["sim2real_pseudospot"] = {
                "sample": sample,
                "region": role,
                "spot_size": spot_size,
            }
            spot_dir = part_dir / f"spot_{spot_size}"
            spot_dir.mkdir()
            spots.write_h5ad(spot_dir / "xenium_spot.h5ad")
            distribution.to_csv(spot_dir / "cell_num_distribution.csv")
        parts[role] = final_dir / part

    temporary_dir.rename(final_dir)
    return BuildResult(sample=sample, output_dir=final_dir, parts=parts)
