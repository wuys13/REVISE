from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad

from revise.io.input_bundle import REVISEDataBundle
from revise.io.spatialdata_service import SpatialDataService


@dataclass
class _LoadedAnnData:
    adata: AnnData
    sdata: Any = None
    metadata: Dict[str, Any] = field(default_factory=dict)


class REVISEInputService:
    """Normalize user inputs into the AnnData contract used by REVISE.

    The service deliberately keeps H5AD as the default read path. SpatialData is
    activated only for ST data when requested explicitly, or when `auto` sees a
    `.zarr` SpatialData store. This preserves existing benchmark/application
    behavior while adding a clean compatibility point for the spatial ecosystem.
    """

    def __init__(self, io_config: Optional[Dict[str, Any]] = None, logger=None):
        self.io_config = dict(io_config or {})
        self.logger = logger

    @classmethod
    def from_context(cls, ctx) -> "REVISEInputService":
        return cls(io_config=ctx.io, logger=ctx.logger)

    def read_st_adata(self, path: Union[str, Path]) -> AnnData:
        return self._read_role(path, role="st").adata

    def read_sc_ref_adata(self, path: Union[str, Path]) -> AnnData:
        return self._read_role(path, role="sc_ref").adata

    def read_real_adata(self, path: Union[str, Path]) -> AnnData:
        return self._read_role(path, role="gt").adata

    def load_bundle(
        self,
        *,
        st_path: Union[str, Path],
        sc_ref_path: Union[str, Path],
        real_st_path: Optional[Union[str, Path]] = None,
    ) -> REVISEDataBundle:
        st_loaded = self._read_role(st_path, role="st")
        sc_loaded = self._read_role(sc_ref_path, role="sc_ref")
        real_loaded = self._read_role(real_st_path, role="gt") if real_st_path is not None else None
        source_report = {
            "st": st_loaded.metadata,
            "sc_ref": sc_loaded.metadata,
        }
        if real_loaded is not None:
            source_report["gt"] = real_loaded.metadata
        return REVISEDataBundle(
            st_adata=st_loaded.adata,
            sc_ref_adata=sc_loaded.adata,
            real_st_adata=real_loaded.adata if real_loaded is not None else None,
            sdata=st_loaded.sdata,
            coordinate_system=self.io_config.get("spatialdata_coordinate_system", "global"),
            source_report=source_report,
        )

    def preflight(self, specs, *, runtime, columns) -> Dict[str, Any]:
        """Validate route-required input metadata without loading expression data."""
        opened: Dict[str, AnnData] = {}
        reports = []
        with ExitStack() as stack:
            for spec in specs:
                role = str(spec.role)
                path = Path(spec.path)
                if not path.exists():
                    raise FileNotFoundError(
                        f"Invalid input: role={role}; path={path}; "
                        "expected=readable input; actual=missing"
                    )
                try:
                    loaded = self._read_role(path, role=role, backed=True)
                except Exception as exc:
                    raise ValueError(
                        f"Invalid input: role={role}; path={path}; "
                        f"expected=open input; actual={type(exc).__name__}: {exc}"
                    ) from exc
                adata = loaded.adata
                input_format = str(loaded.metadata["format"])
                backed = bool(getattr(adata, "isbacked", False))
                if backed:
                    stack.callback(adata.file.close)

                self._validate_metadata(
                    adata,
                    role=role,
                    path=path,
                    runtime=runtime,
                    columns=columns,
                )
                opened[role] = adata
                input_report = {
                    "role": role,
                    "path": str(path),
                    "format": input_format,
                    "backed": backed,
                    "shape": [int(adata.n_obs), int(adata.n_vars)],
                }
                if role == "st" and str(runtime.get("task")) == "sc_svc_sr":
                    raw_mapping = adata.uns.get("all_cells_in_spot")
                    if raw_mapping is not None:
                        mapping_source = "embedded"
                    elif str(runtime.get("mode")) == "benchmark":
                        mapping_source = "generated_from_ground_truth_coordinates"
                    else:
                        mapping_source = "generated_from_transcript_counts"
                    input_report["sr_mapping"] = {
                        "source": mapping_source,
                        "validation": "pre_allocation",
                    }
                if role == "gt" and str(runtime.get("task")) == "sc_svc_sr":
                    input_report["ground_truth_label_source"] = (
                        self._resolve_sr_ground_truth_label_key(
                            adata,
                            str(columns.get("cell_type_col", "Level1")),
                        )
                    )
                reports.append(input_report)

            st = opened.get("st")
            sc_ref = opened.get("sc_ref")
            if st is None or sc_ref is None:
                raise ValueError("Resolved inputs must include st and sc_ref roles")
            overlap = st.var_names.intersection(sc_ref.var_names)
            paths = {str(spec.role): str(spec.path) for spec in specs}
            if overlap.empty:
                raise ValueError(
                    "Invalid input gene overlap: expected=>=1; actual=0; "
                    f"st_path={paths['st']}; sc_ref_path={paths['sc_ref']}"
                )

            self._validate_cross_role_contracts(
                opened,
                runtime=runtime,
                paths=paths,
            )
            return {
                "status": "ready",
                "inputs": reports,
                "gene_overlap": int(len(overlap)),
                "proof_boundary": (
                    "metadata_and_required_arrays_only; "
                    "expression_values_not_fully_scanned"
                ),
            }

    def _validate_metadata(
        self,
        adata: AnnData,
        *,
        role: str,
        path: Path,
        runtime,
        columns,
    ) -> None:
        context = f"role={role}; path={path}"
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(
                f"Invalid input: {context}; field=shape; expected=nonempty; "
                f"actual={adata.shape}"
            )
        if not adata.obs_names.is_unique:
            raise ValueError(
                f"Invalid input: {context}; field=obs_names; expected=unique"
            )
        if not adata.var_names.is_unique:
            raise ValueError(
                f"Invalid input: {context}; field=var_names; expected=unique"
            )

        mode = str(runtime.get("mode"))
        task = str(runtime.get("task"))
        if role == "st":
            if not (mode == "benchmark" and task == "sc_svc_impute"):
                self._validate_spatial(adata, context=context)
            required_obs = []
            if mode == "benchmark" and task == "sp_svc":
                required_obs.append("seg_error")
            if mode == "benchmark" and task == "sc_svc_impute":
                required_obs.append("transcript_counts")
            self._require_obs(adata, required_obs, context=context)
        elif role == "sc_ref":
            required_obs = [str(columns.get("cell_type_col", "Level1"))]
            if mode == "application" and task == "sc_svc":
                required_obs.append(
                    str(columns.get("sub_cell_type_col", "Level2"))
                )
            self._require_reference_labels(
                adata,
                list(dict.fromkeys(required_obs)),
                context=context,
            )
            patient_key = self.io_config.get("patient_key")
            sample_name = self.io_config.get("sample_name")
            requires_patient_match = not (
                mode == "benchmark"
                and str(runtime.get("confounding")) == "batch_effect"
            )
            if (
                requires_patient_match
                and patient_key
                and patient_key in adata.obs
                and sample_name is not None
            ):
                patient_sample = str(sample_name)
                if mode == "benchmark":
                    patient_sample = patient_sample.split("/", 1)[0]
                if not adata.obs[patient_key].astype(str).eq(patient_sample).any():
                    raise ValueError(
                        f"Invalid input: {context}; field=obs[{patient_key!r}]; "
                        f"expected=at least one row for sample {patient_sample!r}"
                    )
        elif role == "gt" and task == "sc_svc_sr":
            label_key = self._resolve_sr_ground_truth_label_key(
                adata,
                str(columns.get("cell_type_col", "Level1")),
            )
            self._require_obs(
                adata,
                ["cell_id", "x", "y"],
                context=context,
            )
            self._require_reference_labels(
                adata,
                [label_key],
                context=context,
            )
            cell_ids = adata.obs["cell_id"]
            normalized_ids = cell_ids.astype(str)
            if cell_ids.isna().any() or normalized_ids.str.strip().eq("").any():
                raise ValueError(
                    f"Invalid input: {context}; field=obs['cell_id']; "
                    "expected=nonempty values"
                )
            if normalized_ids.duplicated().any():
                raise ValueError(
                    f"Invalid input: {context}; field=obs['cell_id']; "
                    "expected=unique"
                )
            if adata.obs[label_key].isna().any():
                raise ValueError(
                    f"Invalid input: {context}; field=obs[{label_key!r}]; "
                    "expected=non-null"
                )
            try:
                coordinates = adata.obs.loc[:, ["x", "y"]].to_numpy(
                    dtype=np.float64
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid input: {context}; field=obs[['x', 'y']]; "
                    "expected=numeric finite coordinates"
                ) from exc
            if not np.all(np.isfinite(coordinates)):
                raise ValueError(
                    f"Invalid input: {context}; field=obs[['x', 'y']]; "
                    "expected=finite coordinates"
                )
            if "spatial" in adata.obsm:
                self._validate_spatial(adata, context=context)

    @staticmethod
    def _resolve_sr_ground_truth_label_key(
        adata: AnnData,
        configured_key: str,
    ) -> str:
        if configured_key != "Level1" or configured_key in adata.obs:
            return configured_key
        if "clusters" in adata.obs:
            return "clusters"
        return configured_key

    @staticmethod
    def _require_obs(adata: AnnData, required, *, context: str) -> None:
        missing = [name for name in required if name not in adata.obs]
        if missing:
            raise KeyError(
                f"Invalid input: {context}; field=obs; missing={missing}"
            )

    @staticmethod
    def _require_reference_labels(
        adata: AnnData,
        required,
        *,
        context: str,
    ) -> None:
        REVISEInputService._require_obs(adata, required, context=context)
        for name in required:
            original = adata.obs[name].astype(str)
            normalized = original.str.replace("/", "_", regex=False)
            label_pairs = pd.DataFrame(
                {"original": original, "normalized": normalized}
            ).drop_duplicates()
            collisions = label_pairs.groupby(
                "normalized",
                sort=False,
            )["original"].nunique()
            if (collisions > 1).any():
                names = collisions[collisions > 1].index.tolist()
                raise ValueError(
                    f"Invalid input: {context}; field=obs[{name!r}]; "
                    "expected=labels must not collide after '/' normalization; "
                    f"actual_collisions={names[:5]}"
                )

    @staticmethod
    def _validate_spatial(adata: AnnData, *, context: str) -> None:
        if "spatial" not in adata.obsm:
            raise KeyError(
                f"Invalid input: {context}; field=obsm['spatial']; expected=present"
            )
        coordinates = np.asarray(adata.obsm["spatial"])
        if (
            coordinates.ndim != 2
            or coordinates.shape[0] != adata.n_obs
            or coordinates.shape[1] < 2
        ):
            raise ValueError(
                f"Invalid input: {context}; field=obsm['spatial']; "
                f"expected=shape (n_obs, >=2); actual={coordinates.shape}"
            )
        if not np.all(np.isfinite(coordinates)):
            raise ValueError(
                f"Invalid input: {context}; field=obsm['spatial']; "
                "expected=finite coordinates"
            )

    @staticmethod
    def _validate_cross_role_contracts(opened, *, runtime, paths) -> None:
        mode = str(runtime.get("mode"))
        task = str(runtime.get("task"))
        if mode != "benchmark":
            return
        if task == "sc_svc_sr":
            return
        if task not in {"sp_svc", "sc_svc_impute"}:
            return
        st = opened["st"]
        gt = opened.get("gt")
        if gt is not None and st.obs_names.intersection(gt.obs_names).empty:
            raise ValueError(
                "Invalid benchmark alignment: role=gt; "
                f"path={paths['gt']}; field=obs_names_overlap; expected=>=1; "
                f"actual=0; st_path={paths['st']}"
            )

    def _read_role(
        self,
        path: Union[str, Path],
        *,
        role: str,
        backed: bool = False,
    ) -> _LoadedAnnData:
        input_format = self._resolve_input_format(path, role=role)
        if input_format == "spatialdata":
            return self._read_spatialdata_st(path)
        if input_format != "h5ad":
            raise ValueError(
                f"Unsupported input format {input_format!r} for role {role!r}; "
                "expected 'h5ad', 'spatialdata', or 'auto'."
            )
        adata = read_h5ad(path, backed="r") if backed else read_h5ad(path)
        return _LoadedAnnData(
            adata=adata,
            metadata={
                "format": "h5ad",
                "path": str(path),
                "role": role,
            },
        )

    def _resolve_input_format(self, path: Union[str, Path], *, role: str) -> str:
        if role != "st":
            return str(self.io_config.get(f"{role}_input_format", "h5ad")).lower()

        configured = str(self.io_config.get("input_format", "h5ad")).lower()
        if configured == "auto":
            spatial_path = self._spatialdata_path(path)
            if self._looks_like_spatialdata_store(spatial_path):
                return "spatialdata"
            return "h5ad"
        return configured

    def _read_spatialdata_st(self, path: Union[str, Path]) -> _LoadedAnnData:
        spatial_path = self._spatialdata_path(path)
        reader = str(self.io_config.get("spatialdata_reader", "zarr")).lower()
        if reader not in {"zarr", "auto"}:
            raise ValueError(
                "Only SpatialData Zarr stores are supported in this initial service layer. "
                f"Received io.spatialdata_reader={reader!r}."
            )

        result = SpatialDataService(logger=self.logger).read_table(
            spatial_path,
            table_key=self.io_config.get("spatialdata_table"),
            spatial_element_key=self.io_config.get("spatialdata_spatial_element"),
            coordinate_system=str(self.io_config.get("spatialdata_coordinate_system", "global")),
        )
        return _LoadedAnnData(adata=result.adata, sdata=result.sdata, metadata=result.metadata)

    def _spatialdata_path(self, fallback_path: Union[str, Path]) -> Union[str, Path]:
        return self.io_config.get("spatialdata_path") or fallback_path

    @staticmethod
    def _looks_like_spatialdata_store(path: Union[str, Path]) -> bool:
        path_str = str(path)
        return path_str.endswith(".zarr")
