from __future__ import annotations

import uuid

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.tacco_runtime import require_tacco


class GlobalAnchoringKernel(BaseKernel):
    """Backend-native global anchoring kernel (POT/TACCO)."""

    def __init__(self, config, logger):
        super().__init__(config, logger)
        self.mode = self.config.annotate_mode
        self.cell_type_col = self.config.cell_type_col
        self.confidence_col = self.config.confidence_col
        self.unknown_key = self.config.unknown_key
        self.event_callback = getattr(self.config, "ot_event_callback", None)

    def _record_event(self, status: str) -> None:
        if self.event_callback is not None:
            self.event_callback("ga", self.mode, status)

    @staticmethod
    def _validate_coupling(result, expected_shape) -> np.ndarray:
        coupling = np.asarray(result)
        if coupling.shape != expected_shape:
            raise ValueError(
                f"Global OT coupling has shape {coupling.shape}, expected {expected_shape}"
            )
        try:
            finite = np.isfinite(coupling)
        except TypeError as exc:
            raise ValueError("Global OT coupling values must be numeric") from exc
        if not np.all(finite):
            raise ValueError("Global OT coupling values must be finite")
        if np.any(coupling < 0):
            raise ValueError("Global OT coupling values must be non-negative")
        if np.any(coupling.sum(axis=1) <= 0):
            raise ValueError("Global OT coupling rows must have positive mass")
        return coupling

    def run(self, st_adata: AnnData, sc_ref_adata: AnnData, **kwargs):
        if "cell_type_col" in kwargs:
            self.cell_type_col = kwargs["cell_type_col"]
        st_adata = st_adata.copy()

        if self.mode == "pot":
            required = ("annotate_pot_reg", "annotate_pot_reg_m", "annotate_pot_reg_type")
            if any(k not in kwargs for k in required):
                raise ValueError(f"mode={self.mode} requires {required}")

            overlap_genes = list(st_adata.var_names.intersection(sc_ref_adata.var_names))
            sc_ref_adata_overlap = sc_ref_adata[:, overlap_genes].copy()
            st_adata_overlap = st_adata[:, overlap_genes].copy()

            dums = pd.get_dummies(sc_ref_adata_overlap.obs[self.cell_type_col], dtype=sc_ref_adata_overlap.X.dtype)
            ncats = dums.sum(axis=0)
            dums /= ncats.to_numpy()
            profiles = sc_ref_adata_overlap.X.T @ dums.to_numpy()
            profiles = pd.DataFrame(profiles, index=sc_ref_adata_overlap.var.index, columns=dums.columns)
            sc_ref_adata_overlap.varm[self.cell_type_col] = profiles
            # Keep parity numerical path: force dense ST matrix before
            # Bhattacharyya distance computation.
            st_dense = st_adata_overlap.X.toarray() if hasattr(st_adata_overlap.X, "toarray") else st_adata_overlap.X
            dist = bhattacharyya_distance(profiles.values.T, st_dense)

            cell_profile_mapping = pd.get_dummies(sc_ref_adata_overlap.obs[self.cell_type_col])
            cell_profile_mapping /= cell_profile_mapping.sum(axis=1).to_numpy()[:, None]
            type_prior = np.array(sc_ref_adata_overlap.X.sum(axis=1)).flatten() @ cell_profile_mapping
            spot_prior = pd.Series(np.array(st_adata_overlap.X.sum(axis=1)).flatten(), index=st_adata_overlap.obs.index)

            spot_prior_sum = spot_prior.sum()
            type_prior_sum = type_prior.sum()
            if spot_prior_sum > 0:
                spot_prior /= spot_prior_sum
            else:
                spot_prior = spot_prior / len(spot_prior)
            if type_prior_sum > 0:
                type_prior /= type_prior_sum
            else:
                type_prior = type_prior / len(type_prior)

            dist = np.nan_to_num(dist)
            dist_max = dist.max()
            if dist_max <= 0:
                dist_max = 1.0
            dist = dist.T / dist_max

            self.logger.info(
                "pot params: reg=%s, reg_m=%s, reg_type=%s",
                kwargs["annotate_pot_reg"],
                kwargs["annotate_pot_reg_m"],
                kwargs["annotate_pot_reg_type"],
            )
            self._record_event("attempted")
            import ot

            t_transform = ot.unbalanced.sinkhorn_unbalanced(
                spot_prior.values,
                type_prior.values,
                dist,
                reg=kwargs["annotate_pot_reg"],
                reg_m=kwargs["annotate_pot_reg_m"],
                reg_type=kwargs["annotate_pot_reg_type"],
                verbose=True,
                numItermax=5000,
            )
            t_transform = self._validate_coupling(
                t_transform,
                (len(spot_prior), len(type_prior)),
            )
            cell_type_ot = pd.DataFrame(t_transform, index=spot_prior.index, columns=type_prior.index)
            cell_type_ot *= (cell_type_ot > 0)
            cell_type_ot /= cell_type_ot.sum(axis=1).to_numpy()[:, None]
            st_adata.obsm[self.cell_type_col] = cell_type_ot
            st_adata = self._assign_result(cell_type_ot, st_adata)
            self._record_event("completed")
            return st_adata

        if self.mode == "tacco":
            self._record_event("attempted")
            tc = require_tacco()

            st_adata_raw = st_adata.copy()
            result_key = f"__revise_tacco_{uuid.uuid4().hex}"
            tc.tl.annotate(
                st_adata_raw,
                sc_ref_adata,
                self.cell_type_col,
                result_key=result_key,
            )
            if result_key not in st_adata_raw.obsm:
                raise KeyError(
                    f"TACCO did not write the fresh requested obsm[{result_key!r}]"
                )
            cell_type_ot = self._validate_tacco_result(
                st_adata_raw.obsm[result_key],
                st_adata.obs_names,
                sc_ref_adata.obs[self.cell_type_col],
            )
            del st_adata_raw.obsm[result_key]
            st_adata.obsm = st_adata_raw.obsm.copy()
            st_adata.obsm[self.cell_type_col] = cell_type_ot
            st_adata = self._assign_result(cell_type_ot, st_adata)
            self._record_event("completed")
            return st_adata

        raise NotImplementedError(f"Unsupported annotate_mode={self.mode}")

    def _validate_tacco_result(
        self,
        result,
        expected_obs_names,
        reference_labels,
    ) -> pd.DataFrame:
        if not isinstance(result, pd.DataFrame):
            raise ValueError("TACCO annotation result must be a pandas DataFrame")

        expected_obs = pd.Index(expected_obs_names)
        if expected_obs.hasnans:
            raise ValueError("Expected ST observation names must not contain nulls")
        if not expected_obs.is_unique:
            raise ValueError("Expected ST observation names must be unique")

        result_obs = pd.Index(result.index)
        if result_obs.hasnans:
            raise ValueError("TACCO result observation names must not contain nulls")
        if not result_obs.is_unique:
            raise ValueError("TACCO result observation names must be unique")
        missing_obs = expected_obs.difference(result_obs)
        extra_obs = result_obs.difference(expected_obs)
        if len(missing_obs) or len(extra_obs):
            raise ValueError(
                "TACCO result observation mismatch: "
                f"missing={missing_obs.tolist()}, extra={extra_obs.tolist()}"
            )

        reference_values = pd.Series(reference_labels, copy=False).dropna()
        if reference_values.empty:
            raise ValueError("TACCO reference labels must contain an observed value")
        observed_reference = pd.Index(pd.unique(reference_values))
        reference_columns = pd.Index([str(value) for value in observed_reference])
        if not reference_columns.is_unique:
            raise ValueError(
                "TACCO reference labels collide after conversion to strings"
            )

        result_columns = pd.Index(result.columns)
        if result_columns.hasnans:
            raise ValueError("TACCO result categories must not contain nulls")
        if not result_columns.is_unique:
            raise ValueError("TACCO result categories must be unique")
        string_columns = pd.Index([str(value) for value in result_columns])
        if not string_columns.is_unique:
            raise ValueError(
                "TACCO result categories collide after conversion to strings"
            )

        missing_categories = reference_columns.difference(string_columns)
        extra_categories = string_columns.difference(reference_columns)
        if len(missing_categories) or len(extra_categories):
            raise ValueError(
                "TACCO result category mismatch: "
                f"missing={missing_categories.tolist()}, "
                f"extra={extra_categories.tolist()}"
            )

        validated = result.copy()
        validated.columns = string_columns
        validated = validated.loc[expected_obs, reference_columns].copy()
        values = validated.to_numpy()
        try:
            finite = np.isfinite(values)
        except TypeError as exc:
            raise ValueError("TACCO result values must be numeric") from exc
        if not np.all(finite):
            raise ValueError("TACCO result values must be finite")
        if np.any(values < 0):
            raise ValueError("TACCO result values must be non-negative")
        if np.any(values.sum(axis=1) <= 0):
            raise ValueError("TACCO result rows must have positive mass")
        return validated

    def _assign_result(self, cell_type_probability: pd.DataFrame, st_adata: AnnData):
        max_columns = cell_type_probability.idxmax(axis=1)
        max_values = cell_type_probability.max(axis=1)
        labels = max_columns.to_numpy(dtype=object, copy=True)
        labels[pd.isna(labels)] = self.unknown_key
        st_adata.obs[self.cell_type_col] = labels
        st_adata.obs[self.confidence_col] = max_values.values.copy()
        return st_adata
