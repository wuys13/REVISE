from __future__ import annotations

import uuid

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import ApplicationScConf


class GlobalAnchoringKernel(BaseKernel):
    """Backend-native global anchoring kernel (POT/TACCO)."""

    def __init__(self, config, logger, *, event_phase: str = "ga"):
        super().__init__(config, logger)
        if event_phase not in {"ga", "lr"}:
            raise ValueError("event_phase must be one of ['ga', 'lr']")
        self.mode = self.config.annotate_mode
        self.cell_type_col = self.config.cell_type_col
        self.confidence_col = self.config.confidence_col
        self.unknown_key = self.config.unknown_key
        self.event_callback = getattr(self.config, "ot_event_callback", None)
        self.event_phase = event_phase

    def _record_event(self, status: str) -> None:
        if self.event_callback is not None:
            self.event_callback(self.event_phase, self.mode, status)

    def _tacco_annotate_kwargs(self) -> dict:
        """Return the configured sc-SVC annotation parameters without fallbacks."""
        if not isinstance(self.config, ApplicationScConf):
            return {}
        multi_center = self.config.tacco_annotate_multi_center
        lamb = self.config.tacco_annotate_lamb
        missing = [
            name
            for name, value in (
                ("multi_center", multi_center),
                ("lamb", lamb),
            )
            if value is None
        ]
        if missing:
            raise ValueError(
                "Application sc-SVC TACCO annotation parameters are missing: "
                f"{missing}. Resolve them from sc.tacco_annotate."
            )
        return {"multi_center": multi_center, "lamb": lamb}

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

            overlap_genes = st_adata.var_names.intersection(sc_ref_adata.var_names)
            overlap_expression = st_adata[:, overlap_genes].X
            overlap_mass = np.asarray(overlap_expression.sum(axis=1)).ravel()
            if not np.all(np.isfinite(overlap_mass)):
                raise ValueError(
                    "TACCO target overlap-gene row sums must be finite"
                )
            if np.any(overlap_mass < 0):
                raise ValueError(
                    "TACCO target overlap-gene row sums must be non-negative"
                )
            if not np.any(overlap_mass > 0):
                raise ValueError(
                    "TACCO requires at least one target row with positive "
                    "overlap-gene mass"
                )

            # TACCO 0.5.0 removes target rows with no expression over the genes
            # retained by its own preprocessing, then reindexes the result to the
            # original observations and leaves those rows as NaN. Ask TACCO for
            # its processed reference so recovery uses the exact final gene space
            # and the same read-mass prior as its OT solver.
            st_adata_raw = st_adata.copy()
            tacco_reference_input = sc_ref_adata.copy()
            result_key = f"__revise_tacco_{uuid.uuid4().hex}"
            tacco_output = tc.tl.annotate(
                st_adata_raw,
                tacco_reference_input,
                self.cell_type_col,
                result_key=result_key,
                return_reference=True,
                **self._tacco_annotate_kwargs(),
            )
            if (
                not isinstance(tacco_output, tuple)
                or len(tacco_output) != 2
                or not isinstance(tacco_output[0], AnnData)
                or not isinstance(tacco_output[1], AnnData)
            ):
                raise ValueError(
                    "TACCO return_reference=True must return "
                    "(annotated_target, processed_reference)"
                )
            st_adata_raw, processed_reference = tacco_output
            if result_key not in st_adata_raw.obsm:
                raise KeyError(
                    f"TACCO did not write the fresh requested obsm[{result_key!r}]"
                )
            cell_type_candidate = self._normalize_tacco_result_axes(
                st_adata_raw.obsm[result_key],
                st_adata_raw.obs_names,
                sc_ref_adata.obs[self.cell_type_col],
            )

            final_genes = st_adata_raw.var_names.intersection(
                processed_reference.var_names
            )
            final_expression = st_adata_raw[:, final_genes].X
            final_mass = np.asarray(final_expression.sum(axis=1)).ravel()
            if not np.all(np.isfinite(final_mass)):
                raise ValueError(
                    "TACCO target final-gene row sums must be finite"
                )
            if np.any(final_mass < 0):
                raise ValueError(
                    "TACCO target final-gene row sums must be non-negative"
                )
            if not np.any(final_mass > 0):
                raise ValueError(
                    "TACCO requires at least one target row with positive "
                    "final-gene mass"
                )

            try:
                candidate_values = cell_type_candidate.to_numpy(
                    dtype=np.float64,
                    copy=True,
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "TACCO result values must be numeric"
                ) from exc
            all_nan_rows = np.isnan(candidate_values).all(axis=1)
            recoverable_rows = all_nan_rows & (final_mass == 0)

            if np.any(recoverable_rows):
                reference_labels = pd.Series(
                    processed_reference.obs[self.cell_type_col],
                    copy=False,
                )
                if reference_labels.empty or reference_labels.isna().any():
                    raise ValueError(
                        "TACCO processed reference labels must be non-null"
                    )
                reference_mapping = pd.get_dummies(
                    reference_labels,
                    dtype=np.float64,
                )
                reference_columns = pd.Index(
                    [str(value) for value in reference_mapping.columns]
                )
                if not reference_columns.is_unique:
                    raise ValueError(
                        "TACCO processed reference labels collide after "
                        "conversion to strings"
                    )
                reference_mapping.columns = reference_columns
                reference_mass = np.asarray(
                    processed_reference.X.sum(axis=1)
                ).ravel()
                if not np.all(np.isfinite(reference_mass)):
                    raise ValueError(
                        "TACCO processed reference row sums must be finite"
                    )
                if np.any(reference_mass < 0):
                    raise ValueError(
                        "TACCO processed reference row sums must be "
                        "non-negative"
                    )
                reference_prior = pd.Series(
                    reference_mass @ reference_mapping.to_numpy(),
                    index=reference_mapping.columns,
                    dtype=np.float64,
                ).reindex(cell_type_candidate.columns)
                prior_sum = reference_prior.sum()
                if (
                    reference_prior.isna().any()
                    or not np.isfinite(prior_sum)
                    or prior_sum <= 0
                ):
                    raise ValueError(
                        "TACCO processed reference prior must have positive "
                        "finite mass for every result category"
                    )
                reference_prior /= prior_sum
                cell_type_candidate.loc[recoverable_rows, :] = (
                    reference_prior.to_numpy()
                )
                self.logger.warning(
                    "TACCO recovered %d non-finite target rows with zero "
                    "final-gene mass using its processed-reference read-mass "
                    "prior",
                    int(recoverable_rows.sum()),
                )

            cell_type_ot = self._validate_tacco_result(
                cell_type_candidate,
                st_adata_raw.obs_names,
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
        validated = self._normalize_tacco_result_axes(
            result,
            expected_obs_names,
            reference_labels,
        )
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

    def _normalize_tacco_result_axes(
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
        return validated

    def _assign_result(self, cell_type_probability: pd.DataFrame, st_adata: AnnData):
        max_columns = cell_type_probability.idxmax(axis=1)
        max_values = cell_type_probability.max(axis=1)
        labels = max_columns.to_numpy(dtype=object, copy=True)
        labels[pd.isna(labels)] = self.unknown_key
        st_adata.obs[self.cell_type_col] = labels
        st_adata.obs[self.confidence_col] = max_values.values.copy()
        return st_adata
