from __future__ import annotations

import uuid
import warnings

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.ops.assignment import (
    GlobalAssignment,
    validate_global_assignment,
)
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.tacco_runtime import require_tacco
from revise.config.runner_conf import ApplicationScConf


class GlobalAnchoringKernel(BaseKernel):
    """Backend-native global anchoring kernel (POT/TACCO)."""

    def __init__(self, config, logger):
        super().__init__(config, logger)
        self.mode = self.config.annotate_mode
        self.cell_type_col = self.config.cell_type_col
        self.confidence_col = self.config.confidence_col

    @staticmethod
    def _reference_categories(reference_labels) -> pd.Index:
        return pd.Index(pd.unique(pd.Series(reference_labels, copy=False)))

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

            expected_categories = self._reference_categories(
                sc_ref_adata_overlap.obs[self.cell_type_col]
            )
            dums = pd.get_dummies(
                sc_ref_adata_overlap.obs[self.cell_type_col],
                dtype=sc_ref_adata_overlap.X.dtype,
            ).loc[:, expected_categories]
            ncats = dums.sum(axis=0)
            dums /= ncats.to_numpy()
            profiles = sc_ref_adata_overlap.X.T @ dums.to_numpy()
            profiles = pd.DataFrame(profiles, index=sc_ref_adata_overlap.var.index, columns=dums.columns)
            sc_ref_adata_overlap.varm[self.cell_type_col] = profiles
            # Keep parity numerical path: force dense ST matrix before
            # Bhattacharyya distance computation.
            st_dense = st_adata_overlap.X.toarray() if hasattr(st_adata_overlap.X, "toarray") else st_adata_overlap.X
            dist = bhattacharyya_distance(profiles.values.T, st_dense)

            cell_profile_mapping = pd.get_dummies(
                sc_ref_adata_overlap.obs[self.cell_type_col]
            ).loc[:, expected_categories]
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
            assignment = validate_global_assignment(
                GlobalAssignment(
                    labels=cell_type_ot.idxmax(axis=1),
                    posterior=cell_type_ot,
                ),
                expected_observations=st_adata.obs_names,
                expected_categories=expected_categories,
            )
            st_adata = self._publish_assignment(assignment, st_adata)
            return st_adata

        if self.mode == "tacco":
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
            st_adata_raw, _processed_reference = tacco_output
            if result_key not in st_adata_raw.obsm:
                raise KeyError(
                    f"TACCO did not write the fresh requested obsm[{result_key!r}]"
                )
            assignment = self._validate_tacco_result(
                st_adata_raw.obsm[result_key],
                st_adata_raw.obs_names,
                sc_ref_adata.obs[self.cell_type_col],
            )
            del st_adata_raw.obsm[result_key]
            st_adata.obsm = st_adata_raw.obsm.copy()
            st_adata = self._publish_assignment(assignment, st_adata)
            return st_adata

        raise NotImplementedError(f"Unsupported annotate_mode={self.mode}")

    def _validate_tacco_result(
        self,
        result,
        expected_obs_names,
        reference_labels,
    ) -> GlobalAssignment:
        if isinstance(result, pd.DataFrame):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                try:
                    labels = result.idxmax(axis=1)
                except (TypeError, ValueError):
                    labels = pd.Series(None, index=result.index, dtype=object)
        else:
            labels = pd.Series(dtype=object)
        return validate_global_assignment(
            GlobalAssignment(labels=labels, posterior=result),
            expected_observations=expected_obs_names,
            expected_categories=self._reference_categories(reference_labels),
        )

    def _publish_assignment(
        self,
        assignment: GlobalAssignment,
        st_adata: AnnData,
    ):
        st_adata.obsm[self.cell_type_col] = assignment.posterior
        st_adata.obs[self.cell_type_col] = assignment.labels.to_numpy(copy=True)
        st_adata.obs[self.confidence_col] = (
            assignment.posterior.max(axis=1).to_numpy(copy=True)
        )
        return st_adata
