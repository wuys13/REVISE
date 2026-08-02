from __future__ import annotations

import copy

import numpy as np
import pandas as pd
from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.kernels.global_anchoring import GlobalAnchoringKernel
from revise.backend.ops.assignment import AssignmentState
from revise.backend.ops.distance import bhattacharyya_distance
from revise.backend.ops.local_ot import solve_local_ot
from revise.config.runner_conf import ApplicationScConf


class LocalAnchoringKernel(BaseKernel):
    """Annotate one sc-SVC local unit through the configured local OT solver."""

    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.method = str(self.config.rec_ot_method).strip().lower()
        self.confidence_col = self.config.confidence_col
        self.unknown_key = self.config.unknown_key

    def run(self, target: AnnData, reference: AnnData, **kwargs) -> AnnData:
        cell_type_col = kwargs.get("cell_type_col", self.config.cell_type_col)
        if isinstance(self.config, ApplicationScConf) and self.method == "tacco":
            delegate_config = copy.copy(self.config)
            delegate_config.annotate_mode = "tacco"
            return GlobalAnchoringKernel(
                delegate_config,
                self.logger,
                event_phase="lr",
            ).run(
                target,
                reference,
                cell_type_col=cell_type_col,
            )
        target = target.copy()
        overlap_genes = list(target.var_names.intersection(reference.var_names))
        if not overlap_genes:
            raise ValueError("Local anchoring inputs have no overlapping genes")

        target_overlap = target[:, overlap_genes]
        reference_overlap = reference[:, overlap_genes]
        labels = pd.get_dummies(
            reference_overlap.obs[cell_type_col],
            dtype=reference_overlap.X.dtype,
        )
        label_counts = labels.sum(axis=0)
        profiles = reference_overlap.X.T @ (labels / label_counts).to_numpy()
        target_matrix = target_overlap.X
        if hasattr(target_matrix, "toarray"):
            target_matrix = target_matrix.toarray()
        cost = bhattacharyya_distance(profiles.T, target_matrix)
        cost = np.nan_to_num(cost)
        cost_max = float(cost.max()) if cost.size else 0.0
        if cost_max <= 0:
            cost_max = 1.0

        target_mass = np.asarray(target_overlap.X.sum(axis=1)).ravel()
        reference_mass = np.asarray(reference_overlap.X.sum(axis=1)).ravel()
        label_prior = np.asarray(reference_mass @ labels).ravel()
        target_total = float(target_mass.sum())
        label_total = float(label_prior.sum())
        if target_total <= 0 or label_total <= 0:
            raise ValueError("Local anchoring marginals must have positive total mass")
        target_mass = target_mass / target_total
        label_prior = label_prior / label_total
        coupling = solve_local_ot(
            target_mass,
            label_prior,
            cost.T / cost_max,
            method=self.method,
            pot_reg=self.config.rec_pot_reg,
            pot_reg_m=self.config.rec_pot_reg_m,
            pot_reg_type=self.config.rec_pot_reg_type,
            pot_verbose=False,
            pot_num_iter_max=5000,
            event_callback=getattr(self.config, "ot_event_callback", None),
        )

        row_mass = coupling.sum(axis=1)
        if np.any(row_mass <= 0):
            raise ValueError("Local anchoring coupling rows must have positive mass")
        probabilities = pd.DataFrame(
            coupling / row_mass[:, None],
            index=target.obs_names,
            columns=labels.columns,
        )
        target.obsm[cell_type_col] = probabilities
        target.obs[cell_type_col] = probabilities.idxmax(axis=1).to_numpy(dtype=object)
        target.obs[self.confidence_col] = probabilities.max(axis=1).to_numpy()
        target.obs[cell_type_col] = target.obs[cell_type_col].replace(
            {np.nan: self.unknown_key}
        )
        return target

    @staticmethod
    def assignment_state(
        annotated: AnnData,
        key: str,
    ) -> AssignmentState | None:
        """Package this anchoring call's labeled soft output for policy validation."""
        if key not in annotated.obsm:
            return None
        probabilities = annotated.obsm[key]
        if not isinstance(probabilities, pd.DataFrame):
            return None
        return AssignmentState(
            values=probabilities.to_numpy(dtype=np.float64),
            observation_labels=probabilities.index,
            category_labels=probabilities.columns,
            source=f"local_anchoring:obsm[{key}]",
            level=str(key),
            value_semantics="soft",
            lineage=[
                {
                    "operation": "local_anchoring",
                    "container": "obsm",
                    "key": str(key),
                }
            ],
        )
