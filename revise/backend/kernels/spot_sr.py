import os

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from revise.backend.kernels.base import BaseKernel


class SpotSrKernel(BaseKernel):
    """Assign spot-level cell-type quotas to virtual cells."""

    def __init__(self, config, logger):
        super().__init__(config, logger)
        if os.path.exists(self.config.pm_on_cell_file):
            self.pm_on_cell = pd.read_csv(self.config.pm_on_cell_file, index_col=0)
        else:
            self.pm_on_cell = None

    def run(self, sc_svc):
        cell_type_col = str(getattr(self.config, "cell_type_col", "Level1"))
        spot_cell_distribution = self.get_spot_cell_distribution(
            cell_contributions=sc_svc.st_adata.obsm[cell_type_col],
            SVC_obs=sc_svc.svc_obs,
        )
        if self.pm_on_cell is not None:
            assigned = self.assign_cell_types(sc_svc.svc_obs, spot_cell_distribution)
        else:
            self.logger.warning("No pm_on_cell file; assigning quota slots randomly")
            assigned = self.assign_cell_types_random(sc_svc.svc_obs, spot_cell_distribution)

        assigned["match"] = assigned["cell_type"] == assigned["true_cell_type"]
        match_rate = float(assigned["match"].mean())
        self.logger.info("%s \n %s", assigned["match"].value_counts(), match_rate)
        sc_svc.svc_obs = assigned

    def get_spot_cell_distribution(self, cell_contributions, SVC_obs):
        """Convert validated spot proportions into exact integer cell quotas."""
        if self.config.svc_completeness is not True:
            raise ValueError("svc_completeness must be exactly true")
        if not isinstance(cell_contributions, pd.DataFrame):
            raise TypeError("cell_contributions must be a pandas DataFrame")
        if cell_contributions.index.has_duplicates:
            raise ValueError("cell_contributions must have unique spot labels")
        if cell_contributions.columns.has_duplicates:
            raise ValueError("cell_contributions must have unique category labels")
        if not isinstance(SVC_obs, pd.DataFrame) or not {"spot_name", "cell_id"} <= set(SVC_obs.columns):
            raise ValueError("SVC_obs must contain spot_name and cell_id columns")

        contribution_spots = set(cell_contributions.index)
        svc_spots = set(pd.unique(SVC_obs["spot_name"]))
        if contribution_spots != svc_spots:
            missing = sorted(svc_spots - contribution_spots, key=str)
            extra = sorted(contribution_spots - svc_spots, key=str)
            raise ValueError(
                f"cell_contributions spot set must exactly match SVC_obs; missing={missing}, extra={extra}"
            )

        try:
            values = cell_contributions.to_numpy(dtype=np.float64, copy=True)
        except (TypeError, ValueError) as exc:
            raise ValueError("cell_contributions values must be convertible to float64") from exc
        if not np.isfinite(values).all():
            raise ValueError("cell_contributions values must be finite")
        if np.any(values < 0):
            raise ValueError("cell_contributions values must be nonnegative")
        row_sums = values.sum(axis=1)
        if not np.all(np.isclose(row_sums, 1.0)):
            raise ValueError("each cell_contributions row must sum to 1")

        spot_cell_counts = SVC_obs.groupby("spot_name", sort=False)["cell_id"].size().to_dict()
        quotas = np.empty(values.shape, dtype=np.int64)
        for row_idx, spot_name in enumerate(cell_contributions.index):
            total_cells = int(spot_cell_counts[spot_name])
            counts = np.round(values[row_idx] * total_cells).astype(np.int64)
            adjustment = total_cells - int(counts.sum())
            if adjustment > 0:
                counts[int(np.argmax(counts))] += adjustment
                adjustment = 0
            elif adjustment < 0:
                while adjustment < 0:
                    previous_adjustment = adjustment
                    positive = np.flatnonzero(counts > 0)
                    order = positive[np.argsort(counts[positive], kind="stable")]
                    for column_idx in order:
                        counts[column_idx] -= 1
                        adjustment += 1
                        if adjustment == 0:
                            break
                    if adjustment == previous_adjustment:
                        raise RuntimeError(
                            f"failed to reduce excess quota for spot {spot_name!r}"
                        )
            if adjustment != 0 or np.any(counts < 0) or int(counts.sum()) != total_cells:
                raise RuntimeError(f"failed to construct an exact nonnegative quota for spot {spot_name!r}")
            quotas[row_idx] = counts

        return pd.DataFrame(
            quotas,
            index=cell_contributions.index.copy(),
            columns=cell_contributions.columns.copy(),
        )

    def assign_cell_types(self, SVC_obs, spot_cell_distribution):
        """Maximize raw PM score subject to each spot's exact cell-type quota."""
        assigned, quotas, type_list = self._validate_assignment_inputs(
            SVC_obs,
            spot_cell_distribution,
        )
        pm_on_cell = self._validate_and_align_pm(assigned, type_list)
        self.pm_on_cell = pm_on_cell
        assigned["cell_type"] = "Unknown"
        cell_type_column = assigned.columns.get_loc("cell_type")

        for spot_name in quotas.index:
            positions = np.flatnonzero(assigned["spot_name"].to_numpy() == spot_name)
            cell_ids = assigned.iloc[positions]["cell_id"].tolist()
            counts = quotas.loc[spot_name].to_numpy(dtype=np.int64, copy=False)
            slots = np.repeat(np.asarray(type_list, dtype=object), counts)
            scores = pm_on_cell.loc[cell_ids, type_list].to_numpy(dtype=np.float64, copy=False)
            slot_columns = np.repeat(np.arange(len(type_list)), counts)
            row_indices, column_indices = linear_sum_assignment(
                scores[:, slot_columns],
                maximize=True,
            )
            labels = np.empty(len(positions), dtype=object)
            labels[row_indices] = slots[column_indices]
            assigned.iloc[positions, cell_type_column] = labels

        return assigned

    def assign_cell_types_random(self, SVC_obs, spot_cell_distribution):
        """Randomly place exact quota slots using one deterministic generator."""
        assigned, quotas, type_list = self._validate_assignment_inputs(
            SVC_obs,
            spot_cell_distribution,
        )
        assigned["cell_type"] = "Unknown"
        cell_type_column = assigned.columns.get_loc("cell_type")
        rng = np.random.default_rng(int(self.config.sr_assignment_seed))

        for spot_name in quotas.index:
            positions = np.flatnonzero(assigned["spot_name"].to_numpy() == spot_name)
            counts = quotas.loc[spot_name].to_numpy(dtype=np.int64, copy=False)
            slots = np.repeat(np.asarray(type_list, dtype=object), counts)
            assigned.iloc[positions, cell_type_column] = rng.permutation(slots)

        return assigned

    @staticmethod
    def _validate_assignment_inputs(SVC_obs, spot_cell_distribution):
        if not isinstance(SVC_obs, pd.DataFrame) or not {"spot_name", "cell_id"} <= set(SVC_obs.columns):
            raise ValueError("SVC_obs must contain spot_name and cell_id columns")
        assigned = SVC_obs.copy()
        assigned["cell_id"] = assigned["cell_id"].astype(str)
        if assigned["cell_id"].duplicated().any():
            raise ValueError("SVC_obs cell_id values must be unique after string conversion")
        if not isinstance(spot_cell_distribution, pd.DataFrame):
            raise TypeError("spot_cell_distribution must be a pandas DataFrame")
        if spot_cell_distribution.index.has_duplicates:
            raise ValueError("spot_cell_distribution must have unique spot labels")
        if spot_cell_distribution.columns.has_duplicates:
            raise ValueError("spot_cell_distribution must have unique category labels")
        if set(spot_cell_distribution.index) != set(pd.unique(assigned["spot_name"])):
            raise ValueError("spot_cell_distribution spot set must exactly match SVC_obs")

        normalized_types = [str(value).replace("/", "_") for value in spot_cell_distribution.columns]
        if pd.Index(normalized_types).has_duplicates:
            raise ValueError("contribution classes collide after slash normalization")
        try:
            quota_values = spot_cell_distribution.to_numpy(dtype=np.float64, copy=True)
        except (TypeError, ValueError) as exc:
            raise ValueError("spot_cell_distribution values must be numeric") from exc
        if not np.isfinite(quota_values).all():
            raise ValueError("spot_cell_distribution values must be finite")
        if np.any(quota_values < 0) or np.any(quota_values != np.round(quota_values)):
            raise ValueError("spot_cell_distribution must contain nonnegative integer quotas")

        quotas = pd.DataFrame(
            quota_values.astype(np.int64),
            index=spot_cell_distribution.index.copy(),
            columns=normalized_types,
        )
        spot_counts = assigned.groupby("spot_name", sort=False)["cell_id"].size()
        if any(int(quotas.loc[spot].sum()) != int(spot_counts.loc[spot]) for spot in quotas.index):
            raise ValueError("each spot quota must equal its SVC cell count")
        return assigned, quotas, normalized_types

    def _validate_and_align_pm(self, SVC_obs, type_list):
        if not isinstance(self.pm_on_cell, pd.DataFrame):
            raise TypeError("pm_on_cell must be a pandas DataFrame")
        pm_on_cell = self.pm_on_cell.copy()
        if pm_on_cell.index.has_duplicates:
            raise ValueError("pm_on_cell has a duplicate row label")
        if pm_on_cell.columns.has_duplicates:
            raise ValueError("pm_on_cell has a duplicate column label")

        string_index = pm_on_cell.index.astype(str)
        if string_index.has_duplicates:
            raise ValueError("pm_on_cell row labels collide after string conversion")
        normalized_columns = [str(value).replace("/", "_") for value in pm_on_cell.columns]
        if pd.Index(normalized_columns).has_duplicates:
            raise ValueError("pm_on_cell columns collide after slash normalization")
        pm_on_cell.index = string_index
        pm_on_cell.columns = normalized_columns

        svc_cell_ids = SVC_obs["cell_id"].tolist()
        if set(pm_on_cell.index) != set(svc_cell_ids):
            missing = sorted(set(svc_cell_ids) - set(pm_on_cell.index))
            extra = sorted(set(pm_on_cell.index) - set(svc_cell_ids))
            raise ValueError(f"pm_on_cell cell IDs must match exactly; missing={missing}, extra={extra}")
        if set(pm_on_cell.columns) != set(type_list):
            missing = sorted(set(type_list) - set(pm_on_cell.columns))
            extra = sorted(set(pm_on_cell.columns) - set(type_list))
            raise ValueError(f"pm_on_cell classes must match exactly; missing={missing}, extra={extra}")
        try:
            values = pm_on_cell.to_numpy(dtype=np.float64, copy=True)
        except (TypeError, ValueError) as exc:
            raise ValueError("pm_on_cell values must be convertible to float64") from exc
        if not np.isfinite(values).all():
            raise ValueError("pm_on_cell values must be finite")

        numeric = pd.DataFrame(values, index=pm_on_cell.index, columns=pm_on_cell.columns)
        return numeric.loc[svc_cell_ids, type_list]
