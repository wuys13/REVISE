import os

import numpy as np
import pandas as pd
from tqdm import tqdm

from revise.methods.base_method import BaseMethod


class SpotSr(BaseMethod):
    """
    Spot-level cell type allocation method for virtual cells.
    
    This class assigns cell types to each virtual cell within spots based on
    cell type probability matrices and spot-level cell type distributions.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)
        if os.path.exists(self.config.pm_on_cell_file):
            self.pm_on_cell = pd.read_csv(self.config.pm_on_cell_file, index_col=0)
        else:
            self.pm_on_cell = None

    def run(self, sc_svc):
        """
        Assign cell types to virtual cells based on spot-level profiles.
        
        Args:
            sc_svc: ScSVCSr instance containing:
                - st_adata: Spatial data with cell type contributions in obsm
                - svc_obs: DataFrame with spot_name and cell_id columns
                
        This method:
        1. Calculates optimal cell type distribution for each spot
        2. Assigns cell types to cells based on probabilities and distributions
        3. Logs matching statistics between assigned and true cell types
        """
        sc_svc.svc_obs = self.assign_cell_types_easy(sc_svc.svc_obs, sc_svc.st_adata.obsm["Level1"], "max")
        sc_svc.svc_obs['match'] = sc_svc.svc_obs['cell_type'] == sc_svc.svc_obs['true_cell_type']
        max_match = sum(sc_svc.svc_obs['match'] == True) / len(sc_svc.svc_obs)
        self.logger.info(f"{sc_svc.svc_obs['match'].value_counts()} \n {max_match}")

        spot_cell_distribution = self.get_spot_cell_distribution(
            cell_contributions=sc_svc.st_adata.obsm["Level1"], SVC_obs=sc_svc.svc_obs
        )
        if self.pm_on_cell is not None:
            sc_svc.svc_obs = self.assign_cell_types(
                SVC_obs=sc_svc.svc_obs, spot_cell_distribution=spot_cell_distribution)
        else:
            self.logger.warning(f"No pm_on_cell file: {self.pm_on_cell}, random assign!")
            sc_svc.svc_obs = self.assign_cell_types_easy(
                SVC_obs=sc_svc.svc_obs, cell_contributions=sc_svc.st_adata.obsm["Level1"], mode="random"
            )
        sc_svc.svc_obs['match'] = sc_svc.svc_obs['cell_type'] == sc_svc.svc_obs['true_cell_type']
        enhance_match = sum(sc_svc.svc_obs['match'] == True) / len(sc_svc.svc_obs)
        self.logger.info(f"{sc_svc.svc_obs['match'].value_counts()} \n {enhance_match}")

    def assign_cell_types(self, SVC_obs, spot_cell_distribution):
        """
        Args:
            SVC_obs: DataFrame containing 'spot_name' and 'cell_id' columns
            spot_cell_distribution: DataFrame, number of cells of each cell type in each spot

        Returns:
            SVC_obs: Updated SVC_obs containing 'cell_type' column
        """
        SVC_obs = SVC_obs.copy()
        type_list = list(spot_cell_distribution.columns)
        self.pm_on_cell = self.pm_on_cell.loc[SVC_obs['cell_id'].values, type_list]
        self.logger.info(f"pm_on_cell shape: {self.pm_on_cell.shape}")

        if 'cell_type' not in SVC_obs.columns:
            SVC_obs['cell_type'] = "Unknown"

        spot_groups = SVC_obs.groupby('spot_name')

        for spot_name in tqdm(spot_cell_distribution.index, desc="Assigning cell types"):
            spot_cells_df = spot_groups.get_group(spot_name)
            valid_cells = spot_cells_df[spot_cells_df['cell_id'].isin(self.pm_on_cell.index)]
            if len(valid_cells) == 0:
                self.logger.warn(f"Warning: No valid cells found for spot {spot_name}")
                continue

            cell_probs = self.pm_on_cell.loc[valid_cells['cell_id']].values
            target_counts = spot_cell_distribution.loc[spot_name].astype(int)

            cell_type_indices = np.argmax(cell_probs, axis=1)
            initial_types = np.array(type_list)[cell_type_indices]

            type_counts = pd.Series(0, index=type_list)
            for t in initial_types:
                type_counts[t] += 1

            adjustments = []
            for cell_type in type_list:
                target = int(target_counts[cell_type])
                current = type_counts[cell_type]
                if current != target:
                    adjustments.append({
                        'cell_type': cell_type,
                        'difference': current - target,
                        'target': target
                    })

            adjustments.sort(key=lambda x: abs(x['difference']), reverse=True)

            for adj in adjustments:
                cell_type = adj['cell_type']
                difference = adj['difference']

                if difference > 0:
                    mask = initial_types == cell_type
                    cells_of_type = np.where(mask)[0]

                    if len(cells_of_type) > 0:
                        probs = cell_probs[cells_of_type]
                        probs[:, type_list.index(cell_type)] = -np.inf

                        best_alternative_scores = np.max(probs, axis=1)
                        cells_to_change = cells_of_type[np.argsort(best_alternative_scores)[-difference:]]

                        for cell_idx in cells_to_change:
                            new_type_idx = np.argmax(cell_probs[cell_idx])
                            initial_types[cell_idx] = type_list[new_type_idx]

                elif difference < 0:
                    mask = initial_types != cell_type
                    other_cells = np.where(mask)[0]

                    if len(other_cells) > 0:
                        probs = cell_probs[other_cells]
                        type_idx = type_list.index(cell_type)

                        best_cells = other_cells[np.argsort(probs[:, type_idx])[-abs(difference):]]
                        initial_types[best_cells] = cell_type

            SVC_obs.loc[valid_cells.index, 'cell_type'] = initial_types

        return SVC_obs

    def get_spot_cell_distribution(self, cell_contributions, SVC_obs):
        """
        Optimize cell type proportions for each spot to match actual cell counts.
        
        This function converts continuous cell type contribution proportions
        to integer cell counts that sum to the actual number of cells in each spot.
        
        Args:
            cell_contributions: DataFrame or array of shape (n_spots, n_cell_types)
                containing cell type contribution proportions for each spot
            SVC_obs: DataFrame with 'spot_name' and 'cell_id' columns
                used to determine actual cell counts per spot
                
        Returns:
            pd.DataFrame: DataFrame of shape (n_spots, n_cell_types) with
                integer cell counts for each cell type in each spot
        """
        spot_cell_counts = SVC_obs.groupby('spot_name')['cell_id'].count().to_dict()
        spot_names = cell_contributions.index

        spot_cell_distribution = []

        for spot_idx, spot_name in enumerate(spot_names):
            total_cells = spot_cell_counts[spot_name]
            original_contrib = cell_contributions.loc[spot_name].values
            cell_counts = original_contrib * total_cells

            if self.config.svc_completeness:
                cell_counts = np.round(cell_counts)

                current_sum = np.sum(cell_counts)
                if current_sum != total_cells:
                    adjustment = total_cells - current_sum
                    if adjustment > 0:
                        max_idx = np.argmax(cell_counts)
                        cell_counts[max_idx] += adjustment
                    elif adjustment < 0:
                        print(spot_name, adjustment, total_cells, current_sum, cell_counts)
                        non_zero_indices = np.where(cell_counts >= 1)[0]
                        sorted_indices = non_zero_indices[np.argsort(cell_counts[non_zero_indices])]

                        for idx in sorted_indices:
                            if adjustment >= 0:
                                break
                            cell_counts[idx] -= 1
                            adjustment += 1
                        print(spot_name, adjustment, total_cells, current_sum, cell_counts)

                spot_cell_distribution.append(cell_counts)

        spot_cell_distribution = pd.DataFrame(np.array(spot_cell_distribution), index=spot_names,
                                              columns=cell_contributions.columns)

        return spot_cell_distribution

    def assign_cell_types_easy(self, SVC_obs, cell_contributions, mode):
        """
        Args:
            SVC_obs: DataFrame，包含'spot_name'和'cell_id'列
            cell_contributions: DataFrame，每个spot中存在的cell类型贡献

        Returns:
            SVC_obs: 更新后的SVC_obs，包含'cell_type'列
        """
        SVC_obs = SVC_obs.copy()
        cell_contributions = cell_contributions.copy()

        assert set(cell_contributions.index) == set(
            SVC_obs['spot_name'].unique()), "cell_contributions的index与SVC_obs中的spot_name不匹配"

        if 'cell_type' not in SVC_obs.columns:
            SVC_obs['cell_type'] = "Unknown"

        spot_groups = SVC_obs.groupby('spot_name')

        for spot_name in tqdm(cell_contributions.index, desc="Assigning cell types"):
            spot_cells_df = spot_groups.get_group(spot_name)
            spot_cells = spot_cells_df['cell_id'].values

            spot_contributions = cell_contributions.loc[spot_name]

            if mode == "max":
                max_type = spot_contributions.idxmax()
                SVC_obs.loc[spot_cells_df.index, 'cell_type'] = max_type
            elif mode == "random":
                if len(spot_cells) == 1:
                    max_type = spot_contributions.idxmax()
                    SVC_obs.loc[spot_cells_df.index, 'cell_type'] = max_type
                else:
                    sorted_types = spot_contributions.sort_values(ascending=False).index
                    if len(sorted_types) >= 2:
                        max_type = sorted_types[0]
                        second_type = sorted_types[1]
                        SVC_obs.loc[spot_cells_df.index, 'cell_type'] = max_type
                        random_index = np.random.choice(spot_cells_df.index)
                        SVC_obs.loc[random_index, 'cell_type'] = second_type
                    else:
                        max_type = sorted_types[0]
                        SVC_obs.loc[spot_cells_df.index, 'cell_type'] = max_type
            else:
                raise ValueError("mode has to be 'random' or 'max'.")

        return SVC_obs
