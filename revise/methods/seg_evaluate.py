import logging
import os

import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData

from revise.methods.base_method import BaseMethod


class SegEvaluate(BaseMethod):
    """
    Segmentation error evaluation and cell classification.
    
    This class evaluates segmentation errors in spatial transcriptomics data
    and classifies cells as needing correction (no_effect=False) or not
    (no_effect=True) based on confidence scores and total counts.
    """
    def __init__(self, config, logger):
        super().__init__(config, logger)
        self.plot_conf = {
            "palettes": {
                'Diminishing': '#8ECFC9',
                'Expanding': '#FFBE7A',
                'Unchanged': '#82B0D2'
            },
            "bins": 20
        }
        
    def run(self, adata: AnnData, logger: logging.Logger):
        """
        Evaluate segmentation errors and classify cells.
        
        Args:
            adata: AnnData object with segmentation error labels in obs['seg_error']
                and confidence scores in obs['Confidence']
            logger: Logger instance for logging results
            
        Returns:
            AnnData: Updated AnnData with 'no_effect' column indicating
                whether each cell needs correction (False) or can remain
                unchanged (True)
                
        This method:
        1. Calculates recovery rates for different error types
        2. Generates visualization plots
        3. Identifies cells that don't need correction based on thresholds
        """
        adata = adata.copy()
        ## Reset cells that are Diminishing and Expanding
        dropout_total_counts = self.config.dropout_total_counts
        swapping_total_counts = self.config.swapping_total_counts
        lower_ts = self.config.lower_ts
        upper_ts = self.config.upper_ts

        adata.obs[['x', 'y']] = adata.obsm['spatial']
        adata.obs['total_counts'] = adata.X.sum(axis=1).A1
        b = adata.obs.copy()
        c = b[b["total_counts"] > swapping_total_counts]
        d = b[b["total_counts"] < swapping_total_counts]
        # Calculate Diminishing Ratio
        diminishing_condition_1 = (b['total_counts'] < dropout_total_counts) & (b['seg_error'] == "Diminishing")
        diminishing_condition_2 = (b['total_counts'] > dropout_total_counts) & (b['Confidence'] < lower_ts) & (
                    b['seg_error'] == "Diminishing")
        diminishing_ratio = (diminishing_condition_1.sum() + diminishing_condition_2.sum()) / (
                    b['seg_error'] == "Diminishing").sum()

        # Calculate Expanding Ratio
        expanding_condition = (b['total_counts'] > swapping_total_counts) & (b['Confidence'] < upper_ts) & (
                    b['seg_error'] == "Expanding")
        expanding_ratio = expanding_condition.sum() / (b['seg_error'] == "Expanding").sum()

        # Calculate Unchanged Ratio
        unchanged_condition_1 = (b['total_counts'] > swapping_total_counts) & (b['Confidence'] > upper_ts) & (
                    b['seg_error'] == "Unchanged")
        unchanged_condition_2 = (b['total_counts'] < swapping_total_counts) & (b['Confidence'] > lower_ts) & (
                    b['seg_error'] == "Unchanged")
        unchanged_ratio = (unchanged_condition_1.sum() + unchanged_condition_2.sum()) / (
                    b['seg_error'] == "Unchanged").sum()

        # Print results
        logger.info("Final recovery rates:")
        logger.info(f"Diminishing Ratio: {diminishing_ratio}")
        logger.info(f"Expanding Ratio: {expanding_ratio}")
        logger.info(f"Unchanged Ratio: {unchanged_ratio}")

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))


        sns.histplot(data=b, x='Confidence', hue='seg_error', bins=self.plot_conf["bins"], palette=self.plot_conf["palettes"], kde=True, ax=ax1)
        sns.histplot(data=c, x='Confidence', hue='seg_error', bins=self.plot_conf["bins"], palette=self.plot_conf["palettes"], kde=True, ax=ax2)
        sns.histplot(data=d, x='Confidence', hue='seg_error', bins=self.plot_conf["bins"], palette=self.plot_conf["palettes"], kde=True, ax=ax3)

        ax1.set_title('All Cells')
        ax2.set_title(f'Cells with Total Counts > {swapping_total_counts}')
        ax3.set_title(f'Cells with Total Counts < {swapping_total_counts}')

        plt.tight_layout()
        fig_file = os.path.join(self.config.result_dir, "segmentation.pdf")
        plt.savefig(fig_file)
        logger.info(f"save segmentation plot in {fig_file}")

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        sns.histplot(data=b, x='total_counts', hue='seg_error', bins=self.plot_conf["bins"], palette=self.plot_conf["palettes"], kde=True, ax=ax)
        plt.tight_layout()
        plt.savefig(os.path.join(self.config.result_dir, "total_counts.pdf"))

        no_effect_indices = ((b['total_counts'] > swapping_total_counts) & (b['Confidence'] > upper_ts)) | (
                    (b['total_counts'] < swapping_total_counts) & (b['Confidence'] > lower_ts))
        logger.info(f"no_effect_indices ratio : {no_effect_indices.sum()} / {len(b)}")
        adata.obs['no_effect'] = no_effect_indices
        return adata
