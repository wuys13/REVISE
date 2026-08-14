from __future__ import annotations

from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.kernels.ot import OTKernel


class LocalAnchoringKernel(BaseKernel):
    """Annotate one sc-SVC local unit through the configured local OT solver."""

    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.method = str(self.config.rec_ot_method).strip().lower()
        self.confidence_col = self.config.confidence_col
        self.unknown_key = self.config.unknown_key

    def run(self, target: AnnData, reference: AnnData, **kwargs) -> AnnData:
        cell_type_col = kwargs.get("cell_type_col", self.config.cell_type_col)
        multi_center = None
        lamb = None
        if self.method == "tacco":
            multi_center = self.config.tacco_annotate_multi_center
            lamb = self.config.tacco_annotate_lamb
            missing = [
                name
                for name, value in (("multi_center", multi_center), ("lamb", lamb))
                if value is None
            ]
            if missing:
                raise ValueError(
                    "Application sc-SVC TACCO annotation parameters are missing: "
                    f"{missing}. Resolve them from sc.tacco_annotate."
                )
        return OTKernel.annotate(
            target,
            reference,
            method=self.method,
            annotation_key=cell_type_col,
            confidence_key=self.confidence_col,
            pot_reg=self.config.rec_pot_reg,
            pot_reg_m=self.config.rec_pot_reg_m,
            pot_reg_type=self.config.rec_pot_reg_type,
            pot_verbose=False,
            pot_num_iter_max=5000,
            multi_center=multi_center,
            lamb=lamb,
            unknown_key=self.unknown_key,
        )
