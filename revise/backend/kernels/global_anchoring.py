from __future__ import annotations

from anndata import AnnData

from revise.backend.kernels.base import BaseKernel
from revise.backend.kernels.ot import OTKernel
from revise.config.runner_conf import ApplicationScConf


class GlobalAnchoringKernel(BaseKernel):
    """Global-anchoring facade over the shared annotation OT contract."""

    def __init__(self, config, logger):
        super().__init__(config, logger)
        self.mode = self.config.annotate_mode
        self.cell_type_col = self.config.cell_type_col
        self.confidence_col = self.config.confidence_col

    def _tacco_annotate_kwargs(self) -> dict:
        if not isinstance(self.config, ApplicationScConf):
            return {}
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
        return {"multi_center": multi_center, "lamb": lamb}

    def run(self, st_adata: AnnData, sc_ref_adata: AnnData, **kwargs) -> AnnData:
        cell_type_col = kwargs.get("cell_type_col", self.cell_type_col)
        normalized_mode = str(self.mode).strip().lower()
        tacco_kwargs = (
            self._tacco_annotate_kwargs() if normalized_mode == "tacco" else {}
        )
        return OTKernel.annotate(
            st_adata,
            sc_ref_adata,
            method=normalized_mode,
            annotation_key=cell_type_col,
            confidence_key=self.confidence_col,
            pot_reg=kwargs.get("annotate_pot_reg"),
            pot_reg_m=kwargs.get("annotate_pot_reg_m"),
            pot_reg_type=kwargs.get("annotate_pot_reg_type"),
            pot_verbose=True,
            pot_num_iter_max=5000,
            multi_center=tacco_kwargs.get("multi_center"),
            lamb=tacco_kwargs.get("lamb"),
            unknown_key=kwargs.get(
                "unknown_key",
                getattr(self.config, "unknown_key", "Unknown"),
            ),
        )
