from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import revise.analysis as analysis


ROOT = Path(__file__).resolve().parents[2]


def test_sc_svc_analysis_service_is_removed_from_public_api():
    assert "ScSVCAnalysisService" not in analysis.__all__
    with pytest.raises(AttributeError):
        getattr(analysis, "ScSVCAnalysisService")


def test_sp_svc_analysis_service_remains_compatible():
    adata = AnnData(
        X=np.ones((4, 1)),
        obs=pd.DataFrame(
            {
                "pred": ["a", "a", "b", "b"],
                "truth": ["x", "x", "y", "y"],
            },
            index=["c0", "c1", "c2", "c3"],
        ),
    )
    service = analysis.SpSVCAnalysisService(SimpleNamespace(spatial=adata))

    assert service.clustering_metrics("pred", "truth") == {"ari": 1.0, "nmi": 1.0}


def test_api_docs_do_not_reference_removed_service():
    index = (ROOT / "docs/source/api/index.rst").read_text(encoding="utf-8")
    generated = (
        ROOT
        / "docs/source/api/generated/revise.analysis.ScSVCAnalysisService.rst"
    )

    assert "ScSVCAnalysisService" not in index
    assert not generated.exists()
