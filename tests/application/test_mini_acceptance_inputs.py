"""The acceptance subset keeps a contiguous original coordinate/ID scope."""
from pathlib import Path
import runpy

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
HELPERS = runpy.run_path(str(ROOT / "scripts/prepare_mini_acceptance.py"))


def test_square_roi_preserves_order_and_all_boundary_ties():
    xy = np.asarray([[2, 0], [0, 0], [-2, 0], [0, 2], [0, -2], [10, 10], [-10, -10]])
    before = xy.copy()
    positions, selection = HELPERS["square_roi"](xy, 0.2)
    assert positions.tolist() == [0, 1, 2, 3, 4]
    assert selection["target_n_obs"] == 2
    assert selection["selected_n_obs"] == 5
    np.testing.assert_array_equal(xy, before)


@pytest.mark.parametrize("fraction", [0, -0.1, 1.1, float("nan")])
def test_roi_rejects_invalid_fraction(fraction):
    with pytest.raises(ValueError, match="fraction"):
        HELPERS["square_roi"](np.zeros((3, 2)), fraction)


def test_sst_request_uses_real_p2_counts_without_mouse_pm():
    doc = HELPERS["request"](ROOT, "P2CRC_Visium", "mini.h5ad", "output/test", "P2CRC_Visium_mini", "default")
    assert doc["inputs"]["reference"]["filter_value"] == "P2CRC"
    assert doc["inputs"]["st"]["expression"]["identity"] == "raw_counts"
    assert "pm_on_cell" not in doc["inputs"]
    assert doc["delivery"]["coordinates"]["microns_per_coordinate"] == 0.73
    assert doc["algorithm"]["sr_cell_count_method"] == "cyto_linear_v1"


def test_hst_request_does_not_invent_expression_history_or_change_reference_scope():
    doc = HELPERS["request"](ROOT, "P1CRC_HD", "mini.h5ad", "output/test", "P1CRC_HD_mini", "default")
    assert doc["inputs"]["st"]["expression"]["identity"] == "unknown"
    assert "filter_value" not in doc["inputs"]["reference"]
    assert doc["algorithm"]["ot_method"] == "pot"


def test_preparation_never_retargets_existing_requests(tmp_path):
    configs = tmp_path / "configs"
    existing = configs / "mini/P2CRC_Xenium-random.yaml"
    existing.parent.mkdir(parents=True)
    existing.write_text("original request")
    output = tmp_path / "new-run"
    with pytest.raises(FileExistsError, match="Existing request is immutable"):
        HELPERS["prepare"](tmp_path, output, configs, 0.01)
    assert existing.read_text() == "original request"
    assert not output.exists()
