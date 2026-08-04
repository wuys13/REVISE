from revise.config import runner_conf


def test_pm_on_cell_path_uses_fixed_data_root_filename(tmp_path):
    assert callable(getattr(runner_conf, "pm_on_cell_path_from_data_root", None))
    assert runner_conf.pm_on_cell_path_from_data_root(str(tmp_path)) == str(
        tmp_path / "PM_on_cell.csv"
    )


def test_pm_on_cell_path_does_not_depend_on_st_filename_or_spot_size(tmp_path):
    data_root = tmp_path / "sample"

    assert runner_conf.pm_on_cell_path_from_data_root(str(data_root)) == str(
        data_root / "PM_on_cell.csv"
    )
