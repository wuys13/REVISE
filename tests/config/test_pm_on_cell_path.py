from revise.config import runner_conf


ApplicationScSrConf = runner_conf.ApplicationScSrConf
BenchmarkSrConf = runner_conf.BenchmarkSrConf


def _base_kwargs(tmp_path):
    return {
        "sample_name": "sample",
        "raw_data_path": str(tmp_path),
        "result_root_path": str(tmp_path / "out"),
        "cell_type_col": "Level1",
        "confidence_col": "confidence",
        "unknown_key": "Unknown",
        "st_file": "spots.h5ad",
        "sc_ref_file": "sc.h5ad",
    }


def test_pm_on_cell_path_is_derived_from_resolved_application_st_path(tmp_path):
    conf = ApplicationScSrConf(**_base_kwargs(tmp_path))

    assert conf.pm_on_cell_file == str(tmp_path / "sample_spots_PM_on_cell.csv")
    assert callable(getattr(runner_conf, "pm_on_cell_path_from_st_path", None))
    assert conf.pm_on_cell_file == runner_conf.pm_on_cell_path_from_st_path(
        conf.st_file_path
    )


def test_pm_on_cell_path_keeps_benchmark_spot_sizes_distinct(tmp_path):
    kwargs = _base_kwargs(tmp_path)
    kwargs["gt_svc_file"] = "gt.h5ad"
    kwargs["annotate_mode"] = "pot"
    spot_25 = BenchmarkSrConf(**kwargs, spot_size=25)
    spot_50 = BenchmarkSrConf(**kwargs, spot_size=50)

    assert spot_25.pm_on_cell_file == str(
        tmp_path / "sample" / "spot_25" / "spots_PM_on_cell.csv"
    )
    assert spot_50.pm_on_cell_file == str(
        tmp_path / "sample" / "spot_50" / "spots_PM_on_cell.csv"
    )
    assert spot_25.pm_on_cell_file != spot_50.pm_on_cell_file
