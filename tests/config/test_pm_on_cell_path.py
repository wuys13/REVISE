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

    assert callable(getattr(runner_conf, "pm_on_cell_path_from_st_path", None))
    assert runner_conf.pm_on_cell_path_from_st_path(conf.st_file_path) == str(
        tmp_path / "sample_spots_PM_on_cell.csv"
    )


def test_pm_on_cell_path_keeps_benchmark_spot_sizes_distinct(tmp_path):
    kwargs = _base_kwargs(tmp_path)
    kwargs["gt_svc_file"] = "gt.h5ad"
    kwargs["annotate_mode"] = "pot"
    spot_25 = BenchmarkSrConf(**kwargs, spot_size=25)
    spot_50 = BenchmarkSrConf(**kwargs, spot_size=50)

    spot_25_pm = runner_conf.pm_on_cell_path_from_st_path(spot_25.st_file_path)
    spot_50_pm = runner_conf.pm_on_cell_path_from_st_path(spot_50.st_file_path)

    assert spot_25_pm == str(
        tmp_path / "sample" / "spot_25" / "spots_PM_on_cell.csv"
    )
    assert spot_50_pm == str(
        tmp_path / "sample" / "spot_50" / "spots_PM_on_cell.csv"
    )
    assert spot_25_pm != spot_50_pm
