from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace

import pytest


class _Artifact:
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name


def test_run_application_exposes_load_preprocess_and_reconstruct_flow(monkeypatch, tmp_path):
    import reconstruct

    config = SimpleNamespace()
    spatial = object()
    reference = object()
    result = object()
    calls = []

    monkeypatch.setattr(reconstruct, "load_application_yaml", lambda path: (path, {}))
    monkeypatch.setattr(reconstruct, "compile_application_config", lambda doc, source: config)
    monkeypatch.setattr(reconstruct, "load_data", lambda value: calls.append(("load", value)) or (spatial, reference))
    monkeypatch.setattr(
        reconstruct,
        "preprocess_data",
        lambda left, right, value: calls.append(("preprocess", left, right, value)) or (left, right),
    )
    monkeypatch.setattr(
        reconstruct,
        "reconstruct",
        lambda left, right, value: calls.append(("reconstruct", left, right, value)) or result,
    )

    assert reconstruct.run_application(tmp_path / "run.yaml") is result
    assert [call[0] for call in calls] == ["load", "preprocess", "reconstruct"]


def test_main_keeps_reconstruction_output_out_of_success_stdout(monkeypatch):
    import reconstruct

    def run_application(_path):
        print("pipeline noise")
        return _Artifact("AnnData object")

    monkeypatch.setattr(reconstruct, "run_application", run_application)
    stdout = StringIO()

    with redirect_stdout(stdout):
        reconstruct.main(["--config", "run.yaml"])

    assert stdout.getvalue().splitlines() == ["Finished", "AnnData object"]


def test_main_prints_both_sc_svc_result_summaries(monkeypatch):
    import reconstruct

    monkeypatch.setattr(
        reconstruct,
        "run_application",
        lambda _path: (_Artifact("spatial AnnData"), _Artifact("expression AnnData")),
    )
    stdout = StringIO()

    with redirect_stdout(stdout):
        reconstruct.main(["--config", "run.yaml"])

    assert stdout.getvalue().splitlines() == [
        "Finished",
        "spatial: spatial AnnData",
        "expression: expression AnnData",
    ]


def test_main_does_not_print_finished_when_reconstruction_fails(monkeypatch):
    import reconstruct

    def fail(_path):
        raise RuntimeError("publication failed")

    monkeypatch.setattr(reconstruct, "run_application", fail)
    stdout = StringIO()

    with pytest.raises(RuntimeError, match="publication failed"), redirect_stdout(stdout):
        reconstruct.main(["--config", "run.yaml"])

    assert "Finished" not in stdout.getvalue()


def test_sc_svc_sr_preprocessing_ensures_spot_cells_before_generic_steps(monkeypatch):
    import reconstruct

    spatial = object()
    reference = object()
    calls = []
    config = SimpleNamespace(
        svc_type="sc-SVC-sr",
        reference_filter_column="Patient",
        reference_filter_value="P2CRC",
        spatial_min_transcript_counts=60,
        spatial_min_cell_counts=100,
        spatial_min_counts=None,
        reference_min_transcript_counts=None,
        reference_min_cell_counts=100,
        reference_min_genes=None,
        broad_column="Level1",
        subtype_column=None,
    )
    monkeypatch.setattr(
        reconstruct,
        "ensure_all_cells_in_spot",
        lambda value: calls.append(("ensure", value)),
    )
    monkeypatch.setattr(
        reconstruct,
        "filter_reference",
        lambda value, column, selected: calls.append(
            ("filter", value, column, selected)
        )
        or value,
    )
    monkeypatch.setattr(
        reconstruct,
        "preprocess_spatial",
        lambda value, *args, **kwargs: calls.append(("spatial", value)) or value,
    )
    monkeypatch.setattr(
        reconstruct,
        "preprocess_reference",
        lambda value, *args, **kwargs: calls.append(("reference", value)) or value,
    )
    monkeypatch.setattr(
        reconstruct,
        "normalize_reference_labels",
        lambda value, columns: calls.append(("labels", value, tuple(columns))) or value,
    )

    assert reconstruct.preprocess_data(spatial, reference, config) == (
        spatial,
        reference,
    )
    assert [call[0] for call in calls] == [
        "ensure",
        "filter",
        "spatial",
        "reference",
        "labels",
    ]
