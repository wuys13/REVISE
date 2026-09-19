from dataclasses import replace

import pytest

from revise.application.config import ApplicationConfigError, compile_application_config, load_application_yaml
from revise.application.publication import application_metadata, output_paths
from tests.application.test_request import _document, _write_config


@pytest.mark.parametrize("old", [None, {"format": "invalid", "filter_column": "obsolete"},
                                 {"path": "/missing/old.h5ad", "format": "h5ad"}])
def test_override_precedes_old_reference_validation(tmp_path, old):
    document = _document()
    document["inputs"]["reference"] = old
    source, document = load_application_yaml(_write_config(tmp_path, document))
    reference = tmp_path / "selected.h5ad"
    reference.write_bytes(b"input validated by loader later")
    config = compile_application_config(document, source=source, reference_override=reference)
    assert config.reference_path == reference.resolve()
    assert config.reference_filter_column is config.reference_filter_value is None
    assert document["inputs"]["reference"] == old
    with pytest.raises(ApplicationConfigError):
        compile_application_config(document, source=source)


@pytest.mark.parametrize("name", ["missing.h5ad", "wrong.csv"])
def test_invalid_override_never_falls_back(tmp_path, name):
    source, document = load_application_yaml(_write_config(tmp_path, _document()))
    candidate = tmp_path / name
    if name.endswith("csv"):
        candidate.write_text("wrong format")
    with pytest.raises(ApplicationConfigError, match="existing H5AD"):
        compile_application_config(document, source=source, reference_override=candidate)


def test_override_metadata_uses_actual_reference_and_evidence(tmp_path):
    source, document = load_application_yaml(_write_config(tmp_path, _document()))
    reference = tmp_path / "selected.h5ad"
    reference.write_bytes(b"reference")
    config = compile_application_config(document, source=source, reference_override=reference)
    evidence = {"origin": "manual", "reference": {"path": str(reference), "sha256": "sample"}}
    config = replace(config, reference_preparation=evidence)
    metadata = application_metadata(config, paths=output_paths(config))
    assert metadata["resolved_inputs"]["reference"] == str(reference)
    assert metadata["source_sha256"] == config.config_sha256
    assert metadata["effective_request"]["inputs"]["reference_preparation"] == evidence


def test_cli_passes_external_reference_config(monkeypatch):
    import reconstruct

    captured = {}
    monkeypatch.setattr(reconstruct, "run_application", lambda path, **kwargs: captured.update(path=path, **kwargs))
    reconstruct.main(["--config", "run.yaml", "--reference-config", "prepared/reference.yaml"])
    assert captured == {"path": "run.yaml", "select_ct": None,
                        "reference_config": "prepared/reference.yaml"}


def _prepared_reference(tmp_path):
    import numpy as np
    import pandas as pd
    import yaml
    from anndata import AnnData
    from revise.reference_preparation import prepare_reference

    source = tmp_path / "pool.h5ad"
    AnnData(np.ones((2, 2)), obs=pd.DataFrame({"pair": ["a", "b"]}, index=["r1", "r2"]),
            var=pd.DataFrame(index=["g1", "g2"])).write_h5ad(source)
    preparation = tmp_path / "preparation.yaml"
    preparation.write_text(yaml.safe_dump({"schema_version": 1, "mode": "paired",
        "source": "pool.h5ad", "pair_column": "pair", "pair_key": "a", "output_dir": "prepared"}))
    return prepare_reference(preparation)


def test_run_application_consumes_paired_result_without_old_filter(monkeypatch, tmp_path):
    import reconstruct

    monkeypatch.setattr(reconstruct, "bind_expression_sources", lambda config, *_: config)

    prepared = _prepared_reference(tmp_path)
    document = _document()
    document["inputs"]["reference"] = {"filter_column": "nonexistent"}
    application = _write_config(tmp_path, document)
    captured = {}

    def load(config):
        captured["config"] = config
        return object(), object()

    monkeypatch.setattr(reconstruct, "load_data", load)
    monkeypatch.setattr(reconstruct, "is_sample_delivery", lambda config: False)
    monkeypatch.setattr(reconstruct, "preprocess_data", lambda st, ref, config: (st, ref))
    result = object()
    monkeypatch.setattr(reconstruct, "reconstruct", lambda *args: result)
    assert reconstruct.run_application(application, reference_config=prepared.reference_config_path) is result
    config = captured["config"]
    assert config.reference_path == prepared.reference_path
    assert config.reference_filter_column is config.reference_filter_value is None
    assert config.reference_preparation["verification"] == "verified_report"


def test_reference_change_during_load_prevents_reconstruction(monkeypatch, tmp_path):
    import reconstruct

    prepared = _prepared_reference(tmp_path)
    application = _write_config(tmp_path, _document())

    def load(config):
        config.reference_path.write_bytes(b"changed")
        return object(), object()

    monkeypatch.setattr(reconstruct, "load_data", load)
    monkeypatch.setattr(reconstruct, "is_sample_delivery", lambda config: False)
    monkeypatch.setattr(reconstruct, "reconstruct", lambda *args: pytest.fail("must not reconstruct"))
    with pytest.raises(ValueError, match="Reference file SHA-256"):
        reconstruct.run_application(application, reference_config=prepared.reference_config_path)


def test_reference_change_during_compute_preserves_previous_delivery(monkeypatch, tmp_path):
    import reconstruct
    from revise.reference_preparation.evidence import resolve_reference_input
    from tests.application.test_sample_delivery import delivery_fixture

    prepared = _prepared_reference(tmp_path)
    selected, evidence = resolve_reference_input(prepared.reference_config_path)
    config, raw, ctx = delivery_fixture(tmp_path)
    config = replace(config, reference_path=selected, reference_preparation=evidence)
    paths = output_paths(config)
    config.output_dir.mkdir(parents=True)
    previous = {name: ("previous-" + name).encode() for name in paths}
    for name, path in paths.items():
        path.write_bytes(previous[name])

    def run(_self, **kwargs):
        selected.write_bytes(b"changed during calculation")
        kwargs["finalize_callback"](ctx)

    monkeypatch.setattr(reconstruct.REVISEPipeline, "run", run)
    with pytest.raises(ValueError, match="Reference file SHA-256"):
        reconstruct.reconstruct(raw.copy(), raw.copy(), config, raw_adata=raw)
    assert {name: path.read_bytes() for name, path in paths.items()} == previous


def test_publication_cannot_overwrite_external_reference_config(tmp_path):
    from revise.application.publication import publish_outputs
    from revise.reference_preparation.evidence import resolve_reference_input
    from tests.application.test_sample_delivery import delivery_fixture

    prepared = _prepared_reference(tmp_path)
    selected, evidence = resolve_reference_input(prepared.reference_config_path)
    config, raw, ctx = delivery_fixture(tmp_path)
    config = replace(config, reference_path=selected, reference_preparation=evidence)
    paths = output_paths(config)
    paths['sample_config'] = prepared.reference_config_path
    previous = prepared.reference_config_path.read_bytes()
    with pytest.raises(ValueError, match="aliases reference preparation input"):
        publish_outputs(config, paths, ctx, raw=raw)
    assert prepared.reference_config_path.read_bytes() == previous
