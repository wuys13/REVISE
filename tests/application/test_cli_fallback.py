from __future__ import annotations

from pathlib import Path

import pytest


def test_existing_application_path_wins_over_packaged_template(tmp_path):
    from revise.application.cli import _resolve_application_source

    config = tmp_path / "Xenium_T.yaml"
    config.write_text("custom: true\n", encoding="utf-8")

    source = _resolve_application_source(str(config))

    assert source.path == config.resolve()
    assert source.label == str(config.resolve())
    assert source.payload == b"custom: true\n"


def test_official_bare_name_falls_back_to_cwd_mirror(tmp_path, monkeypatch):
    from revise.application.cli import _resolve_application_source

    mirror = tmp_path / "configs" / "application" / "Xenium_T.yaml"
    mirror.parent.mkdir(parents=True)
    mirror.write_bytes(b"mirror: true\n")
    monkeypatch.chdir(tmp_path)

    source = _resolve_application_source("Xenium_T.yaml")

    assert source.path == mirror.resolve()
    assert source.payload == b"mirror: true\n"


def test_official_path_falls_back_to_packaged_template(tmp_path, monkeypatch):
    from revise.application.cli import _resolve_application_source

    monkeypatch.chdir(tmp_path)

    source = _resolve_application_source("configs/application/Xenium_T.yaml")

    assert source.path is None
    assert source.label == "package:revise.application.templates/Xenium_T.yaml"
    assert b"select_cell_type: T" in source.payload


@pytest.mark.parametrize(
    "value",
    [
        "custom/Xenium_T.yaml",
        "Xenium_T.ymll",
        str(Path("/tmp") / "Xenium_T.yaml"),
    ],
)
def test_missing_noncanonical_path_does_not_fallback(value, tmp_path, monkeypatch):
    from revise.application.cli import ApplicationConfigError, _resolve_application_source

    monkeypatch.chdir(tmp_path)

    with pytest.raises(ApplicationConfigError, match="Cannot read application config"):
        _resolve_application_source(value)
