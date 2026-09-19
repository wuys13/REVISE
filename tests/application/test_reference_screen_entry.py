"""Independent screen command uses the real host callback and emits consumable YAML."""
import json
import shutil

import yaml

from tests.application.test_reference_ga import _write_config, _write_inputs


def test_screen_command_real_ga_produces_verified_reference(tmp_path):
    from revise.reference_preparation.__main__ import main
    from revise.reference_preparation.evidence import resolve_reference_input

    _, reference = _write_inputs(tmp_path)
    application = _write_config(tmp_path, svc_type="sp-SVC", mode=None)
    candidates = tmp_path / "candidates"
    candidates.mkdir()
    for name in ("a", "b"):
        shutil.copyfile(reference, candidates / f"{name}.h5ad")
    config = tmp_path / "prepare.yaml"
    config.write_text(yaml.safe_dump({
        "schema_version": 1, "mode": "screen",
        "reconstruction_config": str(application),
        "candidates": {"directory": "candidates"}, "output_dir": "prepared",
    }))
    assert main(["--config", str(config)]) == 0
    selected, evidence = resolve_reference_input(tmp_path / "prepared" / "reference.yaml")
    assert selected == candidates / "a.h5ad"
    assert evidence["verification"] == "verified_report"
    report = json.loads((tmp_path / "prepared" / "report.json").read_text())
    assert [row["status"] for row in report["candidates"]] == ["success", "success"]
    assert report["candidates"][1]["comparability"]["st_unit_ids"] == "equal"
    assert not list(tmp_path.rglob("SVC.h5ad"))
