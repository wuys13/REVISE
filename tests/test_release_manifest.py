from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import jsonschema
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE_ROOT = ROOT / "release" / "0.1.0rc1"
SCHEMA_PATH = RELEASE_ROOT / "release-manifest.schema.json"
TEMPLATE_PATH = RELEASE_ROOT / "release-manifest.template.json"
GATE_SCHEMA_PATH = RELEASE_ROOT / "gate-report.schema.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _populated_manifest(tmp_path: Path) -> dict:
    template = json.loads(TEMPLATE_PATH.read_text())
    wheel = tmp_path / "revise_svc-0.1.0rc1-py3-none-any.whl"
    sdist = tmp_path / "revise_svc-0.1.0rc1.tar.gz"
    constraint_310 = ROOT / "constraints" / "python-3.10.txt"
    constraint_311 = ROOT / "constraints" / "python-3.11.txt"
    wheel.write_bytes(b"wheel fixture")
    sdist.write_bytes(b"sdist fixture")
    gate_ids = [
        "unit",
        "integration",
        "installed_cli",
        "package",
        "docs",
        "repository",
        "tacco",
    ]
    reports = {}
    for gate_id in gate_ids:
        report = tmp_path / f"{gate_id}.json"
        report.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "gate": gate_id,
                    "conclusion": "success",
                    "head_sha": "a" * 40,
                    "run_id": 123,
                    "run_attempt": 1,
                    "run_url": "https://github.com/wuys13/REVISE/actions/runs/123",
                    "tested_python": (
                        ["3.10", "3.11"]
                        if gate_id in {"unit", "installed_cli"}
                        else ["3.11"]
                    ),
                    "runner": "ubuntu-latest",
                    "os": "Linux",
                    "architecture": "X64",
                }
            )
            + "\n"
        )
        reports[gate_id] = report

    manifest = copy.deepcopy(template)
    manifest["release"]["tag"] = "v0.1.0rc1"
    manifest["source"].update({"commit": "a" * 40, "tree": "b" * 40})
    manifest["build"].update(
        {
            "workflow_run": "https://github.com/wuys13/REVISE/actions/runs/123",
            "python": {"implementation": "CPython", "version": "3.11.13"},
            "tools": {"setuptools": "80.9.0", "wheel": "0.45.1"},
            "constraints": [
                {
                    "python": "3.10",
                    "path": "constraints/python-3.10.txt",
                    "sha256": _sha256(constraint_310),
                },
                {
                    "python": "3.11",
                    "path": "constraints/python-3.11.txt",
                    "sha256": _sha256(constraint_311),
                },
            ],
        }
    )
    manifest["artifacts"] = [
        {
            "role": "wheel",
            "filename": wheel.name,
            "size": wheel.stat().st_size,
            "sha256": _sha256(wheel),
        },
        {
            "role": "sdist",
            "filename": sdist.name,
            "size": sdist.stat().st_size,
            "sha256": _sha256(sdist),
        },
    ]
    manifest["ci"] = {
        "head_sha": "a" * 40,
        "gates": [
            {
                "id": gate_id,
                "conclusion": "success",
                "run_id": 123,
                "run_attempt": 1,
                "head_sha": "a" * 40,
                "report": reports[gate_id].name,
                "sha256": _sha256(reports[gate_id]),
            }
            for gate_id in gate_ids
        ],
    }
    manifest["docs"].update(
        {
            "source_commit": "a" * 40,
            "workflow_run": "https://github.com/wuys13/REVISE/actions/runs/123",
        }
    )
    return manifest


def test_tracked_template_is_incomplete_but_populates_into_strict_final_schema(
    tmp_path,
):
    schema = json.loads(SCHEMA_PATH.read_text())
    template = json.loads(TEMPLATE_PATH.read_text())

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(template, schema)

    manifest = _populated_manifest(tmp_path)
    jsonschema.validate(manifest, schema, format_checker=jsonschema.FormatChecker())
    assert "${" not in json.dumps(manifest)


def test_detached_manifest_cross_checks_declared_files(tmp_path):
    from tools.release_manifest import validate_declared_files

    manifest = _populated_manifest(tmp_path)
    validate_declared_files(manifest, source_root=ROOT, assets_root=tmp_path)

    (tmp_path / manifest["artifacts"][0]["filename"]).write_bytes(b"changed")
    with pytest.raises(ValueError, match="artifact (size|hash) mismatch"):
        validate_declared_files(manifest, source_root=ROOT, assets_root=tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("gate", "docs", "content mismatch"),
        ("conclusion", "failure", "content mismatch"),
        ("head_sha", "c" * 40, "content mismatch"),
        ("run_id", 124, "content mismatch"),
        ("run_attempt", 2, "content mismatch"),
        ("run_url", "https://example.test/wrong", "run URL mismatch"),
    ],
)
def test_gate_report_content_cannot_disagree_with_manifest(
    tmp_path, field, value, message
):
    from tools.release_manifest import validate_declared_files

    manifest = _populated_manifest(tmp_path)
    unit_gate = next(gate for gate in manifest["ci"]["gates"] if gate["id"] == "unit")
    report = tmp_path / unit_gate["report"]
    payload = json.loads(report.read_text())
    payload[field] = value
    report.write_text(json.dumps(payload))
    unit_gate["sha256"] = _sha256(report)

    with pytest.raises(ValueError, match=message):
        validate_declared_files(manifest, source_root=ROOT, assets_root=tmp_path)


def test_final_manifest_cannot_declare_itself_as_an_artifact(tmp_path):
    schema = json.loads(SCHEMA_PATH.read_text())
    manifest = _populated_manifest(tmp_path)
    manifest["artifacts"][0]["filename"] = "release-manifest.json"

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(manifest, schema)


def test_artifact_roles_require_unique_matching_filenames(tmp_path):
    from tools.release_manifest import validate_declared_files

    schema = json.loads(SCHEMA_PATH.read_text())
    manifest = _populated_manifest(tmp_path)
    manifest["artifacts"][1]["filename"] = manifest["artifacts"][0]["filename"]
    manifest["artifacts"][1]["size"] = manifest["artifacts"][0]["size"]
    manifest["artifacts"][1]["sha256"] = manifest["artifacts"][0]["sha256"]
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(manifest, schema)
    with pytest.raises(ValueError, match="filenames must be unique"):
        validate_declared_files(manifest, source_root=ROOT, assets_root=tmp_path)


def test_schema_rejects_non_string_pattern_fields_and_duplicate_limitations(tmp_path):
    from tools.release_manifest import validate_declared_files

    schema = json.loads(SCHEMA_PATH.read_text())
    manifest = _populated_manifest(tmp_path)
    manifest["docs"]["version"] = 7
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(manifest, schema)

    manifest = _populated_manifest(tmp_path)
    manifest["limitations"].append(copy.deepcopy(manifest["limitations"][0]))
    with pytest.raises(ValueError, match="limitation IDs must be unique"):
        validate_declared_files(manifest, source_root=ROOT, assets_root=tmp_path)


def test_creator_writes_a_valid_manifest_only_outside_the_source_tree(
    tmp_path, monkeypatch
):
    from tools import create_release_manifest as creator

    expected = _populated_manifest(tmp_path)
    build_report = tmp_path / "build-report.json"
    build_report.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "build": {
                    key: expected["build"][key]
                    for key in ["python", "backend", "tools", "isolation"]
                },
                "artifacts": expected["artifacts"],
                "source": {
                    "commit": expected["source"]["commit"],
                    "tree": expected["source"]["tree"],
                },
            }
        )
    )
    monkeypatch.setattr(
        creator,
        "_git",
        lambda _repo, value: "b" * 40 if value.endswith("^{tree}") else "a" * 40,
    )
    monkeypatch.setattr(creator, "_git_status", lambda _repo: "")
    output = tmp_path / "release-manifest.json"
    args = SimpleNamespace(
        source_root=ROOT,
        template=TEMPLATE_PATH,
        schema=SCHEMA_PATH,
        gate_schema=GATE_SCHEMA_PATH,
        assets_dir=tmp_path,
        build_report=build_report,
        output=output,
        tag="v0.1.0rc1",
        run_id=123,
        run_attempt=1,
        run_url="https://github.com/wuys13/REVISE/actions/runs/123",
        docs_run_url=None,
        gate=[
            (gate_id, tmp_path / f"{gate_id}.json")
            for gate_id in [
                "unit",
                "integration",
                "installed_cli",
                "package",
                "docs",
                "repository",
                "tacco",
            ]
        ],
    )
    manifest = creator.create_manifest(args)
    assert json.loads(output.read_text()) == manifest
    assert manifest["source"] == expected["source"]

    unit_report = tmp_path / "unit.json"
    broken = json.loads(unit_report.read_text())
    broken["head_sha"] = "c" * 40
    unit_report.write_text(json.dumps(broken))
    with pytest.raises(ValueError, match="head SHA mismatch"):
        creator.create_manifest(args)
    unit_report.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "gate": "unit",
                "conclusion": "success",
                "head_sha": "a" * 40,
                "run_id": 123,
                "run_attempt": 1,
                "run_url": "https://github.com/wuys13/REVISE/actions/runs/123",
                "tested_python": ["3.10", "3.11"],
                "runner": "ubuntu-latest",
                "os": "Linux",
                "architecture": "X64",
            }
        )
        + "\n"
    )

    build_payload = json.loads(build_report.read_text())
    build_payload["source"]["tree"] = "c" * 40
    build_report.write_text(json.dumps(build_payload))
    with pytest.raises(ValueError, match="build report source"):
        creator.create_manifest(args)
    build_payload["source"]["tree"] = "b" * 40
    build_report.write_text(json.dumps(build_payload))

    monkeypatch.setattr(
        creator,
        "_git",
        lambda _repo, value: (
            "b" * 40
            if value.endswith("^{tree}")
            else ("c" * 40 if value == "HEAD" else "a" * 40)
        ),
    )
    with pytest.raises(ValueError, match="HEAD does not match"):
        creator.create_manifest(args)
    monkeypatch.setattr(
        creator,
        "_git",
        lambda _repo, value: "b" * 40 if value.endswith("^{tree}") else "a" * 40,
    )
    monkeypatch.setattr(creator, "_git_status", lambda _repo: " M tracked-file\n")
    with pytest.raises(ValueError, match="checkout must be clean"):
        creator.create_manifest(args)
    monkeypatch.setattr(creator, "_git_status", lambda _repo: "")

    args.tag = "v0.1.0rc2"
    with pytest.raises(ValueError, match="tag does not match"):
        creator.create_manifest(args)
    args.tag = "v0.1.0rc1"

    invalid_report = json.loads(unit_report.read_text())
    invalid_report["run_url"] = "not a uri"
    unit_report.write_text(json.dumps(invalid_report))
    with pytest.raises(jsonschema.ValidationError):
        creator.create_manifest(args)
    invalid_report["run_url"] = args.run_url
    unit_report.write_text(json.dumps(invalid_report))

    template_copy = json.loads(TEMPLATE_PATH.read_text())
    modified_template = tmp_path / "template.json"
    template_copy["docs"]["version"] = "0.1.0rc2"
    modified_template.write_text(json.dumps(template_copy))
    args.template = modified_template
    with pytest.raises(ValueError, match="docs version"):
        creator.create_manifest(args)

    template_copy = json.loads(TEMPLATE_PATH.read_text())
    template_copy["support"]["tested_extras"].append("cci")
    modified_template.write_text(json.dumps(template_copy))
    with pytest.raises(ValueError, match="corresponding release gate"):
        creator.create_manifest(args)

    template_copy = json.loads(TEMPLATE_PATH.read_text())
    template_copy["support"]["tested_platforms"].append(
        {"runner": "imaginary", "os": "ImaginaryOS", "architecture": "ARM64"}
    )
    modified_template.write_text(json.dumps(template_copy))
    with pytest.raises(ValueError, match="tested platforms"):
        creator.create_manifest(args)
    args.template = TEMPLATE_PATH

    args.output = tmp_path / "unit.json"
    with pytest.raises(ValueError, match="assets directory"):
        creator.create_manifest(args)

    args.output = output
    original_gates = args.gate
    args.gate = [("unit", output), *original_gates[1:]]
    with pytest.raises(ValueError, match="collides with declared evidence"):
        creator.create_manifest(args)
    args.gate = original_gates

    args.output = ROOT / "release" / "0.1.0rc1" / "release-manifest.json"
    with pytest.raises(ValueError, match="outside the source tree"):
        creator.create_manifest(args)


def test_constraint_matrix_and_workflows_encode_the_release_boundary():
    constraints = {
        path.name: path.read_text()
        for path in sorted((ROOT / "constraints").glob("python-*.txt"))
    }
    assert set(constraints) == {"python-3.10.txt", "python-3.11.txt"}
    for content in constraints.values():
        assert "setuptools==80.9.0" in content
        assert "wheel==0.45.1" in content
        assert "POT==0.9.5" in content
        assert "zarr==2.18.7" in content

    ci = yaml.safe_load((ROOT / ".github/workflows/ci.yml").read_text())
    matrix = ci["jobs"]["test"]["strategy"]["matrix"]["python-version"]
    assert matrix == ["3.10", "3.11"]
    ci_text = (ROOT / ".github/workflows/ci.yml").read_text()
    assert "constraints/python-${{ matrix.python-version }}.txt" in ci_text

    release_text = (ROOT / ".github/workflows/release.yml").read_text()
    lowered = release_text.lower()
    assert ci_text.count("tools/build_release_artifacts.py") == 1
    assert {
        "test",
        "unit-report",
        "distribution",
        "installed-cli",
        "installed-cli-report",
        "package",
        "docs",
        "repository",
        "tacco-smoke",
    } <= set(ci["jobs"])
    assert "--ignore=tests/integration" in ci_text
    assert ci_text.count("tools/create_ci_gate_reports.py") == 7
    for gate_id in [
        "unit",
        "integration",
        "installed_cli",
        "package",
        "docs",
        "repository",
        "tacco",
    ]:
        assert ci_text.count(f"--gate {gate_id}") == 1
        assert f'name: gate-{gate_id.replace("_", "-")}' in ci_text
        assert f'--gate "{gate_id}=${{ASSET_DIR}}/gate-{gate_id}.json"' in release_text
    assert "uses: ./.github/workflows/ci.yml" in release_text
    assert "tools/create_ci_gate_reports.py" not in release_text
    assert "pattern: gate-*" in release_text
    assert "tools/create_release_manifest.py" in release_text
    assert release_text.count("--gate ") == 7
    assert "git commit" not in lowered
    assert "git push" not in lowered
    assert "release-manifest.json" not in (ROOT / "MANIFEST.in").read_text()
    assert not (RELEASE_ROOT / "release-manifest.json").exists()
    tracked = subprocess.run(
        ["git", "ls-files", "release/0.1.0rc1/release-manifest.json"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    assert tracked.stdout == ""
