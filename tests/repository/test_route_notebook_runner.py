from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import runpy
import subprocess

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = runpy.run_path(str(ROOT / "reproduce/case/reconstruction_impact/run_route_notebook.py"))



@pytest.fixture(autouse=True)
def forbid_canonical_writes_from_main(monkeypatch):
    """Patch function globals, not the copy returned by runpy.run_path."""
    original = RUNNER["route_paths"]
    def guarded(*args, **kwargs):
        paths = original(*args, **kwargs)
        assert not paths[1].is_relative_to(ROOT / "output"), "Runner tests must use a temporary output root"
        return paths
    monkeypatch.setitem(RUNNER["main"].__globals__, "route_paths", guarded)


def write_pair(tmp_path, *, execution_count=1, cell_id="calculation", output=None):
    source = {"cells": [{"cell_type": "code", "id": "calculation", "source": ["x = 1"], "execution_count": None, "outputs": []}]}
    executed = {"cells": [{"cell_type": "code", "id": cell_id, "source": ["x = 1"], "execution_count": execution_count, "outputs": output or []}]}
    paths = tmp_path / "source.ipynb", tmp_path / "executed.ipynb"
    for path, notebook in zip(paths, (source, executed)):
        path.write_text(json.dumps(notebook))
    return paths


def test_default_output_uses_schema1_sample_root():
    _, output_root, final = RUNNER["route_paths"]("visiumhd")
    assert output_root == ROOT / "output/reconstruction_impact/P1CRC_VisiumHD"
    assert final.parent == output_root / "notebook"


def test_explicit_output_root_overrides_schema1_output_dir(tmp_path):
    _, output_root, final = RUNNER["route_paths"]("visiumhd", str(tmp_path / "override"))
    assert output_root == (tmp_path / "override").resolve()
    assert final == output_root / "notebook/VisiumHD_sp_SVC_Reconstruction_Impact.ipynb"


@pytest.mark.parametrize("count,cell_id", [(None, "calculation"), (1, "changed")])
def test_verification_rejects_unexecuted_or_reidentified_cells(tmp_path, count, cell_id):
    source, executed = write_pair(tmp_path, execution_count=count, cell_id=cell_id)
    with pytest.raises(RuntimeError):
        RUNNER["verify_executed_copy"](source, executed)


def test_verification_preserves_stderr_for_review(tmp_path):
    stderr = {"output_type": "stream", "name": "stderr", "text": ["UserWarning: low support\n"]}
    paths = write_pair(tmp_path, output=[stderr])
    audit = RUNNER["verify_executed_copy"](*paths)
    assert audit["executed_code_cells"] == 1
    assert audit["stderr"] == [{"cell_id": "calculation", "text": "UserWarning: low support\n"}]


def test_verification_rejects_notebook_error(tmp_path):
    error = {"output_type": "error", "ename": "ValueError", "evalue": "bad input"}
    paths = write_pair(tmp_path, output=[error])
    with pytest.raises(RuntimeError, match="error output"):
        RUNNER["verify_executed_copy"](*paths)


def test_verification_rejects_cell_type_and_position_swap(tmp_path):
    source = tmp_path / "source.ipynb"
    executed = tmp_path / "executed.ipynb"
    source.write_text(json.dumps({
        "cells": [
            {"cell_type": "code", "id": "a", "source": ["x = 1"], "execution_count": None, "outputs": []},
            {"cell_type": "markdown", "id": "b", "source": ["explanation"]},
        ]
    }))
    executed.write_text(json.dumps({
        "cells": [
            {"cell_type": "markdown", "id": "a", "source": ["changed explanation"]},
            {"cell_type": "code", "id": "b", "source": ["x = 1"], "execution_count": 1, "outputs": []},
        ]
    }))

    with pytest.raises(RuntimeError):
        RUNNER["verify_executed_copy"](source, executed)


def test_failed_rerun_cannot_leave_a_successful_execution_audit(tmp_path, monkeypatch):
    source, executed = write_pair(tmp_path)
    audit_path = tmp_path / "execution.json"
    audit_path.write_text('{"status": "completed"}')
    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, tmp_path, executed))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd"])

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, "nbconvert")

    monkeypatch.setattr(main.__globals__["subprocess"], "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        main()
    assert json.loads(audit_path.read_text())["status"] == "failed"


def test_interrupted_run_cannot_leave_a_running_execution_audit(tmp_path, monkeypatch):
    source, executed = write_pair(tmp_path)
    audit_path = tmp_path / "execution.json"
    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, tmp_path, executed))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd"])

    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(main.__globals__["subprocess"], "run", interrupt)
    with pytest.raises(KeyboardInterrupt):
        main()
    audit = json.loads(audit_path.read_text())
    assert audit["status"] == "failed"
    assert "interrupted" in audit["reason"].lower()


def test_artifact_index_failure_cannot_leave_a_completed_execution_audit(tmp_path, monkeypatch):
    source, executed = write_pair(tmp_path)
    output_root = tmp_path / "sample"
    final = output_root / "notebook" / "FormalRoute.ipynb"
    final.parent.mkdir(parents=True)
    final.write_text(executed.read_text(encoding="utf-8"), encoding="utf-8")
    analysis_dir = output_root / "analysis"
    analysis_dir.mkdir(parents=True)
    (analysis_dir / "artifacts.csv").write_text("path,size_bytes,sha256\n", encoding="utf-8")

    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, output_root, final))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd"])
    monkeypatch.setattr(main.__globals__["subprocess"], "run", lambda *args, **kwargs: None)
    monkeypatch.setitem(main.__globals__, "publish_report", lambda *args: {})
    monkeypatch.setitem(main.__globals__, "verify_fresh_evidence", lambda *args: [])

    def fail_index(*args, **kwargs):
        raise RuntimeError("artifact index failed")

    monkeypatch.setitem(main.__globals__, "append_execution_artifacts", fail_index)

    with pytest.raises(RuntimeError, match="artifact index failed"):
        main()

    assert json.loads((final.parent / "execution.json").read_text())["status"] == "failed"


def test_run_start_and_failure_hide_previous_report_conclusions(tmp_path, monkeypatch):
    source, executed = write_pair(tmp_path)
    output_root = tmp_path / "sample"
    final = output_root / "notebook" / "FormalRoute.ipynb"
    final.parent.mkdir(parents=True)
    final.write_text(executed.read_text(encoding="utf-8"), encoding="utf-8")
    old_report = output_root / "report.html"
    old_report.parent.mkdir(parents=True, exist_ok=True)
    old_report.write_text("<p>OLD_CONCLUSION</p>", encoding="utf-8")
    old_records = output_root / "analysis" / "conclusions" / "report_records.json"
    old_records.parent.mkdir(parents=True)
    old_records.write_text('{"conclusion": "OLD_CONCLUSION"}', encoding="utf-8")

    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, output_root, final))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd"])

    def fail_after_start(*args, **kwargs):
        if old_report.is_file():
            assert "OLD_CONCLUSION" not in old_report.read_text(encoding="utf-8")
        raise subprocess.CalledProcessError(1, "nbconvert")

    monkeypatch.setattr(main.__globals__["subprocess"], "run", fail_after_start)
    with pytest.raises(subprocess.CalledProcessError):
        main()

    assert "OLD_CONCLUSION" not in old_report.read_text(encoding="utf-8")
    assert json.loads(old_records.read_text())["review"]["not_current"] is True


def test_successful_run_appends_final_notebook_and_execution_audit_to_index(tmp_path, monkeypatch):
    source = tmp_path / "source.ipynb"
    output_root = tmp_path / "sample"
    final = output_root / "notebook" / "FormalRoute.ipynb"
    final.parent.mkdir(parents=True)
    source_notebook = {
        "cells": [
            {
                "cell_type": "code",
                "id": "calculation",
                "source": ["x = 1"],
                "execution_count": None,
                "outputs": [],
            }
        ]
    }
    executed_notebook = {
        **source_notebook,
        "cells": [
            {
                **source_notebook["cells"][0],
                "execution_count": 1,
            }
        ],
    }
    source.write_text(json.dumps(source_notebook), encoding="utf-8")
    final.write_text(json.dumps(executed_notebook), encoding="utf-8")
    analysis_dir = output_root / "analysis"
    analysis_dir.mkdir(parents=True)
    (analysis_dir / "artifacts.csv").write_text(
        "path,size_bytes,sha256\nanalysis/prior.csv,3,prior-hash\n",
        encoding="utf-8",
    )

    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, output_root, final))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd"])
    monkeypatch.setattr(main.__globals__["subprocess"], "run", lambda *args, **kwargs: None)
    monkeypatch.setitem(main.__globals__, "publish_report", lambda *args: {})
    monkeypatch.setitem(main.__globals__, "verify_fresh_evidence", lambda *args: [])

    main()

    audit_path = final.parent / "execution.json"
    rows = list(csv.DictReader((analysis_dir / "artifacts.csv").open(newline="", encoding="utf-8")))
    assert [row["path"] for row in rows] == [
        "analysis/prior.csv",
        "notebook/FormalRoute.ipynb",
        "notebook/execution.json",
        "report.html",
    ]
    indexed = {row["path"]: row for row in rows}
    for path in (final, audit_path):
        relative = path.relative_to(output_root).as_posix()
        assert indexed[relative]["size_bytes"] == str(path.stat().st_size)
        assert indexed[relative]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert "analysis/artifacts.csv" not in indexed


def test_report_placeholders_only_change_presentation(tmp_path):
    """Attaching saved conclusions must preserve executable code and cell IDs."""
    import importlib.util
    path = Path(__file__).resolve().parents[2] / 'reproduce/case/reconstruction_impact/run_route_notebook.py'
    spec = importlib.util.spec_from_file_location('report_runner_test', path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    assert hasattr(runner, 'attach_conclusions')
    notebook = {'cells': [
        {'id': 'overview', 'cell_type': 'markdown', 'metadata': {'impact_records': 'overview'}, 'source': ['placeholder']},
        {'id': 'science', 'cell_type': 'code', 'source': ['x = 1'], 'execution_count': 1, 'outputs': []},
    ]}
    target = tmp_path / 'executed.ipynb'
    target.write_text(json.dumps(notebook))
    runner.attach_conclusions(target, {}, formatter=lambda package, layer=None: '<p>Saved evidence</p>')
    updated = json.loads(target.read_text())
    assert updated['cells'][1] == notebook['cells'][1]
    assert updated['cells'][0]['id'] == 'overview'
    assert 'Saved evidence' in ''.join(updated['cells'][0]['source'])


def test_stale_figure_cannot_be_registered_as_current_evidence(tmp_path):
    directory = tmp_path / "figures"; directory.mkdir()
    figure = directory / "old.png"; figure.write_bytes(b"old")
    with pytest.raises(RuntimeError, match="old.png"):
        RUNNER["verify_fresh_evidence"](tmp_path, figure.stat().st_mtime_ns + 1)


def test_route_configuration_changes_invalidate_source_identity(tmp_path, monkeypatch):
    function = RUNNER["route_source_digests"]
    monkeypatch.setitem(function.__globals__, "ROOT", tmp_path)
    path = tmp_path / RUNNER["ROUTES"]["visiumhd"][0]
    path.parent.mkdir(parents=True); path.write_text("method: original")
    before = function("visiumhd")
    path.write_text("method: changed")
    assert function("visiumhd") != before


def test_report_source_identity_includes_consumed_template(tmp_path, monkeypatch):
    main = RUNNER['main']
    monkeypatch.setitem(main.__globals__, 'ROOT', tmp_path)
    source = tmp_path / 'reproduce/case/reconstruction_impact'
    source.mkdir(parents=True)
    for name in ('report_records.py', 'report_html.py'):
        (source / name).write_text('# reader\n')
    template = tmp_path / 'docs/design/reconstruction-impact/templates/report.html'
    template.parent.mkdir(parents=True)
    template.write_text('<main>$body</main>')
    identity = main.__globals__['report_source_digests']()
    relative = template.relative_to(tmp_path).as_posix()
    assert identity[relative] == hashlib.sha256(template.read_bytes()).hexdigest()
    template.write_text('<main class="evidence">$body</main>')
    assert main.__globals__['report_source_digests']() != identity


def test_executed_notebook_method_links_follow_output_location(tmp_path, monkeypatch):
    attach = RUNNER['attach_conclusions']
    monkeypatch.setitem(attach.__globals__, 'ROOT', tmp_path)
    final = tmp_path / 'output/reconstruction_impact/sample/notebook/result.ipynb'
    final.parent.mkdir(parents=True)
    notebook = {'cells': [
        {'id': 'method', 'cell_type': 'markdown', 'source': '[Method](../../../docs/design/reconstruction-impact/gene-and-function.md#moran)'},
        {'id': 'code', 'cell_type': 'code', 'source': 'x = 1'},
    ]}
    final.write_text(json.dumps(notebook))
    attach(final, {}, formatter=lambda *a, **k: '')
    result = json.loads(final.read_text())
    href = result['cells'][0]['source'].split('](')[1].split(')')[0]
    assert (final.parent / href.split('#')[0]).resolve() == tmp_path / 'docs/design/reconstruction-impact/gene-and-function.md'
    assert href.endswith('#moran')
    assert result['cells'][1] == notebook['cells'][1]
    attach(final, {}, formatter=lambda *a, **k: '')
    assert json.loads(final.read_text()) == result


def test_baseline_change_blocks_publication_before_publish_report(tmp_path, monkeypatch):
    source, executed = write_pair(tmp_path)
    output_root = tmp_path / "P1CRC_VisiumHD"
    final = output_root / "notebook" / "FormalRoute.ipynb"
    final.parent.mkdir(parents=True)
    final.write_text(executed.read_text(encoding="utf-8"), encoding="utf-8")
    table_path = output_root / "analysis" / "scientific.csv"
    table_path.parent.mkdir(parents=True)
    pd.DataFrame({"value": [2]}).to_csv(table_path, index=False)
    baseline = {
        "P1CRC_VisiumHD": {
            "tables": {
                "analysis/scientific.csv": {
                    "columns": ["value"],
                    "rows": 1,
                    "digest": hashlib.sha256(
                        pd.util.hash_pandas_object(
                            pd.DataFrame({"value": [1]}), index=True
                        ).values.tobytes()
                    ).hexdigest(),
                }
            }
        }
    }
    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text(json.dumps(baseline), encoding="utf-8")

    main = RUNNER["main"]
    monkeypatch.setitem(main.__globals__, "ROOT", tmp_path)
    monkeypatch.setitem(main.__globals__, "route_paths", lambda *args: (source, output_root, final))
    monkeypatch.setattr(main.__globals__["sys"], "argv", ["runner", "visiumhd", "--baseline", str(baseline_path)])
    monkeypatch.setattr(main.__globals__["subprocess"], "run", lambda *args, **kwargs: None)
    monkeypatch.setitem(main.__globals__, "local_source_digests", lambda: {})
    monkeypatch.setitem(main.__globals__, "route_source_digests", lambda *args: {})
    monkeypatch.setitem(main.__globals__, "verify_fresh_evidence", lambda *args: [])
    monkeypatch.setitem(main.__globals__, "append_execution_artifacts", lambda *args: None)
    published = []
    monkeypatch.setitem(main.__globals__, "publish_report", lambda *args: published.append(True))

    with pytest.raises(RuntimeError, match="baseline"):
        main()

    assert not published
    assert json.loads((final.parent / "execution.json").read_text())["status"] == "failed"
