"""Analysis code dependency and loaded-source identity contracts."""
import importlib

import pytest

from test_runner import fake_solver as _fake_solver


fake_solver = _fake_solver


def test_declared_helper_change_is_included_in_analysis_identity(tmp_path, monkeypatch):
    from revise.batch import analysis

    helper = tmp_path / "impact_helper.py"
    helper.write_text("VALUE = 'before'\n")
    adapter = tmp_path / "impact_adapter.py"
    adapter.write_text(
        "CODE_DEPENDENCIES = ('impact_helper',)\n"
        "def run(context):\n"
        "    return {'artifacts': {'x': {'path': 'x', 'description': 'x'}},\n"
        "            'calculation': {'input_view': 'x', 'parameters': {},\n"
        "                           'comparison_basis': 'x'}}\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    module = importlib.import_module("impact_adapter")
    dependencies = analysis._adapter_code_dependencies(module)

    before = analysis._analysis_code_identity(adapter, dependencies)
    helper.write_text("VALUE = 'after'\n")
    after = analysis._analysis_code_identity(adapter, dependencies)

    assert before[str(helper)]["sha256"] != after[str(helper)]["sha256"]


def test_loaded_helper_change_fails_before_cached_adapter_can_publish(
    tmp_path, monkeypatch, fake_solver
):
    import yaml

    from revise.batch.analysis import run_analysis_task
    from revise.batch.runner import run_batch
    from test_runner import batch_config, make_sample, result_root

    helper = tmp_path / "impact_runtime_helper.py"
    helper.write_text("VALUE = 'before'\n")
    adapter = tmp_path / "impact_runtime_adapter.py"
    adapter.write_text(
        "import impact_runtime_helper\n"
        "CODE_DEPENDENCIES = ('impact_runtime_helper',)\n"
        "def run(context):\n"
        "    (context.output_dir / 'value.txt').write_text(impact_runtime_helper.VALUE)\n"
        "    return {'artifacts': {'value': {'path': 'value.txt', 'description': 'value'}},\n"
        "            'calculation': {'input_view': 'native_carriers', 'parameters': {},\n"
        "                           'comparison_basis': 'helper fixture'}}\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    sample = tmp_path / "data" / "one"
    make_sample(sample, cell_types=["T"])
    config = batch_config(tmp_path)
    document = yaml.safe_load(config.read_text())
    document["analysis"] = {
        "impact": {
            "entrypoint": "impact_runtime_adapter:run",
            "version": "1",
        }
    }
    config.write_text(yaml.safe_dump(document))
    assert run_batch(config)["summary"]["succeeded"] == 1

    first = run_analysis_task(config, "one", "impact", cell_type="T")
    assert first["status"] == "succeeded"
    assert (result_root(sample) / "T" / "analysis" / "impact" / "value.txt").read_text() == "before"

    helper.write_text("VALUE = 'after'\n")
    second = run_analysis_task(config, "one", "impact", cell_type="T")

    assert second["status"] == "failed"
    assert "restart" in second["error"].lower()
    assert (result_root(sample) / "T" / "analysis" / "impact" / "value.txt").read_text() == "before"


def test_first_analysis_reloads_a_preloaded_helper_constant(tmp_path, monkeypatch, fake_solver):
    import os
    import yaml

    from revise.batch.analysis import run_analysis_task
    from revise.batch.runner import run_batch
    from test_runner import batch_config, make_sample, result_root

    helper = tmp_path / "preloaded_constant_helper.py"
    helper.write_text("VALUE = 'before'\n")
    adapter = tmp_path / "preloaded_constant_adapter.py"
    adapter.write_text(
        "from preloaded_constant_helper import VALUE\n"
        "CODE_DEPENDENCIES = ('preloaded_constant_helper',)\n"
        "def run(context):\n"
        "    (context.output_dir / 'value.txt').write_text(VALUE)\n"
        "    return {'artifacts': {'value': {'path': 'value.txt', 'description': 'value'}},\n"
        "            'calculation': {'input_view': 'native_carriers', 'parameters': {},\n"
        "                           'comparison_basis': 'helper fixture'}}\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.import_module("preloaded_constant_adapter")
    original_times = helper.stat()
    helper.write_text("VALUE = 'after!'\n")
    os.utime(helper, ns=(original_times.st_atime_ns, original_times.st_mtime_ns))

    sample = tmp_path / "data" / "one"
    make_sample(sample, cell_types=["T"])
    config = batch_config(tmp_path)
    document = yaml.safe_load(config.read_text())
    document["analysis"] = {
        "impact": {
            "entrypoint": "preloaded_constant_adapter:run",
            "version": "1",
        }
    }
    config.write_text(yaml.safe_dump(document))
    assert run_batch(config)["summary"]["succeeded"] == 1

    result = run_analysis_task(config, "one", "impact", cell_type="T")

    assert result["status"] == "succeeded"
    assert (result_root(sample) / "T" / "analysis" / "impact" / "value.txt").read_text() == "after!"


def test_changed_declared_dependencies_do_not_attest_the_adapter(tmp_path, monkeypatch):
    from revise.batch import analysis

    first = tmp_path / "attestation_first.py"
    first.write_text("VALUE = 'first'\n")
    second = tmp_path / "attestation_second.py"
    second.write_text("VALUE = 'second'\n")
    adapter = tmp_path / "attestation_adapter.py"
    adapter.write_text("CODE_DEPENDENCIES = ('attestation_first',)\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.import_module("attestation_adapter")
    adapter.write_text("CODE_DEPENDENCIES = ('attestation_second',)\n")

    with pytest.raises(RuntimeError, match="declared dependencies changed"):
        analysis._attest_adapter_module("attestation_adapter")

    assert str(adapter) not in analysis._ADAPTER_SOURCES


def test_declared_dependency_cycles_are_rejected(tmp_path, monkeypatch):
    from revise.batch import analysis

    (tmp_path / "cycle_first.py").write_text("import cycle_second\n")
    (tmp_path / "cycle_second.py").write_text("import cycle_first\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    first = importlib.import_module("cycle_first")
    second = importlib.import_module("cycle_second")

    with pytest.raises(ValueError, match="cycle"):
        analysis._declared_dependency_order((first, second))
