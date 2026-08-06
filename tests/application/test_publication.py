import json
import logging
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
import yaml

from revise.recon.context import PipelineContext
from revise.svc import SVC
from revise.utils.io import build_run_dir


def _adata():
    return AnnData(
        X=np.ones((1, 1)),
        obs=pd.DataFrame(index=["cell-1"]),
        var=pd.DataFrame(index=["g1"]),
    )


class _OutputLifecycleContext(SimpleNamespace):
    """Small lifecycle double for direct-write output tests."""

    def register_output_failure_cleanup(self, cleanup):
        self._output_failure_cleanup = cleanup

    def clear_output_failure_cleanup(self):
        self._output_failure_cleanup = None

    def cleanup_failed_outputs(self):
        cleanup = getattr(self, "_output_failure_cleanup", None)
        self._output_failure_cleanup = None
        if cleanup is not None:
            cleanup()

    def record_artifact(self, artifact):
        self.artifact_records.append(artifact)


def _config(tmp_path, svc_type):
    return SimpleNamespace(
        svc_type=svc_type,
        output_dir=tmp_path / "out",
        output_name="sample_result",
        select_cell_type="T" if svc_type == "sc-SVC" else None,
    )


def _ctx(tmp_path, outputs, *, run_dir=None):
    run_dir = Path(run_dir or (tmp_path / "run"))
    run_dir.mkdir(parents=True, exist_ok=True)
    return _OutputLifecycleContext(
        svc=SimpleNamespace(
            artifacts={"outputs": outputs},
            provenance={},
        ),
        profile="application_test",
        run_dir=run_dir,
        merged_config={"ot": {"ga": {"solver": "pot"}, "lr": {"solver": "pot"}}},
        provenance={},
        artifact_records=[],
        logger=logging.getLogger("test-direct-output-persistence"),
    )


def _pipeline_context(tmp_path):
    return PipelineContext(
        merged_config={
            "io": {"save_outputs": False},
            "ot": {
                "ga": {"solver": "pot"},
                "lr": {"solver": "pot"},
            },
        },
        raw_config={},
        config_path="revise/revise.yaml",
        profile="application_test",
        runtime={"mode": "application", "svc_kind": "sp"},
        route_key="application:sp-SVC",
        run_dir=tmp_path,
        logger=logging.getLogger("application-output-lifecycle-test"),
    )


def test_single_output_is_written_directly_under_unique_run_dir(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sp-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sp-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    adata = _adata()
    written = reconstruct._write_outputs(
        config,
        paths,
        _ctx(tmp_path, {"sp_svc": adata}, run_dir=run_dir),
    )

    assert paths["svc"] == run_dir / "sample_result.h5ad"
    assert paths["svc"].parent == run_dir
    assert paths["svc"].is_file()
    assert read_h5ad(paths["svc"]).uns["revise_reconstruction"]["output_role"] == "svc"
    assert read_h5ad(paths["svc"]).uns["revise_reconstruction"]["run_manifest"] == str(
        run_dir / "provenance.json"
    )
    assert written["svc"].shape == adata.shape
    assert not list(config.output_dir.glob("*.h5ad"))
    assert set(run_dir.glob("*.h5ad")) == set(paths.values())
    assert not list(run_dir.glob("*.backup*"))
    assert not list(run_dir.glob("*.tmp*"))
    assert not any(path.is_dir() for path in run_dir.iterdir())
    assert not any(path.name.startswith(".") for path in run_dir.iterdir())


def test_sc_pair_writes_both_final_artifacts_as_run_dir_children(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sc-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sc-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    spatial = _adata()
    expression = _adata()
    reconstruct._write_outputs(
        config,
        paths,
        _ctx(
            tmp_path,
            {"sc_svc_spatial": spatial, "sc_svc_expr": expression},
            run_dir=run_dir,
        ),
    )

    assert all(path.is_file() and path.parent == run_dir for path in paths.values())
    assert paths["spatial"].name.endswith("_spatial.h5ad")
    assert paths["expression"].name.endswith("_expr.h5ad")
    assert not list(config.output_dir.glob("*.h5ad"))
    assert set(run_dir.glob("*.h5ad")) == set(paths.values())
    assert not list(run_dir.glob("*.backup*"))
    assert not list(run_dir.glob("*.tmp*"))
    assert not any(path.is_dir() for path in run_dir.iterdir())
    assert not any(path.name.startswith(".") for path in run_dir.iterdir())


def test_sc_sr_writes_one_logically_named_artifact_under_run_dir(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sc-SVC-sr")
    run_dir = (
        tmp_path
        / "out"
        / "sample_result"
        / "application__sc-SVC-sr"
        / "run-1"
    )
    paths = reconstruct._output_paths(config, run_dir)

    written = reconstruct._write_outputs(
        config,
        paths,
        _ctx(tmp_path, {"sc_svc_dec": _adata()}, run_dir=run_dir),
    )

    assert paths == {"svc": run_dir / "sample_result.h5ad"}
    assert paths["svc"].is_file()
    assert written["svc"].shape == (1, 1)
    assert not list(config.output_dir.glob("*.h5ad"))
    assert set(run_dir.glob("*.h5ad")) == set(paths.values())


def test_repeated_runs_with_same_output_name_keep_each_run_isolated(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sp-SVC")
    run_dirs = [
        build_run_dir(
            str(config.output_dir),
            config.output_name,
            "application:sp-SVC",
        ),
        build_run_dir(
            str(config.output_dir),
            config.output_name,
            "application:sp-SVC",
        ),
    ]
    values = (1.0, 2.0)
    output_paths = []
    for run_dir, value in zip(run_dirs, values):
        adata = _adata()
        adata.X[:] = value
        paths = reconstruct._output_paths(config, run_dir)
        reconstruct._write_outputs(
            config,
            paths,
            _ctx(tmp_path, {"sp_svc": adata}, run_dir=run_dir),
        )
        output_paths.append(paths["svc"])

    assert run_dirs[0] != run_dirs[1]
    assert run_dirs[0].parent == run_dirs[1].parent
    assert run_dirs[0].parent == config.output_dir / config.output_name / "application__sp-SVC"
    assert re.fullmatch(r"\d{8}_\d{6}_\d{6}_[0-9a-f]{8}", run_dirs[0].name)
    assert re.fullmatch(r"\d{8}_\d{6}_\d{6}_[0-9a-f]{8}", run_dirs[1].name)
    assert all(path.parent == run_dir for path, run_dir in zip(output_paths, run_dirs))
    assert [float(read_h5ad(path).X[0, 0]) for path in output_paths] == [1.0, 2.0]
    assert not list(config.output_dir.glob("*.h5ad"))
    assert {path for path in run_dirs[0].parent.iterdir()} == set(run_dirs)


def test_formal_run_application_uses_one_fresh_run_dir_per_invocation(
    tmp_path,
    monkeypatch,
):
    import reconstruct
    import revise.framework as framework

    st = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["spot-1", "spot-2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    sc_ref = AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"Level1": ["A", "B"]}, index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["g1"]),
    )
    st.write_h5ad(tmp_path / "st.h5ad")
    sc_ref.write_h5ad(tmp_path / "sc.h5ad")
    config_path = tmp_path / "application.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "application": {"svc_type": "sp-SVC"},
                "paths": {"root_dir": str(tmp_path)},
                "algorithm": {},
                "inputs": {
                    "st": {"path": "st.h5ad", "format": "h5ad"},
                    "reference": {"path": "sc.h5ad", "format": "h5ad"},
                },
                "global_anchoring": {"broad_column": "Level1"},
                "local_refinement": {"strength": 0.2},
                "output": {"dir": "out", "name": "sample_result"},
                "execution": {"seed": 42},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    class Strategy:
        strategy_id = "Strategy"

        def prepare_context(self, ctx):
            return None

        def global_anchoring(self, ctx):
            return None

        def prepare_local_units(self, ctx):
            return None

        def build_graph(self, ctx):
            return None

        def build_ot_problem(self, ctx):
            return None

        def solve_ot(self, ctx):
            return None

        def update_expression(self, ctx):
            return None

        def finalize_svc(self, ctx):
            output = _adata()
            return SVC(
                expr=output,
                spatial=output,
                svc_kind="sp",
                artifacts={"outputs": {"sp_svc": output}},
            )

    class Registry:
        def get(self, strategy_id):
            return Strategy()

    monkeypatch.setattr(framework, "build_default_registry", Registry)

    executions = [
        reconstruct.run_application(config_path),
        reconstruct.run_application(config_path),
    ]

    run_dirs = [execution.output_paths["svc"].parent for execution in executions]
    assert run_dirs[0] != run_dirs[1]
    assert run_dirs[0].parent == run_dirs[1].parent
    assert run_dirs[0].parent == (
        tmp_path / "out" / "sample_result" / "application__sp-SVC"
    )
    for execution, run_dir in zip(executions, run_dirs):
        assert execution.status == "succeeded"
        assert set(run_dir.iterdir()) == {
            run_dir / "merged_config.json",
            run_dir / "preflight.json",
            run_dir / "provenance.json",
            run_dir / "run.log",
            run_dir / "sample_result.h5ad",
        }
        manifest = json.loads((run_dir / "provenance.json").read_text())
        assert manifest["run"]["status"] == "succeeded"
        assert manifest["results"]["svc"]["path"] == str(
            run_dir / "sample_result.h5ad"
        )
        assert not (run_dir / "artifacts").exists()
        assert not (run_dir / ".revise-run.lock").exists()


def test_pair_second_file_write_failure_cleans_all_direct_outputs(tmp_path, monkeypatch):
    import reconstruct

    config = _config(tmp_path, "sc-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sc-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    spatial = _adata()
    expression = _adata()
    original_write = AnnData.write_h5ad
    calls = 0

    def fail_second_write(self, path, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected expression write failure")
        return original_write(self, path, *args, **kwargs)

    monkeypatch.setattr(AnnData, "write_h5ad", fail_second_write)
    with pytest.raises(OSError, match="injected expression write failure"):
        reconstruct._write_outputs(
            config,
            paths,
            _ctx(
                tmp_path,
                {"sc_svc_spatial": spatial, "sc_svc_expr": expression},
                run_dir=run_dir,
            ),
        )

    assert calls == 2
    assert not list(run_dir.glob("*.h5ad"))
    assert not list(run_dir.glob(".*"))
    assert not any(path.is_dir() for path in run_dir.iterdir())


def test_pair_write_failure_is_not_masked_when_cleanup_is_refused(
    tmp_path,
    monkeypatch,
):
    import reconstruct

    config = _config(tmp_path, "sc-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sc-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    ctx = _ctx(
        tmp_path,
        {"sc_svc_spatial": _adata(), "sc_svc_expr": _adata()},
        run_dir=run_dir,
    )
    original_write = AnnData.write_h5ad
    original_unlink = Path.unlink
    calls = 0

    def fail_second_write(self, path, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected expression write failure")
        return original_write(self, path, *args, **kwargs)

    spatial_target = paths["spatial"].resolve()

    def refuse_spatial_cleanup(self, *args, **kwargs):
        if self.resolve() == spatial_target:
            raise PermissionError("injected cleanup refusal")
        return original_unlink(self, *args, **kwargs)

    monkeypatch.setattr(AnnData, "write_h5ad", fail_second_write)
    monkeypatch.setattr(Path, "unlink", refuse_spatial_cleanup)

    with pytest.raises(OSError, match="injected expression write failure"):
        reconstruct._write_outputs(config, paths, ctx)

    assert spatial_target.is_file()
    assert ctx.provenance["output_cleanup_errors"] == [
        {
            "path": str(spatial_target),
            "type": "PermissionError",
            "message": "injected cleanup refusal",
        }
    ]


def test_post_write_manifest_failure_rolls_back_run_outputs(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sc-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sc-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    ctx = _ctx(
        tmp_path,
        {"sc_svc_spatial": _adata(), "sc_svc_expr": _adata()},
        run_dir=run_dir,
    )
    returned_outputs = {}

    reconstruct._persist_run_outputs(config, paths, ctx, returned_outputs)
    assert all(path.is_file() for path in paths.values())
    assert set(ctx.provenance["results"]) == {"spatial", "expression"}
    assert len(ctx.artifact_records) == 2
    metric_path = run_dir / "metrics" / "synthetic.csv"
    metric_path.parent.mkdir()
    metric_path.write_text("metric,value\nMSE,0\n", encoding="utf-8")
    metric_artifact = {
        "role": "metric:synthetic",
        "path": str(metric_path),
    }
    ctx.artifact_records.append(metric_artifact)

    # The framework invokes this cleanup when a later evaluation or final
    # manifest write fails. The run envelope itself remains inspectable.
    ctx.cleanup_failed_outputs()

    assert not any(path.exists() for path in paths.values())
    assert not list(run_dir.glob("*.h5ad"))
    assert not list(run_dir.glob(".*"))
    assert returned_outputs == {}
    assert "results" not in ctx.provenance
    assert ctx.artifact_records == [metric_artifact]
    assert metric_path.is_file()
    assert run_dir.is_dir()


def test_success_keeps_direct_outputs_after_failure_cleanup_is_cleared(tmp_path):
    import reconstruct

    config = _config(tmp_path, "sp-SVC")
    run_dir = tmp_path / "out" / "sample_result" / "application__sp-SVC" / "run-1"
    paths = reconstruct._output_paths(config, run_dir)
    ctx = _ctx(tmp_path, {"sp_svc": _adata()}, run_dir=run_dir)
    returned_outputs = {}

    reconstruct._persist_run_outputs(config, paths, ctx, returned_outputs)
    ctx.clear_output_failure_cleanup()

    assert paths["svc"].is_file()
    assert returned_outputs["svc"].shape == (1, 1)
    assert ctx.artifact_records[0]["path"] == str(paths["svc"])


def test_application_writer_does_not_stage_or_replace_h5ad_files():
    import reconstruct
    import revise.recon.pipeline as pipeline

    sources = [
        Path(module.__file__).read_text(encoding="utf-8")
        for module in (reconstruct, pipeline)
    ]

    assert all("NamedTemporaryFile" not in source for source in sources)
    assert all("os.replace" not in source for source in sources)


@pytest.mark.parametrize("interruption", [False, True])
def test_succeeded_manifest_failure_restores_output_cleanup_hook(
    tmp_path,
    interruption,
):
    ctx = _pipeline_context(tmp_path)
    output = tmp_path / "sample_result.h5ad"
    output.write_bytes(b"direct output")
    cleanup_calls = []

    def cleanup():
        cleanup_calls.append(True)
        output.unlink(missing_ok=True)

    ctx.register_output_failure_cleanup(cleanup)
    ctx.skip_pending_stages("test")

    def fail_manifest(current):
        assert current.run_status == "succeeded"
        if interruption:
            raise KeyboardInterrupt("interrupted while writing succeeded manifest")
        raise OSError("succeeded manifest write failed")

    ctx.set_provenance_callback(fail_manifest, notify=False)
    with pytest.raises(
        KeyboardInterrupt if interruption else OSError,
        match="interrupted while writing succeeded manifest"
        if interruption
        else "succeeded manifest write failed",
    ):
        ctx.mark_run_succeeded()

    assert ctx.run_status == "running"
    ctx.cleanup_failed_outputs()
    assert cleanup_calls == [True]
    assert not output.exists()


def test_successful_succeeded_manifest_clears_output_cleanup_hook(tmp_path):
    ctx = _pipeline_context(tmp_path)
    output = tmp_path / "sample_result.h5ad"
    output.write_bytes(b"direct output")
    cleanup_calls = []

    def cleanup():
        cleanup_calls.append(True)
        output.unlink(missing_ok=True)

    ctx.register_output_failure_cleanup(cleanup)
    ctx.skip_pending_stages("test")
    observed = []
    ctx.set_provenance_callback(
        lambda current: observed.append(current.run_status),
        notify=False,
    )

    ctx.mark_run_succeeded()
    ctx.cleanup_failed_outputs()

    assert ctx.run_status == "succeeded"
    assert observed == ["succeeded"]
    assert cleanup_calls == []
    assert output.read_bytes() == b"direct output"
