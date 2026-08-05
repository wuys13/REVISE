from pathlib import Path
from types import SimpleNamespace


def test_run_pipeline_uses_the_application_selector(monkeypatch, tmp_path):
    from revise.application import service

    captured = {}

    class Pipeline:
        def __init__(self, config_path=None):
            captured["config_path"] = config_path

        def _execute_run(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace(provenance={})

    monkeypatch.setattr(service, "REVISEPipeline", Pipeline)
    monkeypatch.setattr(service, "_build_io_overrides", lambda _request: {})
    monkeypatch.setattr(service, "_build_algorithm_overrides", lambda _request: {})
    request = SimpleNamespace(
        svc_type="sp-SVC",
        seed=17,
        effective_action="preflight",
        source_path="application.yaml",
        config_sha256="a" * 64,
        declared_root=".",
        resolved_root=tmp_path,
        cwd=tmp_path,
        resolved_paths={},
        action="preflight",
        dry_run_override=False,
    )

    service._run_pipeline(request)

    assert captured["config_path"] is None
    assert captured["svc_type"] == "sp-SVC"
    assert captured["cf"] is None
    assert captured["runtime_overrides"] == {"seed": 17}


def test_execute_application_returns_engine_preflight_evidence(monkeypatch, tmp_path):
    from revise.application import service

    svc = SimpleNamespace(
        provenance={
            "run_dir": str(tmp_path / "run"),
            "profile": "application_sp",
            "route": {
                "mode": "application",
                "application_route": "sp-SVC",
                "task": "sp_svc",
                "strategy": "SpSvcApplicationStrategy",
            },
        },
        summary=lambda: {"status": "ready"},
    )
    monkeypatch.setattr(service, "_run_pipeline", lambda request: svc)

    execution = service.execute_application(
        SimpleNamespace(svc_type="sp-SVC", effective_action="preflight")
    )

    assert execution.status == "ready"
    assert execution.preflight == Path(tmp_path / "run" / "preflight.json")
    assert execution.pipeline["profile"] == "application_sp"
    assert execution.pipeline["task"] == "sp_svc"
    assert execution.pipeline["strategy"] == "SpSvcApplicationStrategy"
    assert "confounding" not in execution.pipeline["route"]
