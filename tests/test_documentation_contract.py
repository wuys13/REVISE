from __future__ import annotations

import json
import re
import runpy
from pathlib import Path
from types import SimpleNamespace

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]
CLAIM_DOCS = (
    ROOT / "README.md",
    ROOT / "docs/index.rst",
    ROOT / "docs/source/installation.rst",
    ROOT / "docs/source/quickstart.rst",
    ROOT / "docs/source/configuration.rst",
    ROOT / "docs/source/concepts.rst",
    ROOT / "docs/source/benchmark.rst",
    ROOT / "docs/source/case.rst",
    ROOT / "docs/source/architecture.rst",
    ROOT / "docs/source/limitations.rst",
)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _joined(paths=CLAIM_DOCS) -> str:
    return "\n".join(_read(path) for path in paths)


def test_public_docs_distinguish_installed_source_and_paper_entry_paths():
    text = _joined()
    metadata = _read(ROOT / "pyproject.toml")

    assert "revise-reconstruct" in text
    assert 'revise-reconstruct = "revise.cli:main"' in metadata
    assert "python reconstruct.py" in text
    assert (ROOT / "reconstruct.py").is_file()
    assert "python benchmark_main.py" in text
    assert (ROOT / "benchmark_main.py").is_file()
    assert (ROOT / "application_sp_SVC_recon.py").is_file()
    assert (ROOT / "application_sc_SVC_recon.py").is_file()
    assert "installed command" in text.lower()
    assert "compatibility wrapper" in text.lower()
    assert "paper reproduction" in text.lower()


def test_quickstart_matches_application_input_resolution_and_route_fields():
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    normalized = " ".join(quickstart.split())
    runner_contract = runpy.run_path(ROOT / "revise/config/runner_conf.py")

    specs = runner_contract["resolve_input_specs"](
        {"mode": "application", "task": "sp_svc"},
        {
            "data_root": "data",
            "sample_name": "sample",
            "st_file": "st.h5ad",
            "sc_ref_file": "sc_ref.h5ad",
        },
    )
    paths = {spec.role: spec.path for spec in specs}
    st_path = paths["st"]
    sc_path = paths["sc_ref"]
    assert st_path == "data/sample_st.h5ad"
    assert sc_path == "data/sc_ref.h5ad"
    assert st_path in quickstart
    assert sc_path in quickstart
    assert "unique ``obs_names`` and unique ``var_names``" in normalized
    assert "hST requires ``Level1``" in normalized
    assert "iST and sST require both ``Level1`` and ``Level2``" in normalized
    assert "default ``Patient`` column" in normalized
    assert "--spot-size" not in quickstart


def test_installation_distinguishes_wheel_and_checkout_extras():
    installation = _read(ROOT / "docs/source/installation.rst")
    metadata = tomllib.loads(_read(ROOT / "pyproject.toml"))
    template = json.loads(
        _read(ROOT / "release/0.1.0rc1/release-manifest.template.json")
    )
    optional = metadata["project"]["optional-dependencies"]

    assert optional["tacco"] == ["tacco==0.5.0"]
    assert set(optional) == {
        "tacco",
        "pathway",
        "cci",
        "trajectory",
        "spatialdata",
        "dev",
    }
    assert template["support"]["available_extras"] == sorted(optional)
    assert "revise_svc-0.1.0rc1-py3-none-any.whl[tacco]" in installation
    for extra in set(optional) - {"dev"}:
        assert f'python -m pip install ".[{extra}]"' in installation


def test_public_docs_use_the_single_platform_result_contract():
    text = _joined()
    cli = _read(ROOT / "revise/cli.py")
    case = _read(ROOT / "docs/source/case.rst")
    canonical_case, compatibility_case = case.split(
        "Paper notebook compatibility", maxsplit=1
    )

    for filename in ("hST-SVC.h5ad", "iST-SVC.h5ad", "sST-SVC.h5ad"):
        assert filename in text
    for compatibility_filename in (
        "sp_SVC.h5ad",
        "sc_SVC_expr.h5ad",
        "sc_SVC_spatial.h5ad",
    ):
        assert compatibility_filename not in canonical_case
        assert compatibility_filename in compatibility_case
    assert "<output-root>/<sample-name>/<platform>-SVC.h5ad" in text
    assert 'f"{args.platform}-SVC.h5ad"' in cli
    assert "provenance.json" in text


def test_ot_documentation_matches_two_stage_selection_and_failure_semantics():
    text = _joined()
    cli = _read(ROOT / "revise/cli.py")
    runtime = _read(ROOT / "revise/backend/ops/tacco_runtime.py")

    assert "--ot-method" in text
    assert "ot.ga.solver" in text
    assert "ot.lr.solver" in text
    assert "TACCO 0.5.0" in text
    assert "does not fall back" in text.lower()
    assert "annotate.mode" not in text
    assert "local_ot.method" not in text
    assert 'f"ot.ga.solver={args.ot_method}"' in cli
    assert 'f"ot.lr.solver={args.ot_method}"' in cli
    assert 'SUPPORTED_TACCO_VERSION = "0.5.0"' in runtime


def test_benchmark_metric_documentation_states_the_implemented_boundaries():
    text = " ".join(_read(ROOT / "docs/source/benchmark.rst").lower().split())
    implementation = _read(ROOT / "revise/analysis/metrics.py")

    assert "min-max" in text
    assert "total-normalize every observation in each input independently to ``1e4``" in text
    assert "sqrt(mse) / mean(ground truth)" in text
    assert "nrmse is therefore directional" in text
    assert "row order" in text
    assert "a constant normalized gene has undefined pcc and produces ``nan``" in text
    assert (
        "when normalized ground-truth mean and error are both zero, nrmse "
        "produces ``nan``" in text
    )
    assert (
        "when normalized ground-truth mean is zero but error is nonzero, nrmse "
        "produces positive infinity" in text
    )
    assert "they do not prove biological validation" in text
    assert implementation.count("target_sum=1e4") == 2
    assert "nrmse = np.sqrt(mse) / np.mean(expr_gt)" in implementation
    assert "structural_similarity as ssim" in implementation


def test_benchmark_docs_describe_actual_family_cardinality(monkeypatch, tmp_path):
    import benchmark_main

    class FakePipeline:
        def __init__(self, config_path):
            self.config_path = config_path

    monkeypatch.setattr(benchmark_main, "REVISEPipeline", FakePipeline)
    monkeypatch.setattr("builtins.print", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        benchmark_main,
        "_run_case",
        lambda *args, **kwargs: {
            "ok": True,
            "profile": kwargs["profile"],
            "run_dir": "synthetic",
            "summary": {},
            "error": None,
        },
    )

    def args_for(confounding):
        return SimpleNamespace(
            platform="sim2real",
            confounding=confounding,
            data_root=str(tmp_path / "data"),
            dataset_task=None,
            sample_name="sample",
            st_file=None,
            gt_svc_file=None,
            sc_ref_file=None,
            output_root=str(tmp_path / "output"),
            sample_size=None,
            config="revise/revise.yaml",
            seed=42,
            seed_scope="run",
            posterior_mode="off",
            posterior_key=None,
            posterior_beta=None,
            posterior_min_affinity=None,
            posterior_cost_strength=None,
            posterior_strict=False,
            sr_refinement_preset=None,
            set_overrides=[],
        )

    observed = {}
    for confounding in (
        "segmentation",
        "bin2cell",
        "batch_effect",
        "spot_size",
        "gene_panel",
        "gene_dropout",
    ):
        tags = []

        def collect(results, tag, run_result):
            tags.append(tag)
            results.append({"tag": tag, **run_result})

        monkeypatch.setattr(benchmark_main, "_append_result", collect)
        benchmark_main.main(args_for(confounding))
        observed[confounding] = tags
    assert {name: len(tags) for name, tags in observed.items()} == {
        "segmentation": 4,
        "bin2cell": 1,
        "batch_effect": 16,
        "spot_size": 4,
        "gene_panel": 1,
        "gene_dropout": 1,
    }

    for spot_size in (31, 73):
        leaf = tmp_path / "data" / "sample" / f"spot_{spot_size}"
        leaf.mkdir(parents=True)
        (leaf / "xenium_spot.h5ad").touch()
    discovered_tags = []

    def collect_discovered(results, tag, run_result):
        discovered_tags.append(tag)
        results.append({"tag": tag, **run_result})

    monkeypatch.setattr(benchmark_main, "_append_result", collect_discovered)
    benchmark_main.main(args_for("batch_effect"))
    assert len(discovered_tags) == 8
    assert {tag.split(":", 1)[1].split("_", 1)[0] for tag in discovered_tags} == {
        "31",
        "73",
    }
    benchmark = _read(ROOT / "docs/source/benchmark.rst")
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    assert "four segmentation leaves" in benchmark
    assert "four batch levels for every discovered spot size" in benchmark
    assert "four fixed spot-size leaves" in benchmark
    assert "one confounding family" in quickstart
    assert "runs one Sim2Real-ST case" not in quickstart


def test_spot_sr_and_scale_claims_do_not_exceed_candidate_evidence():
    text = " ".join(_joined().lower().split())
    limitations = " ".join(
        _read(ROOT / "docs/source/limitations.rst").lower().split()
    )
    installed_smoke = _read(ROOT / "tests/integration/test_installed_cli.py")
    graph_scale = _read(ROOT / "tests/test_graph_cluster_spatial_score.py")
    quota_scale = _read(ROOT / "tests/test_spot_sr_quota.py")

    assert "np.round" in text
    assert "exact quota" in text
    assert "seeded random" in text
    assert "virtual-cell rows" in text
    assert "not a nucleus or cell-localization result" in text
    assert "small synthetic" in text
    assert "52 observations and 52 genes" in text
    assert "50,000-observation sparse graph" in text
    assert "200,000-cell quota" in text
    assert "do not establish a complete production-scale pipeline limit" in limitations
    assert "real-data end-to-end reconstruction" in limitations
    assert "are deferred" in limitations
    assert "no current gate proves biological validation" in limitations
    assert "size=(52, 52)" in installed_smoke
    assert "n_obs = 50_000" in graph_scale
    assert "n_cells = 200_000" in quota_scale


def test_application_pm_on_cell_contract_names_the_only_resolved_location():
    concepts = _read(ROOT / "docs/source/concepts.rst")
    normalized = " ".join(concepts.split())
    runner_contract = runpy.run_path(ROOT / "revise/config/runner_conf.py")
    config = runner_contract["ApplicationScSrConf"](
        raw_data_path="data",
        result_root_path="output",
        sample_name="sample",
        cell_type_col="Level1",
        confidence_col="confidence",
        unknown_key="unknown",
        st_file="st.h5ad",
        sc_ref_file="sc_ref.h5ad",
    )

    assert config.pm_on_cell_file == "data/PM_on_cell.csv"
    assert "<data-root>/PM_on_cell.csv" in concepts
    assert "case-sensitive" in concepts
    assert "no CLI path override" in concepts
    assert "row labels must exactly match the virtual-cell IDs" in normalized
    assert "columns must exactly match the normalized cell-type labels" in normalized


def test_release_identity_and_external_data_claims_remain_owner_gated():
    text = _joined(CLAIM_DOCS)
    metadata = _read(ROOT / "pyproject.toml")
    project = tomllib.loads(metadata)["project"]
    template = json.loads(
        _read(ROOT / "release/0.1.0rc1/release-manifest.template.json")
    )
    schema = json.loads(
        _read(ROOT / "release/0.1.0rc1/release-manifest.schema.json")
    )
    lower = text.lower()

    assert not re.search(r"zenodo\.org/(?:record|records)/\d+", text)
    assert "pending owner confirmation" in lower
    assert not re.search(r"\b10\.\d{4,9}/[-._;()/:a-z0-9]+", lower)
    assert not re.search(r"\b[\w.+-]+@[\w.-]+\.[a-z]{2,}\b", lower)
    assert "corresponding author:" not in lower
    assert "scientific owner:" not in lower
    assert "zenodo.org" not in metadata
    assert project.get("authors", []) in (
        [],
        [{"name": "Pending owner confirmation"}],
    )
    assert "Data" not in project["urls"]
    source_repository = template["source"]["repository"]
    assert project["urls"]["Repository"] == source_repository
    assert (
        schema["properties"]["source"]["properties"]["repository"]["const"]
        == source_repository
    )


def test_historical_gallery_routes_to_the_legacy_asset_index():
    raw_gallery = _read(ROOT / "docs/source/gallery.rst")
    gallery = " ".join(raw_gallery.lower().split())

    assert "revise-legacy" in gallery
    assert "legacy-assets.json" in gallery
    assert "exact source commit" in gallery
    assert "not included in the clean repository" in gallery
    assert "../case/" not in raw_gallery
    assert "../benchmark/" not in raw_gallery
    assert "zenodo" not in gallery


def test_docs_version_and_support_match_validated_candidate_sources():
    conf = _read(ROOT / "docs/conf.py")
    workflow = _read(ROOT / ".github/workflows/ci.yml")
    version_source = _read(ROOT / "revise/_version.py")
    version = re.search(r'__version__ = "([^"]+)"', version_source).group(1)
    template = json.loads(
        _read(ROOT / "release/0.1.0rc1/release-manifest.template.json")
    )
    metadata = _read(ROOT / "pyproject.toml")
    requirements = _read(ROOT / "docs/requirements.txt").splitlines()

    assert "release-manifest.template.json" not in conf
    assert template["release"]["version"] == version
    assert template["docs"] == {"version": version, "channel": "candidate"}
    assert template["support"]["tested_python"] == ["3.10", "3.11"]
    assert 'requires-python = ">=3.10,<3.12"' in metadata
    assert all((ROOT / f"constraints/python-{minor}.txt").is_file() for minor in ("3.10", "3.11"))
    assert version in _read(ROOT / "README.md")
    assert "sphinx-build -W --keep-going" in workflow
    assert 'author = "Pending owner confirmation"' in conf
    assert "source/limitations" in _read(ROOT / "docs/index.rst")
    assert "nbsphinx" not in conf
    assert "plans/**" in conf
    assert "superpowers/**" in conf
    assert requirements == [
        "sphinx==8.2.3",
        "sphinx-rtd-theme==3.0.2",
        "myst-parser==4.0.1",
        "PyYAML==6.0.3",
    ]
