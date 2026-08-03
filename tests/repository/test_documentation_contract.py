from __future__ import annotations

import json
import re
import runpy
from pathlib import Path
from types import SimpleNamespace

from revise.config import load_raw_config

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]
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
    assert 'revise-reconstruct = "revise.application.cli:main"' in metadata
    assert "python reconstruct.py" in text
    assert (ROOT / "reconstruct.py").is_file()
    assert "python reproduce/benchmark_main.py" in text
    assert (ROOT / "reproduce" / "benchmark_main.py").is_file()
    assert "installed command" in text.lower()
    assert "paper reproduction" in text.lower()


def test_1x_case_notebooks_use_current_routes_and_vocabulary():
    text = []
    for path in sorted((ROOT / "reproduce" / "case").glob("*.ipynb")):
        notebook = json.loads(_read(path))
        for cell in notebook["cells"]:
            text.extend(cell.get("source", []))
            for output in cell.get("outputs", []):
                text.extend(output.get("text", []))
                data = output.get("data", {})
                text.extend(data.get("text/plain", []))
                text.extend(data.get("text/html", []))
    notebooks = "\n".join(text)

    assert "application_sc_sst" not in notebooks
    assert 'route=sST:spot_size' not in notebooks
    assert "legacy_mode" not in notebooks
    for term in ("hST", "iST", "sST"):
        assert term not in notebooks
    assert "algorithm_overrides" not in notebooks


def test_visium_case_custom_config_preserves_profile_defaults():
    notebook = json.loads(
        _read(ROOT / "reproduce/case/sc_SVC_case_Visium_mouse_brain.ipynb")
    )
    source = "\n".join(
        line
        for cell in notebook["cells"]
        for line in cell.get("source", [])
    )

    for section in ("preprocess", "graph", "sc"):
        assert f'case_profile.setdefault("{section}", {{}}).update(' in source


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
    assert "Every route requires the configured broad annotation" in normalized
    assert "Only iST-SVC requires the configured subtype annotation" in normalized
    assert (
        "sST-SVC composition and expression allocation use the broad "
        "assignment and do not require a subtype column"
    ) in normalized
    assert "default ``Patient`` column" in normalized
    assert "--spot-size" not in quickstart


def test_installation_describes_base_and_optional_capability_layers():
    installation = _read(ROOT / "docs/source/installation.rst")
    normalized = " ".join(installation.split())
    metadata = tomllib.loads(_read(ROOT / "pyproject.toml"))
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
    for extra in set(optional) - {"dev"}:
        assert f'python -m pip install ".[{extra}]"' in installation
    assert "releases can lag the repository" in installation
    assert "After a matching package version is published" in installation
    assert (
        "base package contains reconstruction, benchmarking, the POT implementation"
        in normalized
    )
    assert "iST-SVC default solver" in installation
    assert "required by the default iST-SVC route" in normalized
    assert "--ot-method pot" in installation
    assert "never selects POT as an automatic fallback" in normalized


def test_public_docs_use_the_2x_selector_and_single_file_contract():
    text = _joined()
    service = _read(ROOT / "revise/application/service.py")
    case = _read(ROOT / "docs/source/case.rst")
    configuration = _read(ROOT / "docs/source/configuration.rst")
    api_diagram = _read(ROOT / "docs/source/api/classes_revise.svg")
    test_guide = _read(ROOT / "tests/README.md")
    compatibility_heading = "Historical 1.x notebook compatibility"
    assert compatibility_heading in case
    canonical_case, compatibility_case = case.split(compatibility_heading, maxsplit=1)

    for svc_type in ("hST-SVC", "iST-SVC", "sST-SVC"):
        assert svc_type in text
    assert "sp_SVC.h5ad" not in canonical_case
    assert "sp_SVC.h5ad" in compatibility_case
    for sc_filename in ("sc_SVC_expr.h5ad", "sc_SVC_spatial.h5ad"):
        assert sc_filename not in canonical_case
        assert sc_filename in compatibility_case
    assert "<output-root>/<sample-name>/<svc-type>/SVC.h5ad" in text
    assert 'output_dir / "SVC.h5ad"' in service
    assert "Path(args.output_root) / args.sample_name / args.svc_type" in service
    assert "``result`` contains exactly ``filename`` and ``type``" in configuration
    assert "Only iST-SVC adds the top-level ``assembly`` record" in configuration
    assert "--select-ct" in text
    assert "--sc-mapping" not in text
    assert "provenance.json" in text
    assert "one route-qualified SVC.h5ad + manifest" in api_diagram
    assert "single public result file" in configuration
    assert "optional legacy merge" not in test_guide


def test_docs_label_1x_reproduction_material_without_recasting_it_as_2x_output():
    surfaces = (
        _read(ROOT / "README.md"),
        _read(ROOT / "docs/source/gallery.rst"),
        _read(ROOT / "docs/source/case.rst"),
        _read(ROOT / "docs/source/limitations.rst"),
        _read(ROOT / "reproduce/README.md"),
    )

    for surface in surfaces:
        normalized = " ".join(surface.split()).lower()
        assert "1.x historical reproduction material" in normalized
        assert "not current 2.0 output" in normalized


def test_public_docs_route_to_package_owned_application_utilities():
    text = _joined()

    assert "revise-build-histology-priors" in text
    assert "revise-compute-biological-metrics" in text


def test_ot_documentation_matches_two_stage_selection_and_failure_semantics():
    text = _joined()
    readme = _read(ROOT / "README.md")
    installation = _read(ROOT / "docs/source/installation.rst")
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    configuration = _read(ROOT / "docs/source/configuration.rst")
    normalized_docs = tuple(
        " ".join(document.split())
        for document in (readme, installation, quickstart, configuration)
    )
    service = _read(ROOT / "revise/application/service.py")
    runtime = _read(ROOT / "revise/backend/ops/tacco_runtime.py")

    assert "--ot-method" in text
    assert "ot.ga.solver" in text
    assert "ot.lr.solver" in text
    assert "TACCO 0.5.0" in text
    assert "does not fall back" in text.lower()
    assert "annotate.mode" not in text
    assert "local_ot.method" not in text
    assert '"ga": {"solver": args.ot_method}' in service
    assert '"lr": {"solver": args.ot_method}' in service
    assert 'SUPPORTED_TACCO_VERSION = "0.5.0"' in runtime
    for document in normalized_docs:
        assert "--ot-method pot" in document
    (
        normalized_readme,
        normalized_installation,
        normalized_quickstart,
        normalized_config,
    ) = normalized_docs
    assert "never falls back automatically" in normalized_readme
    assert "never selects POT as an automatic fallback" in normalized_installation
    assert "never switches algorithms automatically" in normalized_quickstart
    assert "not an automatic fallback" in normalized_config


def test_benchmark_metric_documentation_states_the_implemented_boundaries():
    text = " ".join(_read(ROOT / "docs/source/benchmark.rst").lower().split())
    implementation = _read(ROOT / "revise/analysis/metrics.py")

    assert "min-max" in text
    assert (
        "total-normalize every observation in each input independently to ``1e4``"
        in text
    )
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


def test_benchmark_docs_describe_dedicated_configuration_controls():
    text = " ".join(_read(ROOT / "docs/source/benchmark.rst").split())

    for option in (
        "--local-refinement-strength",
        "--sr-refinement-preset",
    ):
        assert option in text
    assert "applied after the selected profile or custom ``--config``" in text
    assert "Omitting the strength creates no CLI override" in text
    assert "minimal ``local_refinement`` evidence" in text
    assert "Removed policy and posterior flags are rejected" in text
    assert "--set" not in text


def test_assignment_docs_match_route_specific_runtime_contract():
    configuration = " ".join(
        _read(ROOT / "docs/source/configuration.rst").split()
    )
    architecture = " ".join(
        _read(ROOT / "docs/source/architecture.rst").split()
    )
    limitations = " ".join(
        _read(ROOT / "docs/source/limitations.rst").split()
    )

    assert "The only public local-refinement option" in configuration
    assert "hST-SVC defaults to ``0.2``" in configuration
    assert "sST-SVC defaults to ``0.0``" in configuration
    assert "uses only ``argmax(Q)``" in configuration
    assert "There are no policy" in configuration
    assert "hST-SVC conditions each local OT cost with ``Q``" in architecture
    assert "does not reweight GraphCluster with ``Q``" in architecture
    assert "``route``, ``applied``, and ``strength``" in architecture
    assert "authoritative failure explanation" in architecture
    assert "publication rollback" in architecture
    assert "not that posterior compatibility improves" in limitations


def test_benchmark_docs_describe_actual_family_cardinality(monkeypatch, tmp_path):
    from revise.benchmark import cli as benchmark_main

    class FakePipeline:
        def __init__(self, config_path):
            self.config_path = config_path
            self.raw_config = load_raw_config(ROOT / "revise/revise.yaml")

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


def test_spot_sr_and_scale_claims_do_not_exceed_current_evidence():
    text = " ".join(_joined().lower().split())
    limitations = " ".join(_read(ROOT / "docs/source/limitations.rst").lower().split())
    installed_smoke = _read(ROOT / "tests/integration/application/test_installed_cli.py")
    graph_scale = _read(ROOT / "tests/backend/test_graph_cluster_spatial_score.py")
    quota_scale = _read(ROOT / "tests/backend/test_spot_sr_quota.py")

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
    assert "have not yet been rerun" in limitations
    assert "no current gate proves biological validation" in limitations
    assert "size=(52, 52)" in installed_smoke
    assert "n_obs = 50_000" in graph_scale
    assert "n_cells = 200_000" in quota_scale


def test_application_pm_on_cell_contract_names_location_and_probability_semantics():
    concepts = _read(ROOT / "docs/source/concepts.rst")
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    normalized = " ".join(f"{concepts}\n{quickstart}".split())
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

    assert runner_contract["pm_on_cell_path_from_st_path"](
        config.st_file_path
    ) == "data/sample_st_PM_on_cell.csv"
    assert "<st-parent>/<st-stem>_PM_on_cell.csv" in normalized
    assert "case-sensitive" in concepts
    assert "no CLI path override" in concepts
    assert "sample-local probability prior" in normalized
    assert "rows must exactly equal the active virtual-cell IDs" in normalized
    assert "columns must exactly equal the active normalized cell-type labels" in normalized
    assert "numeric and finite within ``[0, 1]``" in normalized
    assert "absolute tolerance of ``1e-6``" in normalized
    assert "never clips or normalizes" in normalized
    assert "not a case table, cohort registry, or generic assignment posterior" in normalized
    assert "seeded random" in normalized

    notebook = json.loads(
        _read(ROOT / "reproduce/case/sc_SVC_case_Visium_mouse_brain.ipynb")
    )
    source = "\n".join(
        line for cell in notebook["cells"] for line in cell.get("source", [])
    )
    assert (
        'PM_ON_CELL_PATH = PREPARED_ST_PATH.with_name('
        'f"{PREPARED_ST_PATH.stem}_PM_on_cell.csv")'
    ) in source
    pm_cells = [
        cell
        for cell in notebook["cells"]
        if "PM_ON_CELL_PATH" in "".join(cell.get("source", []))
    ]
    assert len(pm_cells) == 3
    assert all(
        cell.get("execution_count") is None and cell.get("outputs") == []
        for cell in pm_cells
    )


def test_docs_state_minimal_manifest_and_publication_guarantees():
    architecture = " ".join(_read(ROOT / "docs/source/architecture.rst").split())
    limitations = " ".join(_read(ROOT / "docs/source/limitations.rst").split())
    combined = f"{architecture} {limitations}"

    assert "``running``, ``succeeded``, and ``failed``" in architecture
    assert "Captured SIGTERM and KeyboardInterrupt become failed" in architecture
    assert "uncatchable termination leaves the last manifest running" in architecture
    assert "``input_identities`` records one content identity per external role" in architecture
    assert "no aggregate data fingerprint" in architecture
    assert "no OT or Assignment event state machine" in architecture
    assert "Software identity is collected once per run" in architecture
    assert "same-directory temporary H5AD" in combined
    assert "reloads it before replacement" in combined
    assert "best-effort caught-exception rollback" in combined
    assert "not reader-atomic or crash-atomic" in combined
    assert "caller must guarantee one writer per stable public target" in combined
    assert "violating that precondition is undefined" in combined
    assert "values that are ``None`` are omitted" in combined
    assert "does not invent sentinel values" in combined
    assert "hST-SVC and sST-SVC store only ``schema_version`` and ``svc_type``" in combined
    assert "iST-SVC mean also stores ``ist_mapping`` and ``expression_source``" in combined
    assert (
        "iST-SVC random additionally stores ``effective_seed``, ``donor_column``, "
        "and ``donor_sha256``"
    ) in combined
    assert "mean assembly may retain explicit ``null`` values" in combined


def test_public_data_and_repository_claims_are_precise():
    text = _joined(CLAIM_DOCS)
    metadata = _read(ROOT / "pyproject.toml")
    project = tomllib.loads(metadata)["project"]
    lower = text.lower()

    assert "https://zenodo.org/records/17705737" in text
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
    assert project["urls"]["Repository"] == "https://github.com/wuys13/REVISE"


def test_gallery_describes_curated_notebooks_and_their_evidence_boundary():
    raw_gallery = _read(ROOT / "docs/source/gallery.rst")
    gallery = " ".join(raw_gallery.lower().split())

    assert "reproduce/benchmark/" in gallery
    assert "reproduce/case/" in gallery
    assert "not part of the installed python package" in gallery
    assert "does not establish" in gallery
    assert "../case/" not in raw_gallery
    assert "../benchmark/" not in raw_gallery
    assert "zenodo" not in gallery


def test_docs_version_and_support_match_project_sources():
    conf = _read(ROOT / "docs/conf.py")
    workflow = _read(ROOT / ".github/workflows/ci.yml")
    metadata = _read(ROOT / "pyproject.toml")
    requirements = _read(ROOT / "docs/requirements.txt").splitlines()

    assert "release-manifest.template.json" not in conf
    assert 'requires-python = ">=3.10,<3.12"' in metadata
    assert all(
        (ROOT / f"constraints/python-{minor}.txt").is_file()
        for minor in ("3.10", "3.11")
    )
    assert 'release = version' in conf
    assert "sphinx-build -W --keep-going" in workflow
    assert 'author = "Pending owner confirmation"' in conf
    assert "source/limitations" in _read(ROOT / "docs/index.rst")
    assert "nbsphinx" not in conf
    assert "plans/**" not in conf
    assert "superpowers/**" not in conf
    assert requirements == [
        "sphinx==8.2.3",
        "sphinx-rtd-theme==3.0.2",
        "myst-parser==4.0.1",
        "PyYAML==6.0.3",
    ]


def test_tacco_smoke_runs_both_solver_contracts_against_candidate_wheel():
    workflow = _read(ROOT / ".github/workflows/ci.yml")
    tacco_job = workflow.split("  tacco-smoke:", maxsplit=1)[1]

    assert "needs: package" in tacco_job
    assert "actions/download-artifact@v4" in tacco_job
    assert "name: candidate-dist" in tacco_job
    assert "revise_svc-2.0.0rc1-py3-none-any.whl" in tacco_job
    assert "working-directory: ${{ runner.temp }}" in tacco_job
    assert 'python -m venv "${RUNNER_TEMP}/candidate-venv"' in tacco_job
    assert '"${RUNNER_TEMP}/candidate-venv/bin/python" -m pip install' in tacco_job
    assert "wheel-smoke-tests/tests/integration/solvers" in tacco_job
    assert "tests/integration/solvers/conftest.py" in tacco_job
    assert "REVISE_EXPECT_INSTALLED_PREFIX" in tacco_job
    assert "--import-mode=importlib" in tacco_job
    assert '"${GITHUB_WORKSPACE}/tests/integration/solvers/' not in tacco_job
    assert "tests/integration/solvers/test_tacco_solver_smoke.py" in tacco_job
    assert (
        "tests/integration/solvers/test_local_refinement_solver_smoke.py"
        in tacco_job
    )
