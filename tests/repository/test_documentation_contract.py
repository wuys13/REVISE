from __future__ import annotations

import json
import re
import runpy
from pathlib import Path
from types import SimpleNamespace

from revise.config.authority import _authority_document

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
    ROOT / "docs/source/concepts.rst",
    ROOT / "docs/source/architecture.rst",
)
APPLICATION_DOCS = (
    ROOT / "README.md",
    ROOT / "docs/index.rst",
    ROOT / "docs/source/installation.rst",
    ROOT / "docs/source/quickstart.rst",
    ROOT / "docs/source/concepts.rst",
)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _joined(paths=CLAIM_DOCS) -> str:
    return "\n".join(_read(path) for path in paths)


def test_public_docs_distinguish_installed_source_and_paper_entry_paths():
    text = _joined()
    metadata = _read(ROOT / "pyproject.toml")

    assert "revise-reconstruct" in text
    assert 'revise-reconstruct = "reconstruct:main"' in metadata
    assert "python reconstruct.py" in text
    assert (ROOT / "reconstruct.py").is_file()
    assert "python reproduce/benchmark_main.py" in text
    assert (ROOT / "reproduce" / "benchmark_main.py").is_file()
    assert "installed command" in text.lower()
    assert "paper" in text.lower()
    assert "reproduce" in text.lower()


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


def test_visium_case_uses_the_canonical_reconstruct_entrypoint():
    notebook = json.loads(
        _read(ROOT / "reproduce/case/sc_SVC_sr_case_Visium_mouse_brain.ipynb")
    )
    source = "\n".join(
        line
        for cell in notebook["cells"]
        for line in cell.get("source", [])
    )

    assert "REVISEPipeline" not in source
    assert "application_sc_sr" not in source
    assert "reconstruct.py --config configs/application/Visium.yaml" in source
    assert "subprocess.run" in source
    assert "cwd=REPO_ROOT" in source
    assert "PM_ON_CELL_PATH" in source


def test_quickstart_matches_application_yaml_entry_and_route_fields():
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    normalized = " ".join(quickstart.split())

    for filename in (
        "Xenium_T.yaml",
        "Xenium_Fib.yaml",
        "Xenium_Mono.yaml",
        "VisiumHD.yaml",
        "Visium.yaml",
    ):
        assert f"configs/application/{filename}" in quickstart
    assert "python reconstruct.py --config" in quickstart
    assert "revise-reconstruct --config" in quickstart
    assert "inputs.st.path" in quickstart
    assert "inputs.reference.path" in quickstart
    assert "unique ``obs_names`` and unique ``var_names``" in normalized
    assert "Every route requires the configured broad annotation" in normalized
    assert "Only standard sc-SVC requires the configured subtype annotation" in normalized
    assert (
        "sc-SVC-sr`` composition and expression allocation use the broad assignment"
        in normalized
    )
    assert "do not require a subtype column" in normalized
    normalized_lower = normalized.lower()
    assert "no implicit patient filter" in normalized_lower
    assert "filter runs only when the yaml declares both" in normalized_lower
    assert "explicitly select" in normalized_lower
    assert "before spatial and reference preprocessing" in normalized_lower
    assert "--spot-size" not in quickstart
    for removed in ("--svc-type", "--st-file", "--select-ct", "--ot-method"):
        assert removed not in quickstart


def test_docs_name_typed_authority_and_all_eleven_route_yamls():
    text = _joined()
    assert "revise.config.authority" in text
    for filename in ("Xenium_T.yaml", "Xenium_Fib.yaml", "Xenium_Mono.yaml", "VisiumHD.yaml", "Visium.yaml"):
        assert filename in text
    for filename in (
        "segmentation.yaml",
        "bin2cell.yaml",
        "batch_effect.yaml",
        "spot_size.yaml",
        "gene_panel.yaml",
        "gene_dropout.yaml",
    ):
        assert (ROOT / "configs/benchmark" / filename).is_file()


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
    assert "Standard sc-SVC default solver" in installation
    assert "required by the default standard sc-SVC route" in normalized
    assert "algorithm.ot_method: pot" in installation
    assert "never selects POT as an automatic fallback" in normalized


def test_public_docs_use_route_specific_single_and_sc_pair_contracts():
    text = _joined()
    publication = _read(ROOT / "revise/application/publication.py")
    api_diagram = _read(ROOT / "docs/source/api/classes_revise.svg")
    test_guide = _read(ROOT / "tests/README.md")

    assert "<output-dir>/<output-name>.h5ad" in text
    assert 'f"{config.output_name}.h5ad"' in publication
    assert 'f"{prefix}expr.h5ad"' in publication
    assert 'f"{prefix}spatial.h5ad"' in publication
    assert "sp-SVC" in text
    assert "sc-SVC" in text
    assert "sc-SVC-sr" in text
    assert "provenance.json" in text
    assert "H5AD result(s) + manifest" in api_diagram
    assert "SVC.h5ad + manifest type" not in api_diagram
    assert "optional legacy merge" not in test_guide


def test_documentation_navigation_has_six_benchmark_notebooks_and_gallery_titles():
    index = _read(ROOT / "docs/index.rst")

    assert index.startswith(
        ".. include:: ../README.md\n   :parser: readme_parser\n"
    )
    assert (ROOT / "docs/readme_parser.py").is_file()
    assert "REVISE Documentation" not in index
    captions = re.findall(r":caption:\s*(.+?)\s*$", index, flags=re.MULTILINE)
    assert captions == ["START HERE", "Sim2Real Benchmark", "Gallery", "Reference"]

    assert "Benchmark Mode" not in index
    assert "source/benchmark_mode" not in index
    assert not (ROOT / "docs/source/benchmark_mode.rst").exists()
    benchmark_block = index.split(":caption: Sim2Real Benchmark", maxsplit=1)[1].split(
        ":caption: Gallery", maxsplit=1
    )[0]
    assert benchmark_block.count("<benchmark/") == 6
    benchmark_entries = [
        line.strip()
        for line in benchmark_block.splitlines()
        if "<benchmark/" in line
    ]
    expected_benchmarks = [
        ("segmentation", "benchmark/segmentation"),
        ("bin2cell", "benchmark/bin2cell"),
        ("batch", "benchmark/batch"),
        ("spot", "benchmark/spot"),
        ("imputation_and_dropout", "benchmark/imputation_and_dropout"),
        ("plot_imputation_case", "benchmark/plot_imputation_case"),
    ]
    assert benchmark_entries == [
        f"{title} <{target}>" for title, target in expected_benchmarks
    ]
    for title, target in expected_benchmarks:
        link = json.loads(_read(ROOT / "docs" / f"{target}.nblink"))
        link_path = ROOT / "docs" / target
        notebook_path = link_path.parent.joinpath(link["path"]).resolve()
        assert notebook_path.suffix == ".ipynb"
        assert title == notebook_path.stem
    assert "source/case" not in index

    gallery_titles = (
        "Visium (sST platform) sc-SVC-sr mouse brain <case/sc_SVC_sr_case_Visium_mouse_brain>",
        "Xenium (iST platform) sc-SVC Fibroblast <case/Xenium_sc_SVC_Fibroblast>",
        "Xenium (iST platform) sc-SVC Mono/Macro <case/Xenium_sc_SVC_Monocyte>",
        "Xenium (iST platform) sc-SVC T cells <case/Xenium_sc_SVC_T>",
        "VisiumHD (hST platform) sp-SVC <case/VisiumHD_sp_SVC>",
    )
    gallery_block = index.split(":caption: Gallery", maxsplit=1)[1].split(
        ":caption: Reference", maxsplit=1
    )[0]
    assert [
        line.strip()
        for line in gallery_block.splitlines()
        if "<case/" in line
    ] == list(gallery_titles)

    assert "source/gallery" not in index
    assert not re.search(r"/Users/|/home/", index)

    for path in CLAIM_DOCS:
        assert not re.search(r"/Users/|/home/", _read(path))


def test_gallery_nblink_inventory_is_exact_and_targets_reproduce_notebooks():
    expected = {
        "docs/benchmark/segmentation.nblink":
        "reproduce/benchmark/segmentation.ipynb",
        "docs/benchmark/bin2cell.nblink":
        "reproduce/benchmark/bin2cell.ipynb",
        "docs/benchmark/batch.nblink":
        "reproduce/benchmark/batch.ipynb",
        "docs/benchmark/spot.nblink":
        "reproduce/benchmark/spot.ipynb",
        "docs/benchmark/imputation_and_dropout.nblink":
        "reproduce/benchmark/imputation_and_dropout.ipynb",
        "docs/benchmark/plot_imputation_case.nblink":
        "reproduce/benchmark/plot_imputation_case.ipynb",
        "docs/case/VisiumHD_sp_SVC.nblink":
        "reproduce/case/VisiumHD_sp_SVC.ipynb",
        "docs/case/Xenium_sc_SVC_T.nblink":
        "reproduce/case/Xenium_sc_SVC_T.ipynb",
        "docs/case/Xenium_sc_SVC_Monocyte.nblink":
        "reproduce/case/Xenium_sc_SVC_Monocyte.ipynb",
        "docs/case/Xenium_sc_SVC_Fibroblast.nblink":
        "reproduce/case/Xenium_sc_SVC_Fibroblast.ipynb",
        "docs/case/sc_SVC_sr_case_Visium_mouse_brain.nblink":
        "reproduce/case/sc_SVC_sr_case_Visium_mouse_brain.ipynb",
    }
    observed = {
        str(path.relative_to(ROOT)): path
        for directory in (ROOT / "docs/benchmark", ROOT / "docs/case")
        for path in sorted(directory.glob("*.nblink"))
    }
    assert set(observed) == set(expected)
    for relative, target in expected.items():
        payload = json.loads(_read(ROOT / relative))
        assert set(payload) == {"path"}
        assert payload["path"] == f"../../{target}"
        assert (ROOT / target).is_file()
        assert (ROOT / relative).parent.joinpath(payload["path"]).resolve() == ROOT / target

    assert not any("_recon" in path.name for path in observed.values())
    assert not any(
        "application_sc_SVC_analysis_case" in path.name for path in observed.values()
    )

def test_public_docs_route_to_package_owned_application_utilities():
    text = _joined()

    assert "revise-build-histology-priors" in text
    assert "revise-compute-biological-metrics" in text


def test_ot_documentation_matches_two_stage_selection_and_failure_semantics():
    text = _joined()
    readme = _read(ROOT / "README.md")
    installation = _read(ROOT / "docs/source/installation.rst")
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    normalized_docs = tuple(
        " ".join(document.split())
        for document in (readme, installation, quickstart)
    )
    engine_compiler = _read(ROOT / "revise/application/config.py")
    runtime = _read(ROOT / "revise/backend/ops/tacco_runtime.py")

    assert "algorithm.ot_method" in text
    assert "ot.ga.solver" in text
    assert "ot.lr.solver" in text
    assert "TACCO 0.5.0" in text
    assert "fallback" in text.lower()
    assert "annotate.mode" not in text
    assert "local_ot.method" not in text
    assert '"ga": {"solver": config.ot_method}' in engine_compiler
    assert '"lr": {"solver": config.ot_method}' in engine_compiler
    assert 'SUPPORTED_TACCO_VERSION = "0.5.0"' in runtime
    for document in normalized_docs:
        assert "algorithm.ot_method" in document
    normalized_readme, normalized_installation, normalized_quickstart = normalized_docs
    assert "never changes solvers automatically" in normalized_readme
    assert "never selects POT as an automatic fallback" in normalized_installation
    assert "never switches algorithms automatically" in normalized_quickstart


def test_assignment_docs_match_route_specific_runtime_contract():
    architecture = " ".join(
        _read(ROOT / "docs/source/architecture.rst").split()
    )

    assert "sp-SVC conditions each local OT cost with ``Q``" in architecture
    assert "does not reweight GraphCluster with ``Q``" in architecture
    assert "``route``, ``applied``, and ``strength``" in architecture
    assert "authoritative failure explanation" in architecture
    assert "same-directory temporary H5AD" in architecture


def test_benchmark_docs_describe_actual_family_cardinality(monkeypatch, tmp_path):
    from revise.benchmark import cli as benchmark_main

    class FakePipeline:
        def __init__(self):
            self.raw_config = _authority_document()

    monkeypatch.setattr(benchmark_main, "REVISEPipeline", FakePipeline)
    monkeypatch.setattr("builtins.print", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        benchmark_main,
        "_run_case",
        lambda *args, **kwargs: {
            "ok": True,
            "profile": "synthetic",
            "seed": kwargs["runtime_seed"],
            "run_dir": "synthetic",
            "manifest_path": "synthetic/provenance.json",
            "local_refinement": None,
            "summary": {},
            "error": None,
        },
    )

    def args_for(route):
        return SimpleNamespace(
            data_root=str(tmp_path / "data"),
            dataset_task=None,
            sample_name="sample",
            output_root=str(tmp_path / "output"),
            sample_size=None,
            config=str(ROOT / "configs" / "benchmark" / f"{route}.yaml"),
            seed=42,
            seed_scope="run",
            local_refinement_strength=None,
            sr_refinement_preset=None,
            evaluate=True,
        )

    observed = {}
    for route in (
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
        benchmark_main.main(args_for(route))
        observed[route] = tags
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
    assert len(discovered_tags) == 16
    assert {tag.split(":", 1)[1].split("_", 1)[0] for tag in discovered_tags} == {
        "50",
        "100",
        "150",
        "200",
    }
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    assert "Benchmark" in quickstart
    assert "runs one Sim2Real-ST case" not in quickstart


def test_spot_sr_and_scale_claims_do_not_exceed_current_evidence():
    text = " ".join(_joined().lower().split())
    installed_smoke = _read(ROOT / "tests/integration/application/test_installed_cli.py")
    graph_scale = _read(ROOT / "tests/backend/test_graph_cluster_spatial_score.py")
    quota_scale = _read(ROOT / "tests/backend/test_spot_sr_quota.py")

    assert "np.round" in text
    assert "spot-level global anchoring posterior proportions" in text
    assert "lr proportions" not in text
    assert "exact quota" in text
    assert "seeded random" in text
    assert "virtual-cell rows" in text
    assert "not a nucleus or cell-localization result" in text
    assert "small synthetic" in text
    assert "size=(52, 52)" in installed_smoke
    assert "n_obs = 50_000" in graph_scale
    assert "n_cells = 200_000" in quota_scale


def test_application_pm_on_cell_contract_names_location_and_probability_semantics():
    concepts = _read(ROOT / "docs/source/concepts.rst")
    quickstart = _read(ROOT / "docs/source/quickstart.rst")
    normalized = " ".join(f"{concepts}\n{quickstart}".split())
    runner_contract = runpy.run_path(ROOT / "revise/config/runner_conf.py")

    assert (
        runner_contract["pm_on_cell_path_from_data_root"]("data")
        == "data/PM_on_cell.csv"
    )
    assert "exact ``inputs.pm_on_cell.path``" in normalized
    assert "no sidecar is probed" in normalized
    assert "sample-local score matrix" in normalized
    assert "rows exactly equal to the active virtual-cell IDs" in normalized
    assert "columns exactly equal to the active normalized broad labels" in normalized
    assert "numeric finite values in ``[0, 1]``" in normalized
    assert "never clips or normalizes" in normalized
    assert "not a case table, cohort registry, or generic assignment posterior" in normalized
    assert "seeded random" in normalized



def test_docs_state_minimal_manifest_and_publication_guarantees():
    architecture = " ".join(_read(ROOT / "docs/source/architecture.rst").split())
    concepts = " ".join(_read(ROOT / "docs/source/concepts.rst").split())
    combined = architecture

    assert "``running``, ``succeeded``, and ``failed``" in architecture
    assert "Captured SIGTERM and KeyboardInterrupt become failed" in architecture
    assert "uncatchable termination leaves the last manifest running" in architecture
    assert "``input_identities`` records one content identity per external role" in architecture
    assert "no aggregate data fingerprint" in architecture
    assert "no OT or Assignment event state machine" in architecture
    assert "solver events" not in concepts
    assert "Software identity is collected once per run" in architecture
    assert "same-directory temporary H5AD" in combined
    assert "does not provide rollback" in combined
    assert "not reader-atomic or crash-atomic" in combined
    assert "caller must guarantee one writer per stable public target" in combined
    assert "violating that precondition is undefined" in combined


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
    assert "not executed during documentation builds" in gallery
    assert ".nblink" in gallery
    assert "does not establish" in gallery
    assert "../case/" not in raw_gallery
    assert "../benchmark/" not in raw_gallery
    assert "zenodo" not in gallery


def test_docs_version_and_support_match_project_sources():
    conf = _read(ROOT / "docs/conf.py")
    workflow = _read(ROOT / ".github/workflows/ci.yml")
    readthedocs = _read(ROOT / ".readthedocs.yaml")
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
    assert 'nbsphinx_execute = "never"' in conf
    assert '"nbsphinx"' in conf
    assert '"nbsphinx_link"' in conf
    assert '"titles_only": True' in conf
    assert "plans/**" not in conf
    assert "superpowers/**" not in conf
    assert "nbsphinx==0.9.8" in requirements
    assert "nbsphinx-link==1.4.1" in requirements
    assert "ipython==8.37.0" in requirements
    assert "sphinx==8.1.3" in requirements
    assert "- requirements: docs/requirements.txt" in readthedocs
    runtime_metadata = metadata.split("[project.scripts]", maxsplit=1)[0].lower()
    for docs_dependency in (
        "sphinx",
        "sphinx-rtd-theme",
        "myst-parser",
        "nbsphinx",
        "nbsphinx-link",
        "ipython",
    ):
        assert docs_dependency not in runtime_metadata


def test_tacco_smoke_runs_both_solver_contracts_against_candidate_wheel():
    workflow = _read(ROOT / ".github/workflows/ci.yml")
    installed_job = workflow.split("  installed-cli:", maxsplit=1)[1].split(
        "  docs:", maxsplit=1
    )[0]
    tacco_job = workflow.split("  tacco-smoke:", maxsplit=1)[1]

    for candidate_job in (installed_job, tacco_job):
        assert 'wheels=("${RUNNER_TEMP}"/candidate/revise_svc-*.whl)' in candidate_job
        assert 'test "${#wheels[@]}" -eq 1' in candidate_job
        assert 'test -f "${wheels[0]}"' in candidate_job
        assert '"${wheels[0]}"' in candidate_job
    assert re.search(r"revise_svc-\d", workflow) is None
    assert 'echo "REVISE_WHEEL=${wheels[0]}" >> "${GITHUB_ENV}"' in installed_job
    assert "needs: package" in tacco_job
    assert "actions/download-artifact@v4" in tacco_job
    assert "name: candidate-dist" in tacco_job
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


def test_release_workflow_uses_oidc_and_separates_release_asset_permissions():
    workflow = _read(ROOT / ".github/workflows/release.yml")
    publish_action = (
        "pypa/gh-action-pypi-publish@"
        "dc37677b2e1c63e2034f94d8a5b11f265b73ba33"
    )
    node24_actions = {
        "actions/checkout@fbc6f3992d24b796d5a048ff273f7fcc4a7b6c09": 5,
        "actions/setup-python@ece7cb06caefa5fff74198d8649806c4678c61a1": 6,
        "actions/upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a": 1,
        "actions/download-artifact@37930b1c2abaa49bbe596cd826c3c89aef350131": 7,
    }

    assert "release:\n    types: [published]" in workflow
    action_refs = re.findall(r"uses:\s+[A-Za-z0-9_.\-/]+@([^\s#]+)", workflow)
    assert action_refs
    assert all(re.fullmatch(r"[0-9a-f]{40}", ref) for ref in action_refs)
    assert len(action_refs) == 21
    for action, count in node24_actions.items():
        assert workflow.count(action) == count
    assert workflow.count("python -m build") == 1
    assert workflow.count(publish_action) == 2
    assert workflow.count("id-token: write") == 2
    assert workflow.count("ref: ${{ github.event.release.tag_name }}") == 5
    assert workflow.count("persist-credentials: false") == 5
    assert "fetch-depth: 0" in workflow
    assert "git merge-base --is-ancestor HEAD refs/remotes/origin/main" in workflow
    assert "name: testpypi" in workflow
    assert "name: pypi" in workflow
    assert "repository-url: https://test.pypi.org/legacy/" in workflow
    assert "password:" not in workflow
    assert "skip-existing" not in workflow
    assert workflow.count("for attempt in range(12)") == 2
    assert "needs: [build, verify-pypi-files]" in workflow
    assert "contents: write" in workflow
    assert 'gh release upload "${RELEASE_TAG}"' in workflow
    assert 'GH_REPO: ${{ github.repository }}' in workflow
    assert "--clobber" not in workflow
