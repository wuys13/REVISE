from __future__ import annotations

import json
import re
from pathlib import Path

import yaml

from revise.config.authority import _authority_document

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]
USER_DOCS = (
    ROOT / "README.md",
    ROOT / "docs/index.rst",
    ROOT / "docs/source/quickstart.rst",
    ROOT / "docs/source/application-reference.rst",
    ROOT / "docs/source/concepts.rst",
    ROOT / "docs/source/installation.rst",
    ROOT / "docs/source/architecture.rst",
    ROOT / "docs/source/api/index.rst",
    ROOT / "docs/source/gallery.rst",
    ROOT / "reproduce/README.md",
    ROOT / "tests/README.md",
)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _joined(paths=USER_DOCS) -> str:
    return "\n".join(_read(path) for path in paths)


def _markdown_section(text: str, heading: str) -> str:
    marker = f"## {heading}\n"
    return text.split(marker, 1)[1].split("\n## ", 1)[0]


def test_readme_first_run_flow_has_resources_install_downloads_and_quick_run():
    readme = _read(ROOT / "README.md")

    assert (
        "[REVISE documentation](https://revise-svc.readthedocs.io/en/latest/) | "
        "[Dataset](https://zenodo.org/records/21921802) | Paper"
    ) in readme
    assert "| ST data | Start from | Public route |" in readme
    assert "| Result |" not in readme
    assert readme.index("## Quick run") < readme.index(
        "### Choose the application template that matches the ST data"
    )

    install = _markdown_section(readme, "Install")
    visible_install, optional_install = install.split("<details>", 1)
    assert visible_install.count("```bash") == 1
    assert "python -m pip install revise-svc" in visible_install
    assert 'python -m pip install "revise-svc[tacco]"' in visible_install
    assert "POT" in visible_install
    assert "TACCO" in visible_install
    assert optional_install.count("```bash") == 1
    assert 'python -m pip install "revise-svc[tacco]"' not in optional_install
    for command in (
        'python -m pip install "revise-svc[pathway]"',
        'python -m pip install "revise-svc[cci]"',
        'python -m pip install "revise-svc[trajectory]"',
        'python -m pip install "revise-svc[spatialdata]"',
    ):
        assert command in optional_install
    assert "[dev]" not in install
    assert "[docs]" not in install
    assert "data-reading or downstream bioinformatics" in optional_install

    assert readme.index("## Install") < readme.index("## Data downloads") < readme.index("## Quick run")
    downloads = _markdown_section(readme, "Data downloads")
    for url, label in (
        ("https://zenodo.org/records/21921802", "Sim2Real-ST benchmark"),
        ("https://zenodo.org/records/21921802", "Reproduced benchmark results"),
        ("https://zenodo.org/records/21921802", "Real-world ST datasets"),
        ("https://zenodo.org/records/18389835", "Reproduced sp-SVC H5AD"),
        ("https://zenodo.org/records/18389211", "Reproduced sc-SVC H5AD"),
    ):
        assert url in downloads
        assert label in downloads
    assert "without configuration" not in downloads.lower()


def test_readme_quick_run_and_outputs_follow_current_cluster_directory_contract():
    readme = _read(ROOT / "README.md")
    quick_run = _markdown_section(readme, "Quick run")

    assert "## Configure" not in readme
    assert quick_run.count("```bash") == 3
    assert "local editable YAML copy" in quick_run
    assert "revise-reconstruct --config VisiumHD.yaml" in quick_run
    assert "revise-reconstruct --config Xenium.yaml --select-ct T" in quick_run
    assert "revise-reconstruct --config Visium.yaml" in quick_run
    assert "Source-checkout equivalent" in quick_run
    for name in ("VisiumHD.yaml", "Xenium.yaml", "Visium.yaml"):
        assert f"configs/application/{name}" in quick_run
    for command in (
        "python reconstruct.py --config configs/application/Xenium.yaml --select-ct T",
        "python reconstruct.py --config configs/application/Xenium.yaml --select-ct Fibroblast",
        'python reconstruct.py --config configs/application/Xenium.yaml --select-ct "Mono/Macro"',
    ):
        assert command in quick_run
    assert (
        "input paths, reference filter, annotation columns, preprocessing thresholds, and output root"
        in " ".join(quick_run.split())
    )
    assert (
        "automatically writes to `T/`, `Fibroblast/`, or `Mono_Macro/`"
        in " ".join(quick_run.split())
    )

    outputs = _markdown_section(readme, "What is written")
    assert "`svc.h5ad` when `output.name` is omitted" in outputs
    assert "`output.dir` is the base directory" in outputs
    assert "final selected-cell-type subdirectory" in " ".join(outputs.split())


def test_application_reference_is_the_only_owner_of_cluster_output_directory_rules():
    reference = _read(ROOT / "docs/source/application-reference.rst")
    normalized_reference = " ".join(reference.split())

    for text in (
        "``preprocessing.spatial.min_counts`` is optional",
        "``preprocessing.reference.min_genes`` is optional",
        "``null`` disables either optional threshold",
        "``execution.seed``",
        "``output.dir`` is the base directory",
        "normalized selected cell-type label",
        "safe for an output directory",
    ):
        assert text in normalized_reference

    for path in (
        ROOT / "docs/source/quickstart.rst",
        ROOT / "docs/source/application-migration.rst",
        ROOT / "docs/source/api/index.rst",
        ROOT / "reproduce/README.md",
    ):
        text = _read(path)
        assert "does not alter ``output.dir``" not in text
        assert "does not change an output path" not in text
        assert "distinct ``output.dir``" not in text


def test_root_reproduction_entry_is_bounded_and_routes_to_detailed_guides():
    readme = _read(ROOT / "README.md")
    reproduce = _markdown_section(readme, "Reproduce published material")

    assert "### Sim2Real-ST Benchmark" in reproduce
    assert "### Real-world ST Application" in reproduce
    assert "BENCHMARK_MAX_JOBS=1" in reproduce
    assert "bash reproduce/benchmark_main.sh" in reproduce
    assert "reproduce/benchmark/" in reproduce
    assert "reproduce/case/" in reproduce
    for link in (
        "[reproduce/README.md](reproduce/README.md)",
        "[Benchmark documentation](https://revise-svc.readthedocs.io/en/latest/)",
        "[Application gallery](https://revise-svc.readthedocs.io/en/latest/source/gallery.html)",
        "[REVISE documentation](https://revise-svc.readthedocs.io/en/latest/)",
    ):
        assert link in reproduce


def test_readme_is_a_first_run_guide_with_the_documentation_link():
    readme = _read(ROOT / "README.md")

    assert "[REVISE documentation](https://revise-svc.readthedocs.io/en/latest/)" in readme
    assert "Choose the application template that matches the ST data" in readme
    assert readme.index("VisiumHD.yaml") < readme.index("Xenium.yaml") < readme.index("Visium.yaml")
    assert "Python 3.10 and 3.11" in readme
    assert "--select-ct T" in readme
    assert "What is written" in readme
    assert "REVISE framework overview" in readme
    assert "Reproduced benchmark results" in readme
    assert "## Python API" not in readme
    assert "## Documentation index" in readme
    assert "## Repository layout" in readme
    assert "## Citation and license" in readme
    assert "source_sha256" not in readme
    assert "authority_hash" not in readme
    assert "engine_defaults_hash" not in readme


def test_readme_documentation_index_points_to_canonical_owners():
    readme = _read(ROOT / "README.md")
    for page in (
        "quickstart",
        "application-reference",
        "concepts",
        "installation",
        "application-migration",
        "architecture",
        "gallery",
    ):
        assert f"https://revise-svc.readthedocs.io/en/latest/source/{page}.html" in readme
    for path in (
        ROOT / "docs/source/application-reference.rst",
        ROOT / "docs/source/application-migration.rst",
        ROOT / "docs/source/api/index.rst",
        ROOT / "reproduce/README.md",
        ROOT / "tests/README.md",
    ):
        assert path.is_file()


def test_application_templates_are_exactly_three_identical_public_requests():
    source_templates = ROOT / "configs/application"
    package_templates = ROOT / "revise/application/templates"
    expected = {"VisiumHD.yaml", "Xenium.yaml", "Visium.yaml"}

    assert {path.name for path in source_templates.glob("*.yaml")} == expected
    assert {path.name for path in package_templates.glob("*.yaml")} == expected
    for name in expected:
        assert (source_templates / name).read_bytes() == (package_templates / name).read_bytes()

    requests = {
        name: yaml.safe_load((source_templates / name).read_text(encoding="utf-8"))
        for name in expected
    }
    assert requests["VisiumHD.yaml"]["application"] == {"svc_type": "sp-SVC"}
    assert requests["Xenium.yaml"]["application"] == {
        "svc_type": "sc-SVC",
        "mode": "cluster",
    }
    assert requests["Visium.yaml"]["application"] == {
        "svc_type": "sc-SVC",
        "mode": "sr",
    }
    assert requests["Xenium.yaml"]["local_refinement"]["select_cell_type"] == "T"


def test_current_user_docs_use_sc_svc_modes_not_a_third_public_category():
    current = _joined()

    assert "sc-SVC-sr" not in current
    assert "Xenium_T.yaml" not in current
    assert "Xenium_Fib.yaml" not in current
    assert "Xenium_Mono.yaml" not in current
    assert "sc_SVC_sr_case_Visium_mouse_brain" not in current
    assert "sc-SVC``, ``cluster`` mode" in _read(ROOT / "docs/source/quickstart.rst")
    assert "sc-SVC``, ``sr`` mode" in _read(ROOT / "docs/source/quickstart.rst")

    migration = _read(ROOT / "docs/source/application-migration.rst")
    for legacy in (
        "sc-SVC-sr",
        "Xenium_T.yaml",
        "application_sc_sr",
        "sc_svc_sr_application",
        "ApplicationScSrConf",
    ):
        assert legacy in migration


def test_quickstart_matches_the_current_cli_and_minimum_fields():
    quickstart = _read(ROOT / "docs/source/quickstart.rst")

    for name in ("VisiumHD.yaml", "Xenium.yaml", "Visium.yaml"):
        assert f"configs/application/{name}" in quickstart
    assert "python reconstruct.py --config" in quickstart
    assert "revise-reconstruct --config" in quickstart
    assert "--select-ct T" in quickstart
    assert "inputs.st.path" in quickstart
    assert "inputs.reference.path" in quickstart
    assert "global_anchoring.broad_column" in quickstart
    assert "local_refinement.subtype_column" in quickstart
    assert "output.dir" in quickstart
    assert "all`` and wildcards are rejected" in quickstart
    checkout = quickstart.index("From a source checkout")
    assert ":ref:`Application templates <application-templates>`" in quickstart
    assert quickstart.index("revise-reconstruct --config VisiumHD.yaml") < checkout
    assert quickstart.index("python reconstruct.py --config") > checkout


def test_application_reference_matches_public_config_boundaries():
    reference = _read(ROOT / "docs/source/application-reference.rst")

    for field in (
        "application.mode",
        "inputs.st.spatialdata.table",
        "filter_column",
        "broad_column",
        "subtype_column",
        "preprocessing.spatial",
        "algorithm.ot_method",
        "inputs.pm_on_cell.path",
        "provenance.json",
    ):
        assert field in reference
    assert "--select-ct VALUE" in reference
    assert "does not substitute" in reference
    assert "one H5AD" in reference
    assert "spatial.h5ad" in reference
    assert "expr.h5ad" in reference


def test_public_docs_distinguish_application_and_benchmark_entrypoints():
    text = _joined()
    metadata = _read(ROOT / "pyproject.toml")

    assert "revise-reconstruct" in text
    assert 'revise-reconstruct = "reconstruct:main"' in metadata
    assert "python reconstruct.py" in text
    assert "python reproduce/benchmark_main.py" in text
    assert "Benchmark YAML" in text
    assert "does not accept Application modes" in text
    authority = _authority_document()
    assert set(authority["router"]["application"]) == {
        "sp-SVC",
        "sc-SVC:cluster",
        "sc-SVC:sr",
    }
    assert authority["router"]["benchmark"]["batch_effect"]["task"] == "sc_svc_sr"


def test_architecture_records_the_application_only_sr_rename_and_failure_contract():
    architecture = _read(ROOT / "docs/source/architecture.rst")

    for value in (
        "application_sc_super_resolution",
        "sc_svc_super_resolution",
        "ScSvcSuperResolutionApplicationStrategy",
        "ApplicationScSuperResolutionConf",
        "sc_svc_super_resolution_application.ScSVCSuperResolution",
        "application_route: sc-SVC",
        "application_mode: cluster",
        "application_mode: sr",
        "Benchmark retains its own ``sc_svc_sr`` task",
        "temporary H5AD",
        "reader-atomic or crash-atomic",
        "``running``, ``succeeded``, or ``failed``",
        "authoritative failure explanation",
    ):
        assert value in architecture


def test_documentation_navigation_has_six_benchmark_pages_and_ordered_gallery():
    index = _read(ROOT / "docs/index.rst")

    assert index.startswith(".. include:: ../README.md\n   :parser: readme_parser\n")
    captions = re.findall(r":caption:\s*(.+?)\s*$", index, flags=re.MULTILINE)
    assert captions == ["START HERE", "Sim2Real Benchmark", "Gallery", "Reference"]
    for target in (
        "source/quickstart",
        "source/application-reference",
        "source/concepts",
        "source/installation",
        "source/application-migration",
    ):
        assert target in index

    benchmark_block = index.split(":caption: Sim2Real Benchmark", 1)[1].split(":caption: Gallery", 1)[0]
    assert [line.strip() for line in benchmark_block.splitlines() if "<benchmark/" in line] == [
        "segmentation <benchmark/segmentation>",
        "bin2cell <benchmark/bin2cell>",
        "batch <benchmark/batch>",
        "spot <benchmark/spot>",
        "imputation_and_dropout <benchmark/imputation_and_dropout>",
        "plot_imputation_case <benchmark/plot_imputation_case>",
    ]

    gallery_block = index.split(":caption: Gallery", 1)[1].split(":caption: Reference", 1)[0]
    assert [line.strip() for line in gallery_block.splitlines() if "<case/" in line] == [
        "VisiumHD (hST platform) sp-SVC <case/VisiumHD_sp_SVC>",
        "Xenium (iST platform) sc-SVC T cells <case/Xenium_sc_SVC_T>",
        "Xenium (iST platform) sc-SVC Fibroblast <case/Xenium_sc_SVC_Fibroblast>",
        "Xenium (iST platform) sc-SVC Mono/Macro <case/Xenium_sc_SVC_Monocyte>",
        "Visium (sST platform) sc-SVC mouse brain <case/Visium_sc_SVC_mouse_brain>",
    ]


def test_nblink_inventory_targets_the_current_notebooks():
    expected = {
        "docs/benchmark/segmentation.nblink": "reproduce/benchmark/segmentation.ipynb",
        "docs/benchmark/bin2cell.nblink": "reproduce/benchmark/bin2cell.ipynb",
        "docs/benchmark/batch.nblink": "reproduce/benchmark/batch.ipynb",
        "docs/benchmark/spot.nblink": "reproduce/benchmark/spot.ipynb",
        "docs/benchmark/imputation_and_dropout.nblink": "reproduce/benchmark/imputation_and_dropout.ipynb",
        "docs/benchmark/plot_imputation_case.nblink": "reproduce/benchmark/plot_imputation_case.ipynb",
        "docs/case/VisiumHD_sp_SVC.nblink": "reproduce/case/VisiumHD_sp_SVC.ipynb",
        "docs/case/Xenium_sc_SVC_T.nblink": "reproduce/case/Xenium_sc_SVC_T.ipynb",
        "docs/case/Xenium_sc_SVC_Fibroblast.nblink": "reproduce/case/Xenium_sc_SVC_Fibroblast.ipynb",
        "docs/case/Xenium_sc_SVC_Monocyte.nblink": "reproduce/case/Xenium_sc_SVC_Monocyte.ipynb",
        "docs/case/Visium_sc_SVC_mouse_brain.nblink": "reproduce/case/Visium_sc_SVC_mouse_brain.ipynb",
    }
    observed = {
        str(path.relative_to(ROOT)): path
        for directory in (ROOT / "docs/benchmark", ROOT / "docs/case")
        for path in directory.glob("*.nblink")
    }
    assert set(observed) == set(expected)
    for relative, target in expected.items():
        payload = json.loads(_read(ROOT / relative))
        assert payload == {"path": f"../../{target}"}
        assert (ROOT / target).is_file()


def test_visium_notebook_uses_the_current_name_and_entrypoint():
    old_path = ROOT / "reproduce/case/sc_SVC_sr_case_Visium_mouse_brain.ipynb"
    notebook_path = ROOT / "reproduce/case/Visium_sc_SVC_mouse_brain.ipynb"
    notebook = json.loads(_read(notebook_path))
    source = "\n".join(
        line for cell in notebook["cells"] for line in cell.get("source", [])
    )

    assert not old_path.exists()
    assert "python reconstruct.py --config configs/application/Visium.yaml" in source
    assert 'SAMPLE_NAME = "REVISEVisiumMouseBrain_sc-SVC"' in source
    assert '("route", "sc-SVC")' in source
    assert '("mode", "sr")' in source
    assert "p reconstruct.py" not in source
    assert "sc-SVC-sr" not in source


def test_reproduce_documentation_uses_the_requested_downloads_and_boundary():
    reproduce = _read(ROOT / "reproduce/README.md")
    readme = _read(ROOT / "README.md")
    urls = (
        "https://zenodo.org/records/21921802",
        "https://zenodo.org/records/18389835",
        "https://zenodo.org/records/18389211",
    )
    for url in urls:
        assert url in reproduce
        assert url in readme
    for label in (
        "Sim2Real-ST benchmark",
        "Reproduced benchmark results",
        "Real-world ST datasets",
        "Reproduced sp-SVC H5AD",
        "Reproduced sc-SVC H5AD",
    ):
        assert label in reproduce
        assert label in readme
    assert "python reproduce/benchmark_main.py" in reproduce
    assert "bash reproduce/benchmark_main.sh" in reproduce
    assert "--select-ct T" in reproduce
    assert (
        "[Application Reference](https://revise-svc.readthedocs.io/en/latest/source/application-reference.html)"
        in reproduce
    )
    assert "not evidence that the current source has been" in reproduce
    assert reproduce.index("VisiumHD sp-SVC") < reproduce.index("Xenium sc-SVC T cells") < reproduce.index("Visium sc-SVC mouse brain")


def test_installation_describes_supported_versions_and_optional_layers():
    installation = _read(ROOT / "docs/source/installation.rst")
    metadata = tomllib.loads(_read(ROOT / "pyproject.toml"))
    optional = metadata["project"]["optional-dependencies"]

    assert "Python 3.10 and 3.11" in installation
    assert optional["tacco"] == ["tacco==0.5.0"]
    for extra in set(optional) - {"dev"}:
        assert f'python -m pip install ".[{extra}]"' in installation
    assert "exactly three maintained files" in installation
    assert "Xenium.yaml" in installation
    assert "POT is the default OT implementation" in installation
    assert "TACCO as an alternative OT method" in installation


def test_docs_have_no_machine_paths_and_documentation_build_is_static():
    for path in (*USER_DOCS, ROOT / "docs/source/application-migration.rst"):
        assert not re.search(r"/Users/|/home/", _read(path))

    conf = _read(ROOT / "docs/conf.py")
    assert 'nbsphinx_execute = "never"' in conf
    assert '"nbsphinx_link"' in conf
    assert 'html_title = "REVISE documentation"' in conf
    assert '"plans/**"' in conf
    assert '"design/**"' in conf


def test_documentation_links_and_failure_contract_are_navigable_and_accurate():
    gallery = _read(ROOT / "docs/source/gallery.rst")
    architecture = " ".join(_read(ROOT / "docs/source/architecture.rst").split())

    assert "`reproduce/README.md`_" in gallery
    assert "https://github.com/wuys13/REVISE/blob/main/reproduce/README.md" in gallery
    assert "not reader-atomic or crash-atomic" in architecture
    assert "catchable replacement failures attempt rollback" in architecture


def test_public_api_diagram_uses_sc_svc_modes():
    diagram = _read(ROOT / "docs/source/api/classes_revise.svg")
    assert "sc-SVC-sr" not in diagram
    assert "sc-SVC (cluster or SR mode)" in diagram
    assert "H5AD result(s) + manifest" in diagram
