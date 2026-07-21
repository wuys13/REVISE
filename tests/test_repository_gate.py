from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / "tools" / "check_repository.py"


def _write(root: Path, path: str, content: str = "fixture\n") -> None:
    target = root / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content, encoding="utf-8")


def _clean_fixture(root: Path) -> None:
    for path in [
        ".github/workflows/ci.yml",
        ".gitignore",
        "LICENSE",
        "MANIFEST.in",
        "MIGRATION.md",
        "README.md",
        "constraints/python-3.11.txt",
        "docs/index.rst",
        "pyproject.toml",
        "revise/__init__.py",
        "tests/test_cli_contract.py",
        "tools/check_repository.py",
    ]:
        _write(root, path)
    (root / "legacy-assets.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "source": {
                    "repository_name": "REVISE-legacy",
                    "repository": "https://github.com/wuys13/REVISE-legacy",
                    "commit": "a" * 40,
                },
                "excluded_assets": [],
            }
        ),
        encoding="utf-8",
    )


def _commit(root: Path, message: str = "clean root") -> None:
    if not (root / ".git").exists():
        subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    subprocess.run(["git", "add", "."], cwd=root, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=Repository Gate Test",
            "-c",
            "user.email=gate@example.invalid",
            "commit",
            "-qm",
            message,
        ],
        cwd=root,
        check=True,
    )


def _run(root: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, CHECKER, "--root", root],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_clean_repository_fixture_passes() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        _clean_fixture(root)
        _commit(root)
        result = _run(root)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "path",
    [
        "case.ipynb",
        "reproduce/case/output.txt",
        "docs/plans/internal.md",
        "docs/superpowers/specs/internal.md",
        "release/legacy/index.json",
        "release/0.1.0rc1/status.json",
        "build/generated.py",
    ],
)
def test_clean_repository_rejects_forbidden_tracked_paths(
    tmp_path: Path, path: str
) -> None:
    _clean_fixture(tmp_path)
    _write(tmp_path, path)
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "forbidden tracked path" in result.stderr


def test_clean_repository_rejects_an_unapproved_large_file(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    large = tmp_path / "large.bin"
    large.write_bytes(b"x" * (5 * 1024 * 1024 + 1))
    _commit(tmp_path)
    large.write_bytes(b"x\n")

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "tracked file exceeds 5 MiB" in result.stderr


def test_clean_repository_rejects_every_symlink(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    link = tmp_path / "docs/source-link.rst"
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to("../index.rst")
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "tracked symlink is forbidden" in result.stderr


def test_clean_repository_rejects_case_insensitive_notebook_suffix(
    tmp_path: Path,
) -> None:
    _clean_fixture(tmp_path)
    _write(tmp_path, "docs/source/HISTORICAL.IPYNB")
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "forbidden tracked path: docs/source/HISTORICAL.IPYNB" in result.stderr


def test_clean_repository_requires_valid_migration_identity(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    payload = json.loads((tmp_path / "legacy-assets.json").read_text())
    payload["source"]["repository_name"] = "REVISE"
    (tmp_path / "legacy-assets.json").write_text(json.dumps(payload))
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "legacy asset index is invalid" in result.stderr


def test_clean_repository_rejects_a_zenodo_migration_plan(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    (tmp_path / "MIGRATION.md").write_text(
        "Upload historical notebooks to Zenodo.\n", encoding="utf-8"
    )
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "Zenodo migration is not part of the clean repository" in result.stderr


def test_clean_repository_allows_a_zenodo_dataset_link(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    (tmp_path / "README.md").write_text(
        "Dataset: https://zenodo.org/records/17705737\n", encoding="utf-8"
    )
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode == 0, result.stderr


def test_clean_repository_allows_curated_notebooks(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    _write(tmp_path, "application_sc_SVC_analysis_case.ipynb", "{}\n")
    _write(tmp_path, "reproduce/benchmark/seg_benchmark.ipynb", "{}\n")
    _write(tmp_path, "reproduce/case/sp_SVC_case.ipynb", "{}\n")
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode == 0, result.stderr


def test_clean_repository_allows_readme_logo_assets(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    logo = tmp_path / "logo" / "REVISE.png"
    logo.parent.mkdir(parents=True)
    logo.write_bytes(b"png fixture\n")
    _commit(tmp_path)

    result = _run(tmp_path)

    assert result.returncode == 0, result.stderr


def test_clean_repository_allows_followup_commits(tmp_path: Path) -> None:
    _clean_fixture(tmp_path)
    _commit(tmp_path)
    (tmp_path / "README.md").write_text("second commit\n", encoding="utf-8")
    _commit(tmp_path, "second commit")

    result = _run(tmp_path)

    assert result.returncode == 0, result.stderr


def test_clean_repository_rejects_objects_reachable_only_from_an_extra_ref(
    tmp_path: Path,
) -> None:
    _clean_fixture(tmp_path)
    _commit(tmp_path)
    tree = subprocess.run(
        ["git", "rev-parse", "HEAD^{tree}"],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    extra = subprocess.run(
        ["git", "commit-tree", tree, "-m", "extra root"],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=True,
        env={
            **os.environ,
            "GIT_AUTHOR_NAME": "Repository Gate Test",
            "GIT_AUTHOR_EMAIL": "gate@example.invalid",
            "GIT_COMMITTER_NAME": "Repository Gate Test",
            "GIT_COMMITTER_EMAIL": "gate@example.invalid",
        },
    ).stdout.strip()
    subprocess.run(
        ["git", "update-ref", "refs/heads/extra", extra],
        cwd=tmp_path,
        check=True,
    )

    result = _run(tmp_path)

    assert result.returncode != 0
    assert "objects outside the HEAD closure" in result.stderr


def test_clean_repository_rejects_a_shallow_nonroot_checkout(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    _clean_fixture(source)
    _commit(source)
    (source / "README.md").write_text("second commit\n", encoding="utf-8")
    _commit(source, "second commit")
    shallow = tmp_path / "shallow"
    subprocess.run(
        ["git", "clone", "-q", "--depth=1", f"file://{source}", shallow],
        check=True,
    )

    result = _run(shallow)

    assert result.returncode != 0
    assert "shallow repositories are forbidden" in result.stderr
