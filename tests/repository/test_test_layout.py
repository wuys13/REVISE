"""Test-suite ownership layout contract.

Covers: top-level test directories and per-file README discoverability.
Proof limit: validates organization, not the behavior covered by each test.
"""

from __future__ import annotations

from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[1]
OWNERS = {
    "analysis",
    "application",
    "backend",
    "benchmark",
    "config",
    "integration",
    "io",
    "preprocess",
    "recon",
    "repository",
}


def test_tests_are_grouped_by_production_owner():
    assert not list(TEST_ROOT.glob("test_*.py"))
    assert OWNERS <= {path.name for path in TEST_ROOT.iterdir() if path.is_dir()}


def test_readme_indexes_every_test_module():
    readme = (TEST_ROOT / "README.md").read_text(encoding="utf-8")
    test_modules = sorted(
        path.relative_to(TEST_ROOT).as_posix()
        for path in TEST_ROOT.rglob("test_*.py")
    )

    for module in test_modules:
        assert f"`{module}`" in readme

