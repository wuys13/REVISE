from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session", autouse=True)
def assert_candidate_wheel_import():
    expected_prefix = os.environ.get("REVISE_EXPECT_INSTALLED_PREFIX")
    if expected_prefix is None:
        return

    import revise

    installed = Path(revise.__file__).resolve()
    prefix = Path(expected_prefix).resolve()
    assert installed.is_relative_to(prefix), (
        f"solver smoke imported {installed}, expected candidate under {prefix}"
    )

    workspace = os.environ.get("GITHUB_WORKSPACE")
    if workspace is not None:
        source = Path(workspace).resolve()
        assert not installed.is_relative_to(source), (
            f"solver smoke imported source checkout {installed}"
        )
