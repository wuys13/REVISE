"""Shared fixtures for the REVISE test suite."""

import pytest

from tests._paths import REPO_ROOT


@pytest.fixture
def repo_root():
    return REPO_ROOT

