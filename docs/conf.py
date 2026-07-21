"""Configuration file for the Sphinx documentation builder."""

from __future__ import annotations

import sys
from pathlib import Path

# -- Path setup --------------------------------------------------------------

# Use the repository root so that `import revise` works when building locally.
# The previous value climbed two directories up, escaping the repo and
# breaking imports from the current `revise` package.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DOCS_SOURCE = Path(__file__).resolve().parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# -- Project information -----------------------------------------------------

project = "revise-svc"
copyright = "Pending owner confirmation"
author = "Pending owner confirmation"
version_namespace: dict[str, str] = {}
exec((PROJECT_ROOT / "revise" / "_version.py").read_text(), version_namespace)
version = version_namespace["__version__"]
release = version


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "myst_parser",
]

templates_path = ["source/_templates"] if (DOCS_SOURCE / "source" / "_templates").exists() else []
exclude_patterns: list[str] = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    # These source-only migration records are absent from the clean export and
    # remain excluded while this legacy checkout is still used to build docs.
    "plans/**",
    "superpowers/**",
    # Ignore stale local copies from the pre-unified API layout. Current
    # autosummary pages are generated from source/api/index.rst during build.
    "source/api/generated/revise.application.*",
    "source/api/generated/revise.benchmark.*",
    "source/api/generated/revise.conf.*",
]

autodoc_mock_imports = [
    "anndata",
    "gseapy",
    "igraph",
    "leidenalg",
    "mkl",
    "networkx",
    "numba",
    "numpy",
    "ot",
    "pandas",
    "scanpy",
    "scipy",
    "seaborn",
    "skimage",
    "sklearn",
    "sparse_dot_mkl",
    "squidpy",
    "statsmodels",
    "torch",
    "tqdm",
]

# Sphinx's generic mock does not support assigning to matplotlib.rcParams in
# revise.analysis.bio. Give that one import the smallest compatible surface;
# autodoc handles every other unavailable scientific dependency itself.
from unittest.mock import MagicMock

mock_matplotlib = MagicMock()
mock_matplotlib.rcParams = {}
sys.modules.setdefault("matplotlib", mock_matplotlib)
sys.modules.setdefault("matplotlib.pyplot", MagicMock())

autosummary_generate = True
# Avoid dumping full Methods/Attributes tables for classes
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': False,
    'member-order': 'bysource',
}

napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["source/_static"] if (DOCS_SOURCE / "source" / "_static").exists() else []
html_css_files = ["revise.css"]
html_title = "REVISE documentation"

# Navigation UX tweaks for Read the Docs theme
html_theme_options = {
    # Keep the sidebar expanded; don't collapse on navigation
    "collapse_navigation": False,
    # Preserve sidebar scroll/position between pages
    "sticky_navigation": True,
    # Control depth of the sidebar tree
    "navigation_depth": 4,
    # Keep prev/next at bottom to reduce header movement
    "prev_next_buttons_location": "bottom",
    # Show only document titles in the sidebar (hides per-page headings)
    "titles_only": False,
}


# Set the root document to a top-level index, so that Sphinx
# generates an HTML `index.html` at the output root (required by RTD).
root_doc = "index"
