"""Configuration file for the Sphinx documentation builder."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pypandoc

# -- Path setup --------------------------------------------------------------

# Use the repository root so that `import revise` works when building locally.
# The previous value climbed two directories up, escaping the repo and
# breaking imports from the current `revise` package.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DOCS_SOURCE = Path(__file__).resolve().parent

# nbsphinx invokes `pandoc` as a subprocess. The docs-only dependency bundles
# that executable for clean CI and Read the Docs environments.
PANDOC_DIR = Path(pypandoc.get_pandoc_path()).parent
os.environ["PATH"] = f"{PANDOC_DIR}{os.pathsep}{os.environ['PATH']}"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(DOCS_SOURCE) not in sys.path:
    sys.path.insert(0, str(DOCS_SOURCE))


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
    "nbsphinx",
    "nbsphinx_link",
]

# Reproduction notebooks are historical snapshots. Keep documentation builds
# deterministic and link every notebook page to the current workflow.
nbsphinx_execute = "never"
nbsphinx_prolog = r"""
.. note::

   **Historical static snapshot.** This notebook is preserved for its
   paper-facing analysis and is not executed during documentation builds. See
   :doc:`Quick Start </source/quickstart>` for the current workflow and input
   contract.
"""

# Keep Gallery thumbnails tied to representative embedded notebook outputs.
# Wide platform-comparison figures are focused on their rightmost REVISE panel
# by the Gallery CSS.
nbsphinx_thumbnails = {
    "case/VisiumHD_sp_SVC": "_images/case_VisiumHD_sp_SVC_36_6.png",
    "case/SlideSeq_mouse_olfactory_bulb_sp_SVC": (
        "_images/case_SlideSeq_mouse_olfactory_bulb_sp_SVC_8_0.png"
    ),
    "case/SlideSeq_mouse_colon_sp_SVC": (
        "_images/case_SlideSeq_mouse_colon_sp_SVC_9_0.png"
    ),
    "case/StereoSeq_zebrafish_5hpf_sp_SVC": (
        "_images/case_StereoSeq_zebrafish_5hpf_sp_SVC_8_0.png"
    ),
    "case/CosMx_SMI_267T_not_sp_SVC": (
        "_images/case_CosMx_SMI_267T_not_sp_SVC_5_0.png"
    ),
    "case/Xenium_sc_SVC_T": "_images/case_Xenium_sc_SVC_T_9_4.png",
    "case/Xenium_sc_SVC_Fibroblast": (
        "_images/case_Xenium_sc_SVC_Fibroblast_16_4.png"
    ),
    "case/Xenium_sc_SVC_Monocyte": (
        "_images/case_Xenium_sc_SVC_Monocyte_9_4.png"
    ),
    "case/osmFISH_sc_SVC_cluster": (
        "_images/case_osmFISH_sc_SVC_cluster_10_0.png"
    ),
    "case/MERFISH_Allen_VISp_sc_SVC_cluster": (
        "_images/case_MERFISH_Allen_VISp_sc_SVC_cluster_8_0.png"
    ),
    "case/Visium_sc_SVC_mouse_brain": (
        "_images/case_Visium_sc_SVC_mouse_brain_29_2.png"
    ),
}

templates_path = ["source/_templates"] if (DOCS_SOURCE / "source" / "_templates").exists() else []
exclude_patterns: list[str] = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    # Historical decision records stay in the repository but are not part of
    # the current user documentation build.
    "plans/**",
    "design/**",
    # Ignore stale local copies from the old internal API reference. Current
    # public autosummary pages are generated from source/api/index.rst.
    "source/api/generated/revise.application.*",
    "source/api/generated/revise.analysis.*Service.rst",
    "source/api/generated/revise.backend.*",
    "source/api/generated/revise.config.*",
    "source/api/generated/revise.recon.*",
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
html_static_path = ["source/_static", "../logo"]
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
    # Keep the sidebar focused on runnable pages and notebook entries.
    "titles_only": True,
}


# Set the root document to a top-level index, so that Sphinx
# generates an HTML `index.html` at the output root (required by RTD).
root_doc = "index"
