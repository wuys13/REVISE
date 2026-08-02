from __future__ import annotations

import importlib
import importlib.metadata


SUPPORTED_TACCO_VERSION = "0.5.0"
_INSTALL_HINT = (
    'run: python -m pip install "revise-svc[tacco]" '
    'or python -m pip install "tacco==0.5.0"; if a different OT algorithm '
    'is acceptable, explicitly select POT (application CLI: "--ot-method pot"); '
    "REVISE does not fall back automatically"
)


def require_tacco():
    """Load the one TACCO runtime version verified by REVISE."""
    try:
        module = importlib.import_module("tacco")
    except ModuleNotFoundError as exc:
        if exc.name != "tacco":
            raise ImportError(
                "TACCO OT could not import tacco=="
                f"{SUPPORTED_TACCO_VERSION} because runtime dependency "
                f"{exc.name!r} is missing; {_INSTALL_HINT}"
            ) from exc
        raise ModuleNotFoundError(
            f"TACCO OT requires tacco=={SUPPORTED_TACCO_VERSION}; {_INSTALL_HINT}",
            name="tacco",
        ) from exc
    except ImportError as exc:
        raise ImportError(
            f"TACCO OT could not import tacco=={SUPPORTED_TACCO_VERSION}; "
            f"{_INSTALL_HINT}"
        ) from exc

    try:
        version = importlib.metadata.version("tacco")
    except importlib.metadata.PackageNotFoundError:
        version = getattr(module, "__version__", "unknown")
    if version != SUPPORTED_TACCO_VERSION:
        raise RuntimeError(
            "REVISE requires tacco=="
            f"{SUPPORTED_TACCO_VERSION} for TACCO OT, found {version!r}; "
            f"{_INSTALL_HINT}"
        )
    return module
