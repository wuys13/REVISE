from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - only for static analyzers
    from revise.application.sc_svc import ScSVC
    from revise.application.sc_svc_sr import ScSVCSr
    from revise.application.sp_svc import SpSVC
    from revise.application.sc_svc import ScSVCAnalysis

def __getattr__(name):
    """Lazily load application classes to avoid heavy imports at package load time."""
    mapping = {
        "SpSVC": "revise.application.sp_svc",
        "ScSVC": "revise.application.sc_svc",
        "ScSVCSr": "revise.application.sc_svc_sr",
        "ScSVCAnalysis": "revise.application.sc_svc"
    }
    if name in mapping:
        module = import_module(mapping[name])
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ['SpSVC', 'ScSVC', 'ScSVCSr', 'ScSVCAnalysis']
