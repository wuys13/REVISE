"""Backend runtime runners.

Keep this package lightweight: do not import heavy modules at package import time.
Import concrete runners from submodules directly, e.g.:
- revise.backend.runners.sp_svc_application.SpSVC
- revise.backend.runners.sc_svc_application.ScSVC
"""

__all__ = []
