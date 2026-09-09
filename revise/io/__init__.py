"""Input normalization services for REVISE.

The reconstruction kernels continue to consume AnnData objects. This package
keeps external input handling isolated so new spatial file ecosystems can be
added without changing algorithm code.
"""

from revise.io.input_service import REVISEInputService

__all__ = ["REVISEInputService"]
