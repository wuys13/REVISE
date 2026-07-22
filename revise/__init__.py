from revise._version import __version__

__all__ = ["REVISEPipeline", "__version__"]


def __getattr__(name):
    if name == "REVISEPipeline":
        from revise.framework import REVISEPipeline

        return REVISEPipeline
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
