from __future__ import annotations

from abc import ABC
from abc import abstractmethod

from revise.utils.logging import ensure_logger


class BaseKernel(ABC):
    """Backend-side algorithm kernel base class."""

    def __init__(self, config, logger=None):
        self.config = config
        self.logger = ensure_logger(logger)

    @abstractmethod
    def run(self, *args, **kwargs):
        raise NotImplementedError("run method not implemented")
