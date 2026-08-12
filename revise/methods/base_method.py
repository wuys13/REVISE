from abc import ABC
from abc import abstractmethod

from revise.tools.log import ensure_logger


class BaseMethod(ABC):
    def __init__(self, config, logger=None):
        self.config = config
        self.logger = ensure_logger(logger)

    @abstractmethod
    def run(self, *args, **kwargs):
        raise NotImplementedError("run method not implemented")
