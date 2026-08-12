import logging
import threading
from typing import Optional


class Logger:
    _instances = {}
    _lock = threading.Lock()  # Ensure thread safety

    def __new__(cls, log_name: str = 'GlobalLogger', log_file='app.log', level=logging.INFO, **kwargs):
        """
        Create or return a logger instance keyed by log_name.

        Calls with the same log_name share an instance; different names get independent loggers.
        """
        name = log_name or 'GlobalLogger'

        if name not in cls._instances:
            with cls._lock:
                if name not in cls._instances:
                    instance = super(Logger, cls).__new__(cls)
                    instance._init_logger(log_name=name, log_file=log_file, level=level)
                    cls._instances[name] = instance
        return cls._instances[name]

    def _init_logger(self, log_name='GlobalLogger', log_file='app.log', level=logging.INFO):
        self.logger = logging.getLogger(log_name)
        self.logger.setLevel(level)

        # If handler already exists, don't add duplicate
        if not self.logger.handlers:
            formatter = logging.Formatter(
                '[%(asctime)s] [%(levelname)s] [%(filename)s:%(lineno)d] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )

            # Console output
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)

            # File output
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        return self.logger


def ensure_logger(logger: Optional[logging.Logger] = None) -> logging.Logger:
    """
    Return the provided logger or fall back to the global singleton.

    This helper avoids scattering `if logger is None` checks throughout the
    codebase and ensures logging always has a valid sink.
    """
    if logger is not None:
        return logger
    return Logger().get_logger()
