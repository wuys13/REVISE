from __future__ import annotations

import logging
import sys
import threading
from pathlib import Path
from typing import Optional


class Logger:
    _instances = {}
    _lock = threading.Lock()

    def __new__(cls, log_name: str = "GlobalLogger", log_file: str = "app.log", level=logging.INFO, **kwargs):
        # Keep compatibility with the original constructor keyword: `name=...`.
        if "name" in kwargs and (not log_name or log_name == "GlobalLogger"):
            log_name = str(kwargs["name"])
        name = log_name or "GlobalLogger"
        if name not in cls._instances:
            with cls._lock:
                if name not in cls._instances:
                    instance = super().__new__(cls)
                    instance._init_logger(log_name=name, log_file=log_file, level=level)
                    cls._instances[name] = instance
        return cls._instances[name]

    def _init_logger(self, log_name: str = "GlobalLogger", log_file: str = "app.log", level=logging.INFO) -> None:
        self.logger = logging.getLogger(log_name)
        self.logger.setLevel(level)

        if not self.logger.handlers:
            formatter = logging.Formatter(
                "[%(asctime)s] [%(levelname)s] [%(filename)s:%(lineno)d] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )

            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)

            file_handler = logging.FileHandler(log_file, encoding="utf-8")
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        return self.logger


def ensure_logger(logger: Optional[logging.Logger] = None) -> logging.Logger:
    if logger is not None:
        return logger
    return Logger().get_logger()


def build_run_logger(run_name: str, run_dir: Path) -> logging.Logger:
    run_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / "run.log"
    return Logger(log_name=run_name, log_file=str(log_file)).get_logger()


def log_exception_to_run_file(logger: logging.Logger, message: str) -> None:
    """Write the active exception traceback to run.log without echoing it."""
    record = logger.makeRecord(
        logger.name,
        logging.ERROR,
        __file__,
        0,
        message,
        (),
        sys.exc_info(),
    )
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.handle(record)
