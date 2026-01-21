import logging
from pathlib import Path
import sys

GREY = "\x1b[38;20m"
BLUE = "\x1b[34;20m"
GREEN = "\x1b[32;20m"
YELLOW = "\x1b[33;20m"
RED = "\x1b[31;20m"
BOLD_RED = "\x1b[31;1m"
RESET = "\x1b[0m"
FORMAT = "%(levelname)-8s :: %(asctime)s :: %(message)s"
EXCEPTION_FORMAT = "\n"


class Formatter(logging.Formatter):

    _FORMATS: dict[int, str] | None = None

    @property
    def FORMATS(self) -> dict[int, str]:

        if self._FORMATS is None:
            raise NotImplementedError(
                f"Formats not set for {type(self).__name__}"
            )

        return self._FORMATS

    def format(self, record) -> str:
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


class FileFormatter(Formatter):

    _FORMATS = {
        logging.DEBUG: FORMAT,
        logging.INFO: FORMAT,
        logging.WARNING: FORMAT,
        logging.ERROR: FORMAT,
        logging.CRITICAL: FORMAT,
    }


class StdoutFormatter(Formatter):

    _FORMATS = {
        logging.DEBUG: GREEN + FORMAT + RESET,
        logging.INFO: BLUE + FORMAT + RESET,
        logging.WARNING: YELLOW + FORMAT + RESET,
        logging.ERROR: RED + FORMAT + RESET,
        logging.CRITICAL: BOLD_RED + FORMAT + RESET,
    }


class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    blue = "\x1b[34;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    _format = "%(levelname)-8s :: %(asctime)s :: %(message)s"

    FORMATS = {
        logging.DEBUG: green + _format + reset,
        logging.INFO: blue + _format + reset,
        logging.WARNING: yellow + _format + reset,
        logging.ERROR: red + _format + reset,
        logging.CRITICAL: bold_red + _format + reset,
    }

    def format(self, record) -> str:
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


class FileHandler(logging.FileHandler):

    @staticmethod
    def from_path_and_level(filename: Path, level: str) -> "FileHandler":

        handler = FileHandler(str(filename), mode="w")
        handler.setLevel(level)
        handler.setFormatter(FileFormatter())

        return handler


class AlternativeFileHandler(logging.StreamHandler):

    @staticmethod
    def from_buffer_and_level(buffer, level: str) -> "AlternativeFileHandler":

        handler = AlternativeFileHandler(buffer)
        handler.setLevel(level)
        handler.setFormatter(FileFormatter())

        return handler


class StdoutHandler(logging.StreamHandler):

    @staticmethod
    def from_level(level: str) -> "StdoutHandler":

        handler = StdoutHandler(stream=sys.stdout)
        handler.setLevel(level)
        handler.setFormatter(StdoutFormatter())

        return handler


log = logging.getLogger(__package__)
log.setLevel(logging.DEBUG)

# stdout_handler = logging.StreamHandler()
# stdout_handler.setLevel(logging.DEBUG)
# stdout_handler.setFormatter(StdoutFormatter())
# # log.addHandler(stdout_handler)

# file_handler = logging.FileHandler("output.log", mode="w")
# file_handler.setLevel(logging.DEBUG)
# file_handler.setFormatter(FileFormatter())
# # log.addHandler(file_handler)
