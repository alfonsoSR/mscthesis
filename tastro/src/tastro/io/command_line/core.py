import argparse
from ...core import AutoDataclass
from pathlib import Path
import abc
from ...logging import log
import traceback
from typing import Any, Type


class CommandLineInput(metaclass=AutoDataclass):

    source_dirs: list[Path]
    verbose: bool

    @property
    def source_dir(self) -> Path:

        return self.source_dirs[0]


class CommandLineParser[T: CommandLineInput](
    argparse.ArgumentParser, metaclass=abc.ABCMeta
):

    namespace: Type[T]

    def __init__(self) -> None:

        # Initialize with constructor of base class
        super().__init__()

        # Add common argument to all command line parsers
        self.add_argument(
            "source_dirs",
            nargs="+",
            help="Directories containing configuration and results",
        )
        self.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            action="store_true",
            help="Print debug information",
        )

        return None

    def _ensure_path_exists(self, path: Path, id: str) -> None:

        if not path.exists():

            log.fatal(f"Invalid {id}: {path}")
            log.fatal(f"{traceback.extract_stack()[-2]}")
            exit(1)

        return None

    @abc.abstractmethod
    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        # Process path to source directories
        arguments["source_dirs"] = []
        for source_dir in defaults.source_dirs:

            arguments["source_dirs"].append(Path(source_dir).absolute())
            self._ensure_path_exists(
                arguments["source_dirs"][-1], "source directory"
            )

        # Process verbosity level
        arguments["verbose"] = defaults.verbose

        return arguments

    def _default_arguments(self) -> argparse.Namespace:

        return super().parse_args()

    def parse_args(self) -> T:

        defaults = self._default_arguments()
        arguments = self.local_parser(defaults, {})

        return self.namespace(**arguments)


class CLParser(CommandLineParser[CommandLineInput]):

    namespace = CommandLineInput

    def local_parser(self, defaults, arguments) -> dict[str, Any]:
        return super().local_parser(defaults, arguments)


class CLInputFigure(CommandLineInput):

    base_dir: Path
    save: bool
    show: bool
    name_modifier: str


class CLParserFigureBase[T: CLInputFigure](CommandLineParser[T]):

    def __init__(self) -> None:

        super().__init__()

        self.figures_group = self.add_argument_group("Figure processing")
        self.figures_group.add_argument(
            "--base-dir",
            dest="base_dir",
            default=".",
            help="Base directory to define IDs for configurations",
        )
        self.figures_group.add_argument(
            "-s",
            dest="save",
            action="store_true",
            help="Save figure",
        )
        self.figures_group.add_argument(
            "-x",
            dest="show",
            action="store_false",
            help="Do not show figure",
        )
        self.figures_group.add_argument(
            "--name-mod",
            dest="name_modifier",
            default="",
            help="Modifier for the file name of the figure",
        )

        return None

    @abc.abstractmethod
    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)

        # Define path to base directory
        base_dir = Path(defaults.base_dir).absolute()
        self._ensure_path_exists(base_dir, "base directory")
        arguments["base_dir"] = base_dir

        # Update with name modifier
        arguments["name_modifier"] = defaults.name_modifier

        # Update with save and show options
        arguments["save"] = defaults.save
        arguments["show"] = defaults.show

        return arguments


class CLParserFigure(CLParserFigureBase[CLInputFigure]):

    namespace = CLInputFigure

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:
        return super().local_parser(defaults, arguments)
