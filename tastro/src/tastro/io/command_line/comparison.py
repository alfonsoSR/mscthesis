from typing import Any
from .core import (
    CommandLineInput,
    CommandLineParser,
    CLInputFigure,
    CLParserFigureBase,
)
from pathlib import Path


class CLInputComparison(CLInputFigure):

    ref_dir: Path


class CLParserComparison(CLParserFigureBase[CLInputComparison]):

    namespace = CLInputComparison

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "-r",
            dest="ref_dir",
            required=True,
            help="Directory with reference configuration and results",
        )

        return None

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)

        # Add reference directory to arguments
        ref_dir = Path(defaults.ref_dir).absolute()
        self._ensure_path_exists(ref_dir, "reference directory")
        arguments["ref_dir"] = ref_dir

        return arguments


class CLInputGeneralComparison(CLInputFigure):

    ref_dirs: list[Path]

    @classmethod
    def from_input_figure(
        cls, user_input: CLInputFigure
    ) -> "CLInputGeneralComparison":

        arguments = user_input.__dict__.copy()
        arguments["ref_dirs"] = []

        return cls(**arguments)


class CLParserGeneralComparison(CLParserFigureBase[CLInputGeneralComparison]):

    namespace = CLInputGeneralComparison

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "-r",
            dest="ref_dirs",
            nargs="+",
            required=True,
            help="List of directories with reference configurations",
        )

        return None

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)

        # Add reference directories to arguments
        ref_dirs: list[Path] = []
        for ref_dir in defaults.ref_dirs:
            ref_dirs.append(Path(ref_dir).absolute())
            self._ensure_path_exists(ref_dirs[-1], "reference directory")

        arguments["ref_dirs"] = ref_dirs
        return arguments
