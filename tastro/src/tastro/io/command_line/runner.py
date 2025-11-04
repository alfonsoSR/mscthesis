from .core import CommandLineInput, CommandLineParser
from pathlib import Path
import typing


class CommandLineInputRunner(CommandLineInput):

    config_files: list[Path]
    metakernels: list[Path]
    propagate: bool
    estimate: bool

    @property
    def config_file(self) -> Path:
        return self.config_files[0]

    @property
    def metakernel(self) -> Path:
        return self.metakernels[0]


class CommandLineParserRunner(CommandLineParser[CommandLineInputRunner]):

    namespace = CommandLineInputRunner

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "-p",
            dest="propagate",
            action="store_true",
            help="Propagate translational dynamics",
        )
        self.add_argument(
            "-e",
            dest="estimate",
            action="store_true",
            help="Perform estimation based on configuration",
        )

        return None

    def local_parser(self, defaults, arguments) -> dict[str, typing.Any]:

        arguments = super().local_parser(defaults, arguments)

        # Define path to configuration file
        arguments["config_files"] = []
        arguments["metakernels"] = []
        for source_dir in arguments["source_dirs"]:

            arguments["config_files"].append(source_dir / "configuration.yaml")
            self._ensure_path_exists(
                arguments["config_files"][-1], "configuration file"
            )

            arguments["metakernels"].append(source_dir / "metak.tm")
            self._ensure_path_exists(arguments["metakernels"][-1], "metakernel")

        # config_file = arguments["source_dirs"][0] / "configuration.yaml"
        # self._ensure_path_exists(config_file, "configuration file")
        # arguments["config_file"] = config_file

        # # Define path to metakernel
        # metakernel = arguments["source_dirs"][0] / "metak.tm"
        # self._ensure_path_exists(metakernel, "metakernel")
        # arguments["metakernel"] = metakernel

        # Define flags: propagate and estimate
        arguments["propagate"] = defaults.propagate
        arguments["estimate"] = defaults.estimate

        return arguments

    # def parse_args(self):

    #     # Get default output of parse
    #     defaults = self._default_arguments()

    #     # Get arguments as dictionary
    #     arguments = self.local_parser(defaults, {})

    #     # Pack output in dataclass
    #     return CommandLineInputRunner(**arguments)
