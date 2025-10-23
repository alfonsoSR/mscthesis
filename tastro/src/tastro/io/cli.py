import numpy as np
from dataclasses import dataclass
from pathlib import Path
import argparse


@dataclass
class CommandLineArguments:

    config_file: Path
    propagate: bool
    estimate: bool
    accelerations: bool


class UserInputParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "configuration_dir", help="Path to configuration file"
        )
        self.add_argument("-p", dest="propagate", action="store_true")
        self.add_argument("-e", dest="estimate", action="store_true")
        self.add_argument("-a", dest="accelerations", action="store_true")

        return None

    def parse_args(self) -> CommandLineArguments:

        # Get default output parse_args
        defaults = super().parse_args()

        # Verify path to configuration file
        config_path = (
            Path(defaults.configuration_dir).absolute() / "configuration.yaml"
        )
        if not config_path.exists():
            raise FileNotFoundError(
                f"Invalid configuration directory: {defaults.configuration_file}"
            )

        # Package output in dataclass
        output = CommandLineArguments(
            config_file=config_path,
            propagate=defaults.propagate,
            estimate=defaults.estimate,
            accelerations=defaults.accelerations,
        )
        return output
