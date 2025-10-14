import numpy as np
from dataclasses import dataclass
from pathlib import Path
import argparse


@dataclass
class CommandLineArguments:

    config_file: Path


class UserInputParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "configuration_file", help="Path to configuration file"
        )

        return None

    def parse_args(self) -> CommandLineArguments:

        # Get default output parse_args
        defaults = super().parse_args()

        # Verify path to configuration file
        config_path = Path(defaults.configuration_file).absolute()
        if not config_path.exists():
            raise FileNotFoundError(
                f"Invalid configuration file: {defaults.configuration_file}"
            )

        # Package output in dataclass
        output = CommandLineArguments(config_file=config_path)
        return output
