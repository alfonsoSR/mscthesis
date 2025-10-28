import numpy as np
from dataclasses import dataclass
from pathlib import Path
import pickle
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


@dataclass
class PropagationOutput:

    epochs: np.ndarray
    cstate_j2000: np.ndarray
    rstate_j2000: np.ndarray
    cstate_rsw: np.ndarray
    rstate_rsw: np.ndarray
    dvars: dict[str, np.ndarray]
    number_of_function_evaluations: float
    ustate: np.ndarray

    def save_to_file(self, file: Path) -> None:

        with file.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @staticmethod
    def from_file(file: Path) -> "PropagationOutput":

        with file.open("rb") as buffer:
            _self = pickle.load(buffer)

        return _self

    @staticmethod
    def from_config_file(config_file: Path) -> "PropagationOutput":

        # Get path to results from path to config file
        results_file = config_file.parent / "results.pkl"
        return PropagationOutput.from_file(results_file)
