from .core import SetupCollectionBase, SetupBase
from dataclasses import dataclass
from .environment import EnvironmentSetup
from .propagation import PropagationSetup
from .estimation import EstimationSetup
from tudatpy.astro import time_representation as ttime
from pathlib import Path
import yaml
import typing

if typing.TYPE_CHECKING:
    from ..io.cli import CommandLineArguments


class SimulationIntervalSetup(SetupBase):

    initial_epoch: ttime.Time = NotImplemented
    final_epoch: ttime.Time = NotImplemented


@dataclass
class CaseSetup(SetupCollectionBase):

    time: SimulationIntervalSetup
    environment: EnvironmentSetup
    propagation: PropagationSetup
    estimation: EstimationSetup

    perform_estimation: bool = False
    perform_propagation: bool = False
    evaluate_accelerations: bool = False

    @classmethod
    def from_config_file(cls, config_path: Path) -> "CaseSetup":

        raw_config = yaml.safe_load(config_path.open("r"))
        return cls.from_raw(raw_config)

    @classmethod
    def from_user_input(cls, user_input: "CommandLineArguments") -> "CaseSetup":

        # Get setup from configuration file
        setup = cls.from_config_file(user_input.config_file)

        # Modify with command line input
        if user_input.propagate:
            setup.perform_propagation = True
        if user_input.estimate:
            setup.perform_estimation = True
        if user_input.accelerations:
            setup.evaluate_accelerations = True

        return setup
