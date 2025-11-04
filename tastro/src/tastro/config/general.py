from .core import SetupBase

from .environment import EnvironmentSetup
from .propagation import PropagationSetup
from .estimation import EstimationSetup
from tudatpy.astro import time_representation as ttime
from pathlib import Path
import yaml
import typing
from ..logging import log

if typing.TYPE_CHECKING:
    from ..io.cli import CommandLineArguments


class SimulationIntervalSetup(SetupBase):

    initial_epoch: ttime.Time
    final_epoch: ttime.Time


class CaseSetup(SetupBase):

    time: SimulationIntervalSetup
    environment: EnvironmentSetup
    propagation: PropagationSetup
    estimation: EstimationSetup

    perform_estimation: bool = False
    perform_propagation: bool = False
    evaluate_accelerations: bool = False

    @classmethod
    def from_config_file(cls, config_path: Path) -> "CaseSetup":

        log.info(f"Loading configuration from {config_path}")

        raw_config = yaml.safe_load(config_path.open("r"))
        output = cls.from_raw(raw_config)

        log.info(f"Finished loading configuration")

        return output

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
