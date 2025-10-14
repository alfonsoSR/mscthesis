from .core import SetupCollectionBase, SetupBase
from dataclasses import dataclass
from .environment import EnvironmentSetup
from .propagation import PropagationSetup
from .estimation import EstimationSetup
from tudatpy.astro import time_representation as ttime
from pathlib import Path
import yaml


class SimulationIntervalSetup(SetupBase):

    initial_epoch: ttime.Time = NotImplemented
    final_epoch: ttime.Time = NotImplemented


@dataclass
class CaseSetup(SetupCollectionBase):

    time: SimulationIntervalSetup
    environment: EnvironmentSetup
    propagation: PropagationSetup
    estimation: EstimationSetup

    @classmethod
    def from_config_file(cls, config_path: Path) -> "CaseSetup":

        raw_config = yaml.safe_load(config_path.open("r"))
        return cls.from_raw(raw_config)
