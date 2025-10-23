from ...core import SetupBase, SetupCollectionBase
from dataclasses import dataclass


class ParametersSetup(SetupBase):

    initial_states: bool = NotImplemented
    drag_coefficient: bool = NotImplemented
    radiation_pressure_coefficient: bool = NotImplemented
