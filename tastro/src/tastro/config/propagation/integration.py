from ..core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from tudatpy.dynamics.propagation_setup import (
    integrator as tigrs,
    propagator as tprops,
)
from tudatpy.astro import time_representation as ttime
import numpy as np


class GeneralIntegrationSetup(SetupBase):

    starting_point: str = NotImplemented
    custom_start_epoch: ttime.Time | None = NotImplemented
    terminate_exactly: bool = NotImplemented
    state_representation: tprops.TranslationalPropagatorType = NotImplemented
    integrator_type: str = NotImplemented


class RKFixedSetup(SetupBase):

    coefficients: tigrs.CoefficientSets = NotImplemented
    step_size: ttime.Time = NotImplemented


class RKVariableSetup(SetupBase):

    coefficients: tigrs.CoefficientSets = NotImplemented
    initial_step: ttime.Time = NotImplemented
    rtols: np.ndarray = NotImplemented
    atols: np.ndarray = NotImplemented
    safety_factor: float = NotImplemented
    min_increment: float = NotImplemented
    max_increment: float = NotImplemented
    min_step: ttime.Time = NotImplemented
    max_step: ttime.Time = NotImplemented


@dataclass
class IntegrationSetup(SetupCollectionBase):

    general: GeneralIntegrationSetup
    rkf_fixed: RKFixedSetup
    rkf_variable: RKVariableSetup
