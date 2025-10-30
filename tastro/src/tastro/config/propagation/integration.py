from ..core import SetupBase

from tudatpy.dynamics.propagation_setup import (
    integrator as tigrs,
    propagator as tprops,
)
from tudatpy.astro import time_representation as ttime
import numpy as np


class GeneralIntegrationSetup(SetupBase):

    starting_point: str
    custom_start_epoch: ttime.Time
    terminate_exactly: bool
    state_representation: tprops.TranslationalPropagatorType
    integrator_type: str


class RKFixedSetup(SetupBase):

    present: bool
    coefficients: tigrs.CoefficientSets
    step_size: ttime.Time


class RKVariableSetup(SetupBase):

    present: bool
    coefficients: tigrs.CoefficientSets
    initial_step: ttime.Time
    rtols: np.ndarray
    atols: np.ndarray
    safety_factor: float
    min_increment: float
    max_increment: float
    min_step: ttime.Time
    max_step: ttime.Time


class IntegrationSetup(SetupBase):

    general: GeneralIntegrationSetup
    rkf_fixed: RKFixedSetup
    rkf_variable: RKVariableSetup
