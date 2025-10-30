from ..core import SetupBase

from .observation_models import ObservationModelSetup
from .observations import ObservationsSetup
from .light_propagation import LightTimeSetup
from .parameters import ParametersSetup


class EstimationSetup(SetupBase):

    light_propagation: LightTimeSetup
    observation_models: ObservationModelSetup
    observations: ObservationsSetup
    parameters: ParametersSetup
