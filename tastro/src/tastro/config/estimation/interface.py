from ..core import SetupCollectionBase
from dataclasses import dataclass
from .observation_models import ObservationModelSetup
from .observations import ObservationsSetup
from .light_propagation import LightTimeSetup
from .parameters import ParametersSetup


@dataclass
class EstimationSetup(SetupCollectionBase):

    light_propagation: LightTimeSetup
    observation_models: ObservationModelSetup
    observations: ObservationsSetup
    parameters: ParametersSetup
