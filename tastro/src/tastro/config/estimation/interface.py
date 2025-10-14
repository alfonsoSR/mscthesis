from ..core import SetupCollectionBase
from dataclasses import dataclass
from .observation_models import ObservationModelSetup
from .observations import ObservationsSetup
from .light_propagation import LightTimeSetup


@dataclass
class EstimationSetup(SetupCollectionBase):

    present: bool
    light_propagation: LightTimeSetup
    observation_models: ObservationModelSetup
    observations: ObservationsSetup
