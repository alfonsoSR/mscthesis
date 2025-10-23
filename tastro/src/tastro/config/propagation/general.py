from ..core import SetupCollectionBase
from dataclasses import dataclass
from .accelerations import TargetAccelerationSetup
from .integration import IntegrationSetup


@dataclass
class PropagationSetup(SetupCollectionBase):

    integrator: IntegrationSetup
    accelerations: dict[str, TargetAccelerationSetup]
