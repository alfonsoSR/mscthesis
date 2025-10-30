from ..core import SetupBase

from .accelerations import TargetAccelerationSetup
from .integration import IntegrationSetup


class PropagationSetup(SetupBase):

    integrator: IntegrationSetup
    accelerations: dict[str, TargetAccelerationSetup]
