from ...core import SetupBase
from .closed_loop import ClosedLoopObservationsSetup
from .cartesian import CartesianObservationsSetup


class ObservationsSetup(SetupBase):

    closed_loop: ClosedLoopObservationsSetup
    cartesian: CartesianObservationsSetup
