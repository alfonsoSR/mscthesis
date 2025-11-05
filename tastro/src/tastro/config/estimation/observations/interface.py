from ...core import SetupBase
from .doppler import ClosedLoopObservationsSetup, OpenLoopObservationsSetup
from .cartesian import CartesianObservationsSetup


class ObservationsSetup(SetupBase):

    closed_loop: ClosedLoopObservationsSetup
    cartesian: CartesianObservationsSetup
    open_loop: OpenLoopObservationsSetup
