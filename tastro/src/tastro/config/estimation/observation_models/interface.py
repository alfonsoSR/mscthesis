from .doppler import ClosedLoopDopplerSetup, OpenLoopDopplerSetup
from .cartesian import CartesianSetup
from ...core import SetupBase


class ObservationModelSetup(SetupBase):

    closed_loop: ClosedLoopDopplerSetup
    cartesian: CartesianSetup
    open_loop: OpenLoopDopplerSetup
