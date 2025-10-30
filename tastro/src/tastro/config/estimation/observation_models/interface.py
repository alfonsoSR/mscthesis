from .closed_loop import ClosedLoopDopplerSetup
from .cartesian import CartesianSetup
from ...core import SetupBase


class ObservationModelSetup(SetupBase):

    closed_loop: ClosedLoopDopplerSetup
    cartesian: CartesianSetup
