from .closed_loop import ClosedLoopDopplerSetup
from .cartesian import CartesianSetup
from ...core import SetupCollectionBase
from dataclasses import dataclass


@dataclass
class ObservationModelSetup(SetupCollectionBase):

    closed_loop: ClosedLoopDopplerSetup
    cartesian: CartesianSetup
