from ...core import SetupCollectionBase
from .closed_loop import ClosedLoopObservationsSetup
from .cartesian import CartesianObservationsSetup
from dataclasses import dataclass


@dataclass
class ObservationsSetup(SetupCollectionBase):

    closed_loop: ClosedLoopObservationsSetup
    cartesian: CartesianObservationsSetup
