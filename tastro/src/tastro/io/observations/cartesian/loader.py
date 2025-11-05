import numpy as np
import pickle
from ..common import RawObservationRecord
from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from tudatpy.astro import time_representation as ttime
    from ....config import CaseSetup


class CartesianObservationRecord(RawObservationRecord):

    @property
    def state_history(self) -> "dict[ttime.Time, np.ndarray]":

        return {
            epoch: obs for (epoch, obs) in zip(self.epochs, self.observations)
        }

    @staticmethod
    def from_config(config: "CaseSetup") -> "CartesianObservationRecord":

        source = config.estimation.observations.cartesian.sources[0]
        with source.open("rb") as buffer:
            state_history: dict[ttime.Time, np.ndarray] = pickle.load(buffer)

        epochs = np.array(list(state_history.keys()))
        observations = np.array(list(state_history.values()))

        return CartesianObservationRecord(epochs, observations)
