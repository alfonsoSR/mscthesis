from ...config import CaseSetup
import numpy as np
import pickle
from tudatpy.astro import time_representation as ttime


class CartesianObservationRecord:

    def __init__(self, epochs: np.ndarray, observations: np.ndarray) -> None:

        self.epochs = epochs
        self.observations = observations
        self.state_history = {
            epoch: obs for (epoch, obs) in zip(self.epochs, self.observations)
        }

        return None

    @classmethod
    def from_config(cls, config: CaseSetup) -> "CartesianObservationRecord":

        source = config.estimation.observations.cartesian.state_history
        with source.open("rb") as buffer:
            state_history: dict[ttime.Time, np.ndarray] = pickle.load(buffer)

        epochs = np.array(list(state_history.keys()))
        observations = np.array(list(state_history.values()))

        return CartesianObservationRecord(epochs, observations)
