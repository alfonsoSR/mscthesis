from tastro.config import CaseSetup
from ..common import RawObservationRecord
import numpy as np
from typing import TYPE_CHECKING
from tudatpy.astro import time_representation as ttime
from ....logging import log
from .fdet_parser import FdetsLoader

if TYPE_CHECKING:
    from ....config import CaseSetup


class OpenLoopDopplerObservationRecord(RawObservationRecord):

    @classmethod
    def from_observation_history(
        cls, observation_history: dict[ttime.Time, float]
    ) -> "OpenLoopDopplerObservationRecord":

        epochs = np.array(list(observation_history.keys()))
        observations = np.array(list(observation_history.values()))

        return OpenLoopDopplerObservationRecord(epochs, observations)


def load_open_loop_doppler_observations_from_config(
    config: "CaseSetup",
) -> "dict[str, OpenLoopDopplerObservationRecord]":

    # Display information
    log.info("Loading open-loop Doppler observations")

    # Initialize container of observations per station
    epochs_per_station: dict[str, list[ttime.Time]] = {}
    observations_per_station: dict[str, list[float]] = {}
    noise_per_station: dict[str, list[float]] = {}

    # Loop over sources
    for fdet in config.estimation.observations.open_loop.sources:

        # Display info
        log.debug(f"Fdet source: {fdet.path.name}")

        # Initialize IFMS loader and extract observations
        loader = FdetsLoader(fdet.path, fdet.station)
        epochs = loader.load_observation_epochs()
        observations = loader.load_observation_values()
        noise = loader.load_observation_noise()

        # Identify station and update container
        if fdet.station not in epochs_per_station:
            epochs_per_station[fdet.station] = epochs
            observations_per_station[fdet.station] = observations
            noise_per_station[fdet.station] = noise
        else:
            epochs_per_station[fdet.station] += epochs
            observations_per_station[fdet.station] += observations
            noise_per_station[fdet.station] += noise

    # Pack observations into data structures and return
    return {
        station: OpenLoopDopplerObservationRecord(
            epochs=np.array(epochs_per_station[station], dtype=ttime.Time),
            observations=np.array(station_observations, dtype=float),
            noise=np.array(noise_per_station[station], dtype=float),
        )
        for station, station_observations in observations_per_station.items()
    }
