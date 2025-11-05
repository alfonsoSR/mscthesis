from tastro.config import CaseSetup
from ..common import RawObservationRecord
import numpy as np
from typing import TYPE_CHECKING, Self
from tudatpy.astro import time_representation as ttime
from ....logging import log
from .ifms_parser import IFMSLoader
from .odf_parser import ODFLoader
from .... import utils
from tudatpy import data as tdata

if TYPE_CHECKING:
    from ....config import CaseSetup


class ClosedLoopDopplerObservationRecord(RawObservationRecord):

    @classmethod
    def from_observation_history(
        cls, observation_history: dict[ttime.Time, float]
    ) -> "ClosedLoopDopplerObservationRecord":

        epochs = np.array(list(observation_history.keys()))
        observations = np.array(list(observation_history.values()))

        return ClosedLoopDopplerObservationRecord(epochs, observations)


def load_closed_loop_doppler_observations_from_config(
    config: "CaseSetup",
) -> "dict[str, ClosedLoopDopplerObservationRecord]":

    # Display information
    log.info("Loading closed-loop Doppler observations")

    # Initialize container of observations per station
    observations_per_station: dict[str, dict[ttime.Time, float]] = {}

    # Loop over IFMS files
    for ifms in config.estimation.observations.closed_loop.sources.ifms:

        # Display info
        log.debug(f"IFMS source: {ifms.name}")

        # Initialize IFMS loader and extract observations
        loader = IFMSLoader(ifms)
        observations = loader.load_raw_observations()

        # Identify station and update container
        station = utils.identify_station_from_id(loader.station_id)
        if station not in observations_per_station:
            observations_per_station[station] = observations
        else:
            observations_per_station[station].update(observations)

    # Loop over ODF files
    for odf_data in config.estimation.observations.closed_loop.sources.odf:

        # Display information
        log.debug(f"ODF source: {odf_data.path.name} - {odf_data.station}")

        # Initialize ODF loader and extract observations
        loader = ODFLoader(
            source_file=odf_data.path,
            rx_station=odf_data.station,
            observable_type="ESTRACK",
        )
        observations = loader.load_raw_observations(
            data_type=tdata.OdfDataType.two_way_doppler,
            time_range=None,
        )

        # Do not update container if no observations were loaded
        if len(observations) == 0:
            continue

        # Update container
        if odf_data.station not in observations_per_station:
            observations_per_station[odf_data.station] = observations
        else:
            observations_per_station[odf_data.station].update(observations)

    # Pack observations into data structures and return
    return {
        station: ClosedLoopDopplerObservationRecord.from_observation_history(
            history
        )
        for station, history in observations_per_station.items()
    }
