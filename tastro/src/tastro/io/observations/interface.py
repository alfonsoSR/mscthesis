from .doppler import (
    DopplerObservationRecord,
    RawDopplerObservationRecord,
    get_metadata_from_ifms_filename,
    identify_station_from_id,
)
from ...config import CaseSetup
from tudatpy import data as tdata
import numpy as np
from ...logging import log


def load_doppler_observations_from_config(
    config: CaseSetup,
) -> dict[str, DopplerObservationRecord]:

    observations_per_station: dict[str, list[RawDopplerObservationRecord]] = {}
    log.info("Loading closed-loop Doppler observations")

    # Load raw data from IFMS files
    for ifms in config.estimation.observations.closed_loop.sources.ifms:

        log.debug(f"IFMS source: {ifms.name}")

        # Identify station from IFMS file name
        metadata = get_metadata_from_ifms_filename(ifms)
        station = identify_station_from_id(int(metadata["station_id"]))

        # Load observations for station
        observations = RawDopplerObservationRecord.from_ifms_file(ifms)

        # Add observations to container
        if station not in observations_per_station:
            observations_per_station[station] = [observations]
        else:
            observations_per_station[station].append(observations)

    # Load raw data from ODF files
    for odf_data in config.estimation.observations.closed_loop.sources.odf:

        log.debug(f"ODF source: {odf_data.path.name} - {odf_data.station}")

        # Load raw observations from ODF
        observations = RawDopplerObservationRecord.from_odf_file(
            odf_data.path,
            odf_data.station,
            data_type=tdata.OdfDataType.two_way_doppler,
        )

        # Add to container if not empty
        if len(observations) != 0:
            if odf_data.station not in observations_per_station:
                observations_per_station[odf_data.station] = [observations]
            else:
                observations_per_station[odf_data.station].append(observations)

    # Blend all the observations of the stations
    output: dict[str, DopplerObservationRecord] = {}
    for station, station_data in observations_per_station.items():

        epochs = np.concatenate(
            [item.epochs for item in station_data]
        ).flatten()
        values = np.concatenate(
            [item.observations for item in station_data]
        ).flatten()

        output[station] = DopplerObservationRecord.from_config(
            epochs, values, station, config
        )

    log.info("Finished loading raw Doppler observations")

    return output
