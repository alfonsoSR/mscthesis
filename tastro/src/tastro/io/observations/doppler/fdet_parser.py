from pathlib import Path
from ....logging import log
import traceback
from tudatpy.astro import time_representation as ttime
from tudatpy import data as tdata
import numpy as np


class FdetsLoader:

    def __init__(self, source_file: Path, station: str) -> None:

        # Ensure file exists
        if not source_file.exists():
            log.fatal(f"Invalid Fdet file: {source_file}")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Load raw data
        columns = [
            "utc_datetime_string",
            "signal_to_noise_ratio",
            "normalised_spectral_max",
            "doppler_measured_frequency_hz",
            "doppler_noise_hz",
        ]
        raw_data = tdata.read_fdets_file(str(source_file), columns)
        self.data = raw_data.raw_datamap

        # Save station name
        self.station = station

        # Initialize time converter
        self.time_converter = ttime.default_time_scale_converter()

        return None

    def load_observation_epochs(self) -> list[ttime.Time]:

        # Get observation epochs in TDB
        log.warning(
            "Ground station position not considered in UTC-TDB of Fdets"
        )
        observation_epochs = [
            self.time_converter.convert_time_object(
                input_scale=ttime.TimeScales.utc_scale,
                output_scale=ttime.TimeScales.tdb_scale,
                input_value=ttime.DateTime.from_iso_string(
                    epoch
                ).to_epoch_time_object(),
            )
            for epoch in self.data["utc_datetime_string"]
        ]
        return observation_epochs

    def load_observation_values(self) -> list[float]:

        return [
            8.412e9 + float(value)
            for value in self.data["doppler_measured_frequency_hz"]
        ]

    def load_observation_noise(self) -> list[float]:

        return [float(value) for value in self.data["doppler_noise_hz"]]

    def load_raw_observations(self) -> dict[ttime.Time, float]:

        # Get observation epochs in TDB
        observation_epochs = self.load_observation_epochs()

        # Return observation history
        log.warning("The reference frequency of the Fdets is hard coded")
        return {
            epoch: value
            for (epoch, value) in zip(
                observation_epochs, self.load_observation_values()
            )
        }
