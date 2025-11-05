from pathlib import Path
from ....logging import log
import traceback
from tudatpy import data as tdata
from tudatpy.astro import time_representation as ttime
from .... import utils
from tudatpy.estimation.observations_setup import ancillary_settings as tancs
from typing import Literal


class ODFLoader:

    def __init__(
        self,
        source_file: Path,
        rx_station: str,
        observable_type: Literal["ESTRACK", "DSN"],
    ) -> None:

        # Ensure ODF file exists
        if not source_file.exists():
            log.fatal(f"Invalid ODF file: {source_file}")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Load raw data map
        self.data = tdata.read_odf_file(str(source_file))

        # Save station name and get station ID
        self.station = rx_station
        self.station_id = utils.get_id_from_station_name(rx_station)

        # Initialize time converter
        self.time_converter = ttime.default_time_scale_converter()

        # Define sign for computation of observables
        match observable_type:

            case "ESTRACK":
                self.observable_sign = -1.0
            case "DSN":
                self.observable_sign = 1.0
            case _:
                log.fatal(f"Invalid observable type: {observable_type}")
                log.fatal(traceback.extract_stack()[-2])
                exit(1)

        return None

    def get_reference_epoch(self) -> "ttime.Time":

        # Fail if not EME50
        if (
            self.data.file_reference_date != 19500101
            or self.data.file_reference_time != 0
        ):
            raise NotImplementedError("ODF reference epoch is not EME50")

        # Return reference epoch as time object
        return ttime.DateTime.from_iso_string(
            "1950-01-01T00:00:00"
        ).to_epoch_time_object()

    def get_turnaround_ratio_for_odf_block(
        self, common: "tdata.OdfCommonDataBlock"
    ) -> float:

        # Identify uplink band
        if common.uplink_band_id == 0:
            uplink_band = tancs.FrequencyBands.ku_band
        else:
            uplink_band = tancs.FrequencyBands(common.uplink_band_id - 1)

        # Identify downlink band
        if common.downlink_band_id == 0:
            downlink_band = tancs.FrequencyBands.ku_band
        else:
            downlink_band = tancs.FrequencyBands(common.downlink_band_id - 1)

        # Return default turnaround ratio based on bands
        return tancs.dsn_default_turnaround_ratios(
            uplink_band=uplink_band, downlink_band=downlink_band
        )

    def load_raw_observations(
        self,
        data_type: "tdata.OdfDataType",
        time_range: "tuple[ttime.Time, ttime.Time] | None" = None,
    ) -> dict[ttime.Time, float]:

        # Retrieve reference epoch from header
        reference_epoch = self.get_reference_epoch()

        # Get Earth-fixed position of ground station
        try:
            station_cstate = utils.get_ground_station_reference_state_itrf(
                station=self.station, source="from_glo"
            )
        except KeyError:
            station_cstate = utils.get_ground_station_reference_state_itrf(
                station=self.station, source="approximated_dsn"
            )

        # Extract observations and epochs
        observation_history: dict[ttime.Time, float] = {}
        for block in self.data.data_blocks:

            # Get interfaces to common and specific data
            common = block.common_data_block
            specific = block.observable_specific_data_block

            # Discard block if data is invalid
            if common.is_invalid:
                continue

            # Discard if not Doppler data type
            if not isinstance(specific, tdata.OdfDopplerDataBlock):
                continue

            # Discard block if undesired data type
            if common.data_type != data_type:
                continue

            # Discard block if incorrect station
            if common.receiving_station_id != self.station_id:
                continue

            # Discard block if out of time range
            epoch = reference_epoch + common.observable_time
            if time_range is not None:
                if (epoch < time_range[0]) or (epoch > time_range[1]):
                    continue

            # Transform epoch to TDB
            tdb_epoch = self.time_converter.convert_time_object(
                input_scale=ttime.TimeScales.utc_scale,
                output_scale=ttime.TimeScales.tdb_scale,
                input_value=epoch,
                earth_fixed_position=station_cstate[:3],
            )

            # Compute observation value from raw data
            turnaround_ratio = self.get_turnaround_ratio_for_odf_block(common)
            observation = turnaround_ratio * specific.reference_frequency + (
                self.observable_sign * common.observable_value
            )

            # Update collection
            observation_history[tdb_epoch] = observation

        return observation_history
