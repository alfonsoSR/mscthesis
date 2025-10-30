from pathlib import Path
from tudatpy.estimation.observations import SingleObservationSet
from tudatpy.astro import time_representation as ttime
from tudatpy import data as tdata
from dataclasses import dataclass
import numpy as np
from tudatpy.estimation.observations_setup import ancillary_settings as tancs
from ...config import CaseSetup
from tudatpy.dynamics.environment import PiecewiseLinearFrequencyInterpolator
import typing
from tudatpy.dynamics.environment_setup import ground_station as tgss


def get_ground_station_reference_state_itrf(
    station: str, source: str = "approximated_dsn"
) -> np.ndarray:

    # Get reference state for station
    match source:

        case "from_glo":

            xpos = tgss.get_vlbi_station_positions()[station]
            xvel = tgss.get_vlbi_station_velocities()[station]
            xvel /= 1000.0 * 365.0 * 86400.0
            return np.array([*xpos, *xvel])

        case "approximated_dsn":

            xpos = tgss.get_approximate_dsn_ground_station_positions()[station]
            xvel = np.zeros(3)
            return np.array([*xpos, *xvel])

        case _:
            raise NotImplementedError(f"Invalid type of available position: {source}")


def get_metadata_from_ifms_filename(source_file: Path) -> dict[str, str]:

    # Extract metadata from file name
    filename_components: dict[str, str] = {
        "mission_id": source_file.stem[0],
        "station_id": source_file.stem[1:3],
        "data_source_id": source_file.stem[3:7],
        "data_level": source_file.stem[7:10],
        "data_type": source_file.stem[11:14],
        "ref_epoch": source_file.stem[15:24],
        "version": source_file.stem[25:],
    }

    return filename_components


def sort_ifms_files_by_epoch(ifms_list: list[Path]) -> list[Path]:

    mask = np.argsort(
        [int(get_metadata_from_ifms_filename(file)["ref_epoch"]) for file in ifms_list]
    )
    return [ifms_list[idx] for idx in mask]


def identify_station_from_id(station_id: int) -> str:

    match station_id:
        case 32:
            station_name = "NWNORCIA"
        case 62:
            station_name = "CEBREROS"
        case 63:
            station_name = "DSS63"
        case 84:
            station_name = "MALARGUE"
        case 14:
            station_name = "DSS-14"
        case 65:
            station_name = "DSS65"
        case _:
            raise ValueError(f"Invalid station ID: {station_id}")

    return station_name


def _define_odf_reference_epoch_from_header(
    odf_content: tdata.OdfRawFileContents,
) -> ttime.Time:

    # Fail if not EME50
    if (
        odf_content.file_reference_date != 19500101
        or odf_content.file_reference_time != 0
    ):
        raise NotImplementedError("ODF reference epoch is not EME50")

    # Return reference epoch as time object
    return ttime.DateTime.from_iso_string("1950-01-01T00:00:00").to_epoch_time_object()


def _get_turnaround_ratio_for_odf_block(
    common: tdata.OdfCommonDataBlock,
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


def _get_odf_id_from_station_name(station_name: str) -> int:

    upper_name = station_name.upper()
    match upper_name:
        case "NWNORCIA":
            return 32
        case "CEBREROS":
            return 62
        case "MALARGUE":
            return 84
        case _:
            if "DSS-" in upper_name:
                return int(upper_name.split("-")[-1])
            elif "DSS" in upper_name:
                return int(upper_name[3:])
            else:
                raise ValueError(f"Invalid station name: {station_name}")


class RawDopplerObservationRecord(dict[ttime.Time, float]):

    def __init__(self, epochs: np.ndarray, values: np.ndarray) -> None:

        # Initialize dictionary
        super().__init__({ti: vi for ti, vi in zip(epochs, values)})

        # Initialize arrays of epochs and values
        self.epochs = epochs
        self.observations = values

        return None

    @staticmethod
    def from_ifms_file(ifms_file: Path) -> "RawDopplerObservationRecord":

        # Load raw data from IFMS file
        raw_data = tdata.read_ifms_file(
            str(ifms_file), apply_tropospheric_correction=False
        ).raw_datamap

        # Load observation epochs and values
        observation_epochs = np.array(
            [ttime.Time(float(epoch)) for epoch in raw_data["tdb_seconds_since_j2000"]]
        )
        observation_values = np.array(
            raw_data["doppler_averaged_frequency_hz"], dtype=float
        )

        return RawDopplerObservationRecord(observation_epochs, observation_values)

    @staticmethod
    def from_odf_file(
        odf_file: Path,
        rx_station: str,
        data_type: tdata.OdfDataType,
        time_range: tuple[ttime.Time, ttime.Time] | None = None,
    ) -> "RawDopplerObservationRecord":

        # Load raw contents
        raw_data = tdata.read_odf_file(str(odf_file))

        # Retrieve reference epoch from header
        reference_epoch = _define_odf_reference_epoch_from_header(raw_data)

        # Retrieve ODF code for station
        station_id = _get_odf_id_from_station_name(rx_station)

        # Get Earth-fixed position of ground station
        try:
            station_cstate = get_ground_station_reference_state_itrf(
                station=rx_station, source="from_glo"
            )
        except KeyError as err:
            station_cstate = get_ground_station_reference_state_itrf(rx_station)

        # Initialize time-scale converter
        time_converter = ttime.default_time_scale_converter()

        # Initialize containers of observation epochs and values
        observation_epochs: list[ttime.Time] = []
        observation_values: list[float] = []

        # Extract observations and observation epochs
        for block in raw_data.data_blocks:

            # Get interfaces to common and specific data
            common = block.common_data_block
            specific = block.observable_specific_data_block

            # Discard block if data is invalid
            if common.is_invalid:
                continue

            # Discard if not Doppler data type
            if not isinstance(specific, tdata.OdfDopplerDataBlock):
                continue

            # Discard block if incorrect data type
            if common.data_type != data_type:
                continue

            # Discard block if incorrect station
            if common.receiving_station_id != station_id:
                continue

            # Discard block if out of time range
            epoch = reference_epoch + common.observable_time
            if time_range is not None:
                if (epoch < time_range[0]) or (epoch > time_range[1]):
                    continue

            # Transform epoch to TDB and add to container
            observation_epochs.append(
                time_converter.convert_time_object(
                    input_scale=ttime.TimeScales.utc_scale,
                    output_scale=ttime.TimeScales.tdb_scale,
                    input_value=epoch,
                    earth_fixed_position=station_cstate[:3],
                )
            )

            # Add observation value to container
            turnaround_ratio = _get_turnaround_ratio_for_odf_block(common)
            observation_values.append(
                turnaround_ratio * specific.reference_frequency
                - common.observable_value
            )

        # Cast containers to numpy arrays and return
        return RawDopplerObservationRecord(
            epochs=np.array(observation_epochs),
            values=np.array(observation_values),
        )


class DopplerObservationRecord:

    def __init__(
        self,
        epochs: np.ndarray,
        observations: np.ndarray,
        ramping_start_epochs: list[ttime.Time],
        ramping_end_epochs: list[ttime.Time],
        ramping_freq: list[float],
    ) -> None:

        # Define observation epochs and values (Sorted by epoch)
        idx_order = np.argsort([ti.to_float() for ti in epochs])
        self.epochs = epochs[idx_order]
        self.observations = observations[idx_order]

        # Define uplink frequency interpolator
        self.uplink_interpolator = PiecewiseLinearFrequencyInterpolator(
            start_times=ramping_start_epochs,
            end_times=ramping_end_epochs,
            ramp_rates=np.zeros(len(ramping_freq)).tolist(),
            start_frequency=ramping_freq,
        )

        return None

    @classmethod
    def from_config(
        cls,
        epochs: np.ndarray,
        observations: np.ndarray,
        station: str,
        config: CaseSetup,
    ) -> "DopplerObservationRecord":

        uplink_config = config.estimation.observations.closed_loop.uplink[station]

        return DopplerObservationRecord(
            epochs=epochs,
            observations=observations,
            ramping_start_epochs=uplink_config.start,
            ramping_end_epochs=uplink_config.end,
            ramping_freq=uplink_config.ref_freq,
        )

    @classmethod
    def from_raw_record_and_config(
        cls,
        raw_record: RawDopplerObservationRecord,
        station: str,
        config: CaseSetup,
    ) -> "DopplerObservationRecord":

        uplink_config = config.estimation.observations.closed_loop.uplink[station]

        return DopplerObservationRecord(
            epochs=raw_record.epochs,
            observations=raw_record.observations,
            ramping_start_epochs=uplink_config.start,
            ramping_end_epochs=uplink_config.end,
            ramping_freq=uplink_config.ref_freq,
        )
