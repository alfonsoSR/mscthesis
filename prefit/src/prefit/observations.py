from tudatpy import data as tdata
from pathlib import Path
from tudatpy.astro import time_representation as ttime
from tudatpy.estimation import observations_setup as tobss
from .utils import transform_utc_epochs_to_tdb
from .io import (
    TwoWayDopplerObservations,
    sort_ifms_files_by_epoch,
    get_metadata_from_ifms_filename,
)
import numpy as np


def __define_odf_reference_epoch_from_header(
    odf_content: tdata.OdfRawFileContents,
) -> ttime.Time:

    # Fail if not EME50
    if (
        odf_content.file_reference_date != 19500101
        or odf_content.file_reference_time != 0
    ):
        raise NotImplementedError("ODF reference epoch is not EME50")

    # Return reference epoch as time object
    return ttime.DateTime.from_iso_string(
        "1950-01-01T00:00:00"
    ).to_epoch_time_object()


def __get_turnaround_ratio_for_odf_block(
    common: tdata.OdfCommonDataBlock,
) -> float:

    # Identify uplink band
    if common.uplink_band_id == 0:
        uplink_band = tobss.ancillary_settings.FrequencyBands.ku_band
    else:
        uplink_band = tobss.ancillary_settings.FrequencyBands(
            common.uplink_band_id - 1
        )

    # Identify downlink band
    if common.downlink_band_id == 0:
        downlink_band = tobss.ancillary_settings.FrequencyBands.ku_band
    else:
        downlink_band = tobss.ancillary_settings.FrequencyBands(
            common.downlink_band_id - 1
        )

    # Return default turnaround ratio based on bands
    return tobss.ancillary_settings.dsn_default_turnaround_ratios(
        uplink_band=uplink_band, downlink_band=downlink_band
    )


def load_two_way_doppler_data_from_odf_file(
    odf_file: Path,
    time_range: tuple[ttime.Time, ttime.Time] | None = None,
    data_container: dict[int, dict[str, list]] | None = None,
) -> dict[int, dict[str, list]]:

    # Load raw data
    content = tdata.read_odf_file(str(odf_file))

    # Get reference epoch from header
    reference_epoch = __define_odf_reference_epoch_from_header(content)

    # Initialize time-converter
    time_converter = ttime.default_time_scale_converter()

    # Initialize dictionary for ODF data
    if data_container is not None:
        data_per_station = data_container
    else:
        data_per_station: dict[int, dict[str, list]] = {}

    # Extract observations and observation epochs
    for block in content.data_blocks:

        # Get interfaces to common and specific data
        common = block.common_data_block
        specific = block.observable_specific_data_block

        # Discard block if data is invalid
        if common.is_invalid:
            continue

        # Discard block if it does not contain Doppler data
        if not isinstance(specific, tdata.OdfDopplerDataBlock):
            continue

        # Discard block if Doppler is not two-way
        if common.data_type != tdata.OdfDataType.two_way_doppler:
            continue

        # Discard block if it is out of time range
        observation_epoch = reference_epoch + common.observable_time
        if time_range is not None:
            if (observation_epoch < time_range[0]) or (
                observation_epoch > time_range[1]
            ):
                continue

        # If station is not in dictionary, initialize entry
        station_id = common.transmitting_station_id
        if station_id not in data_per_station:
            data_per_station[station_id] = {
                "obs_epochs_tdb": [],
                "observations": [],
                "ramp_start_tdb": [],
                "ramp_end_tdb": [],
                "ramp_f0": [],
                "ramp_df": [],
            }

        # Get pointer to dictionary entry for current station
        station_data = data_per_station[station_id]

        # Update observation epochs
        station_data["obs_epochs_tdb"].append(
            time_converter.convert_time_object(
                ttime.TimeScales.utc_scale,
                ttime.TimeScales.tdb_scale,
                observation_epoch,
            )
        )

        # Update observations
        turnaround_ratio = __get_turnaround_ratio_for_odf_block(common)
        station_data["observations"].append(
            turnaround_ratio * specific.reference_frequency
            - common.observable_value
        )

    # Load ramping data
    for station_id, station_blocks in content.ramp_blocks.items():

        # Skip if station has no observations
        if station_id not in data_per_station:
            continue

        # Get pointer to entry for current station
        station_data = data_per_station[station_id]

        # Add contents of blocks to dictionary
        index_of_last_block = len(station_blocks) - 1
        for idx, block in enumerate(station_blocks):

            # Get start and end of ramping interval
            ramp_start = time_converter.convert_time_object(
                ttime.TimeScales.utc_scale,
                ttime.TimeScales.tdb_scale,
                reference_epoch + block.ramp_start_epoch,
            )
            ramp_end = time_converter.convert_time_object(
                ttime.TimeScales.utc_scale,
                ttime.TimeScales.tdb_scale,
                reference_epoch + block.ramp_end_epoch,
            )

            # Update dictionary with initial and final ramping epochs in TDB
            station_data["ramp_start_tdb"].append(ramp_start)
            station_data["ramp_end_tdb"].append(ramp_end)

            # Update with constant and linear ramping terms
            station_data["ramp_f0"].append(block.ramp_start_frequency)
            station_data["ramp_df"].append(block.ramp_rate)

    return data_per_station


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
            station_name = "DSS14"
        case 65:
            station_name = "DSS65"
        case _:
            raise ValueError(f"Invalid station ID: {station_id}")

    return station_name


def load_odf_data_per_station(
    odf_files: list[Path],
    data_per_station: dict[str, TwoWayDopplerObservations] | None = None,
) -> dict[str, TwoWayDopplerObservations]:

    # Load data into dictionaries
    raw_data_per_station: dict[int, dict[str, list]] = {}
    for odf_file in odf_files:
        raw_data_per_station = load_two_way_doppler_data_from_odf_file(
            odf_file=odf_file, data_container=raw_data_per_station
        )

    # Initialize output dictionary
    if data_per_station is None:
        output: dict[str, TwoWayDopplerObservations] = {}
    else:
        output = data_per_station

    # Fill dictionary with observation data
    for station_id, station_data in raw_data_per_station.items():

        # Get name of station from ID
        station_name = identify_station_from_id(station_id)

        # Fail if station is already in dictionary
        if station_name in output:
            raise NotImplementedError(
                f"Station {station_name} already present in dictionary"
            )

        # Generate TwoWayDopplerObservations object from data
        station_observations = TwoWayDopplerObservations(
            np.array(station_data["obs_epochs_tdb"]),
            np.array(station_data["observations"]),
            np.array(station_data["ramp_f0"]),
            np.array(station_data["ramp_df"]),
            np.array(station_data["ramp_start_tdb"]),
            np.array(station_data["ramp_end_tdb"]),
            np.zeros(len(station_data["obs_epochs_tdb"])),
        )

        # Sort observations
        _float_obs_epochs = np.array(
            [ti.to_float() for ti in station_observations.observation_epochs_et]
        )
        _obs_mask = np.argsort(_float_obs_epochs)
        station_observations.observation_epochs_et = (
            station_observations.observation_epochs_et[_obs_mask]
        )
        station_observations.observation_values = (
            station_observations.observation_values[_obs_mask]
        )

        # Sort ramping data
        _float_ramp_epochs = np.array(
            [ti.to_float() for ti in station_observations.ramping_start_tdb]
        )
        _ramp_mask = np.argsort(_float_ramp_epochs)
        station_observations.ramping_start_tdb = (
            station_observations.ramping_start_tdb[_ramp_mask]
        )
        station_observations.ramping_stop_tdb = (
            station_observations.ramping_stop_tdb[_ramp_mask]
        )
        station_observations.ramping_f0 = station_observations.ramping_f0[
            _ramp_mask
        ]
        station_observations.ramping_df = station_observations.ramping_df[
            _ramp_mask
        ]

        # Add entry to dictionary
        output[station_name] = station_observations

    return output


def load_doppler_observations_from_list_of_ifms_files(
    source_files: list[Path], ground_station_config: dict[str, dict[str, bool]]
) -> TwoWayDopplerObservations:

    # Sort files by epoch
    source_files = sort_ifms_files_by_epoch(source_files)

    # Initialize containers
    observation_epochs_et = []
    observation_values = []
    _ramping_f0 = []
    _ramping_df = []
    _ramping_utc0 = []
    residuals = []
    tropo_correction = []

    # Fill containers with data from files
    for ifms in source_files:

        # Load content into object
        content = tdata.read_ifms_file(
            str(ifms), apply_tropospheric_correction=False
        ).raw_datamap
        print(ifms.name)

        # Update containers with observation data
        observation_epochs_et += content["tdb_seconds_since_j2000"]
        observation_values += content["doppler_averaged_frequency_hz"]
        residuals += content["doppler_noise_hz"]
        tropo_correction += content["doppler_troposphere_correction"]

        # Update containers with ramping data
        _ramping_f0.append(content["transmission_frequency_constant_term"])
        _ramping_df.append(content["transmission_frequency_linear_term"])
        _ramping_utc0.append(content["ramp_reference_time"])

    # Fill gaps between observation intervals
    ramping_f0: list[str] = []
    ramping_df: list[str] = []
    ramping_utc0: list[str] = []

    if len(_ramping_f0) == 1:
        ramping_f0 = _ramping_f0[0]
        ramping_df = _ramping_df[0]
        ramping_utc0 = _ramping_utc0[0]
    else:
        for idx, f0 in enumerate(_ramping_f0[:-1]):

            # Get current and next values for reference frequency
            current_f0 = f0
            next_f0 = _ramping_f0[idx + 1]

            # Get current and next values for ramp
            current_df = _ramping_df[idx]
            next_df = _ramping_df[idx + 1]

            # Get current and next sets of reference epochs
            current_ref = _ramping_utc0[idx].copy()
            next_ref = _ramping_utc0[idx + 1].copy()

            # Get intermediate epoch between intervals
            last_epoch_current = ttime.DateTime.from_iso_string(current_ref[-1])
            first_epoch_next = ttime.DateTime.from_iso_string(next_ref[0])
            diff = (
                first_epoch_next.to_epoch_time_object()
                - last_epoch_current.to_epoch_time_object()
            )
            intermediate_epoch = last_epoch_current.to_epoch_time_object() + (
                0.5 * diff
            )
            # print(
            #     current_ref[-1],
            #     next_ref[0],
            #     diff.to_float(),
            #     ttime.DateTime.from_epoch_time_object(
            #         intermediate_epoch
            #     ).to_iso_string(),
            # )
            if diff <= ttime.Time(2):
                raise NotImplementedError(
                    f"We want separation by more than 1 second: {diff.to_float()} : {first_epoch_next.to_iso_string()} : {last_epoch_current.to_iso_string()}"
                )

            # Extend current at intermediate -1
            current_ref.append(
                ttime.DateTime.from_epoch_time_object(
                    intermediate_epoch - ttime.Time(1)
                ).to_iso_string()
            )
            next_ref_first = ttime.DateTime.from_epoch_time_object(
                intermediate_epoch
            ).to_iso_string()
            # print(
            #     current_ref[-1],
            #     _ramping_utc0[idx + 1][-1],
            #     ttime.DateTime.from_epoch_time_object(
            #         intermediate_epoch
            #     ).to_iso_string(),
            # )

            # Add to containers
            ramping_f0 += current_f0
            ramping_f0.append(current_f0[-1])
            ramping_f0.append(next_f0[0])

            ramping_df += current_df
            ramping_df.append(current_df[-1])
            ramping_df.append(next_df[0])

            ramping_utc0 += current_ref
            ramping_utc0.append(next_ref_first)

    # Remove invalid values
    for idx, val in enumerate(ramping_df):
        if val == "-99999.999999":
            ramping_df[idx] = "0"

    # Transform reference ramping epochs to TDB
    # ramping_utc0 = [
    #     ttime.DateTime.from_iso_string(utci).to_epoch_time_object()
    #     for utci in ramping_utc0
    # ]
    ramping_tdb0 = transform_utc_epochs_to_tdb(
        [
            ttime.DateTime.from_iso_string(utci).to_epoch_time_object()
            for utci in ramping_utc0
        ],
        np.array([0.0, 0.0, 0.0]),
    )

    # Get station name from file metadata
    _station_id = get_metadata_from_ifms_filename(source_files[0])["station_id"]
    station_name = identify_station_from_id(int(_station_id))

    # Calculate observation values based on configuration
    obs_values = np.array(observation_values, dtype=float)
    if bool(ground_station_config[station_name]["troposphere"]):
        print(f"Correcting the tropo for {station_name}")
        obs_values -= np.array(tropo_correction, dtype=float)

    return TwoWayDopplerObservations(
        observation_epochs_et=np.array(
            [ttime.Time(float(ti)) for ti in observation_epochs_et]
        ),
        observation_values=obs_values,
        ramping_f0=np.array(ramping_f0, dtype=float)[:-1],
        ramping_df=np.array(ramping_df, dtype=float)[:-1],
        ramping_start_tdb=np.array(ramping_tdb0)[:-1],
        ramping_stop_tdb=np.array(ramping_tdb0)[1:],
        residuals=np.array(residuals, dtype=float),
    )
