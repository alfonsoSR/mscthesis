from dataclasses import dataclass
import numpy as np
from tudatpy.astro import time_representation as ttime
from pathlib import Path
from tudatpy import data as tdata
from .utils import transform_utc_epochs_to_et


@dataclass
class TwoWayDopplerObservations:

    observation_epochs_et: np.ndarray
    observation_values: np.ndarray
    ramping_f0: np.ndarray
    ramping_df: np.ndarray
    ramping_tdb0: list[ttime.Time]
    residuals: np.ndarray


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
        [
            int(get_metadata_from_ifms_filename(file)["ref_epoch"])
            for file in ifms_list
        ]
    )
    return [ifms_list[idx] for idx in mask]


def load_doppler_observations_from_ifms_files(
    source_files: list[Path],
) -> TwoWayDopplerObservations:

    # Sort files by epoch
    source_files = sort_ifms_files_by_epoch(source_files)

    # Initialize containers
    observation_epochs_et = []
    observation_values = []
    ramping_f0 = []
    ramping_df = []
    ramping_utc0 = []
    residuals = []
    tropo_correction = []

    # Fill containers with data from files
    for ifms in source_files:

        # Load content into object
        content = tdata.read_ifms_file(
            str(ifms), apply_tropospheric_correction=False
        ).raw_datamap

        # Update containers with data
        observation_epochs_et += content["tdb_seconds_since_j2000"]
        observation_values += content["doppler_averaged_frequency_hz"]
        ramping_f0 += content["transmission_frequency_constant_term"]
        ramping_df += content["transmission_frequency_linear_term"]
        ramping_utc0 += content["ramp_reference_time"]
        residuals += content["doppler_noise_hz"]
        tropo_correction += content["doppler_troposphere_correction"]

    # Remove invalid values
    for idx, val in enumerate(ramping_df):
        if val == "-99999.999999":
            ramping_df[idx] = "0"

    # Transform reference ramping epochs to TDB
    ramping_utc0 = [
        ttime.DateTime.from_iso_string(utci).to_epoch_time_object()
        for utci in ramping_utc0
    ]
    ramping_tdb0 = transform_utc_epochs_to_et(
        ramping_utc0, np.array([0.0, 0.0, 0.0])
    )

    return TwoWayDopplerObservations(
        observation_epochs_et=np.array(
            [ttime.Time(float(ti)) for ti in observation_epochs_et]
        ),
        observation_values=np.array(observation_values, dtype=float)
        - np.array(tropo_correction, dtype=float),
        ramping_f0=np.array(ramping_f0, dtype=float),
        ramping_df=np.array(ramping_df, dtype=float),
        ramping_tdb0=ramping_tdb0,
        residuals=np.array(residuals, dtype=float),
    )
