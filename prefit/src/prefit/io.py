from dataclasses import dataclass
import numpy as np
from tudatpy.astro import time_representation as ttime
from pathlib import Path
from tudatpy import data as tdata
from .utils import transform_utc_epochs_to_tdb


@dataclass
class TwoWayDopplerObservations:

    observation_epochs_et: np.ndarray
    observation_values: np.ndarray
    ramping_f0: np.ndarray
    ramping_df: np.ndarray
    ramping_start_tdb: np.ndarray
    ramping_stop_tdb: np.ndarray
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
