from dataclasses import dataclass
import numpy as np
from pathlib import Path
import yaml
from typing import Any


@dataclass
class PrefitSettings:

    general: dict[str, str]
    time: dict[str, float]
    ephemerides: dict[str, str]
    bodies: dict[str, dict[str, str]]
    stations: dict[str, dict[str, bool]]
    light_time: dict[str, list[str]]
    observations: dict[str, int] | None = None
    plotting: dict[str, Any] | None = None
    uplink: dict[str, dict[str, list[str]]] | None = None


def load_configuration(config_file: Path) -> PrefitSettings:

    if not config_file.exists():
        raise FileNotFoundError(f"Config file {config_file} does not exist")

    content = yaml.safe_load(config_file.open())

    observations = (
        content["ObservationCollection"]
        if "ObservationCollection" in content
        else None
    )
    plotting = content["Plotting"] if "Plotting" in content else None
    uplink = (
        content["FrequencyRamping"] if "FrequencyRamping" in content else None
    )

    return PrefitSettings(
        general=content["General"],
        time=content["Time"],
        ephemerides=content["Ephemerides"],
        bodies=content["BodySettings"],
        stations=content["GroundStationSettings"],
        light_time=content["LightTimeCorrections"],
        plotting=plotting,
        observations=observations,
        uplink=uplink,
    )


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


class PrefitResults:

    def __init__(self, source_file: Path) -> None:

        content = np.load(source_file)

        self.cepochs: np.ndarray = content[0]
        self.cobs: np.ndarray = content[1]
        self.robs: np.ndarray = content[2]
        self.repochs: np.ndarray = content[3]
        self.rres: np.ndarray = content[4]
        self.cres: np.ndarray = content[5]

        return None
