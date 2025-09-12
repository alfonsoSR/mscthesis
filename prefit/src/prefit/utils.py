import random
from tudatpy.dynamics.environment_setup import BodyListSettings
import spiceypy
from pathlib import Path
from tudatpy.astro import time_representation as ttime
from dataclasses import dataclass
import numpy as np
from tudatpy.dynamics import environment_setup as tenvs


def generate_random_color():
    """Generate a random color in hexadecimal format."""
    return "#{:02x}{:02x}{:02x}".format(
        random.randint(0, 255),  # Red
        random.randint(0, 255),  # Green
        random.randint(0, 255),  # Blue
    )


def display_body_settings(
    environment_settings: BodyListSettings,
    body_name: str,
    settings_type: str,
) -> None:

    complete_settings = environment_settings.get(body_name)
    settings = getattr(complete_settings, settings_type)

    if not isinstance(settings, list):
        settings = [settings]

    if len(settings) == 0:
        print(f"Body {body_name} does not have {settings_type}")
        return None

    for csettings in settings:
        print(f"Displaying {settings_type} for {body_name}: {csettings}")
        print("---------------------------------------")
        for item in dir(csettings):
            if item.startswith("_") or item.endswith("_"):
                continue
            print(f"{item}: {getattr(csettings, item)}")
        print("---------------------------------------")

    return None


def get_loaded_kernels(ktype: str = "all") -> list[Path]:

    total = spiceypy.ktotal(ktype)

    if not total:
        print(f"No kernels of type {ktype} have been loaded")
        return []

    output: list[Path] = []
    for idx in range(total):
        output.append(Path(spiceypy.kdata(idx, ktype)[0]))

    return output


def get_epochs_in_et(
    initial_utc_epoch: ttime.DateTime, step: ttime.Time, duration_in_steps: int
) -> list[ttime.Time]:

    # Time conversion settings
    time_converter = ttime.default_time_scale_converter()
    j2000_tt = ttime.DateTime.from_iso_string(
        "2000-01-01T12:00:00"
    ).to_epoch_time_object()
    j2000_tdb = time_converter.convert_time_object(
        ttime.TimeScales.tt_scale, ttime.TimeScales.tdb_scale, j2000_tt
    )

    # Time interval (UTC)
    initial_epoch = initial_utc_epoch.to_epoch_time_object()
    epochs_utc = [
        initial_epoch + step * i for i in range(duration_in_steps + 1)
    ]

    # Time interval ET (TDB)
    epochs_tdb = [
        time_converter.convert_time_object(
            ttime.TimeScales.utc_scale, ttime.TimeScales.tdb_scale, epoch
        )
        for epoch in epochs_utc
    ]
    epochs_et = [epoch - j2000_tdb for epoch in epochs_tdb]

    return epochs_et


def transform_utc_epochs_to_tdb(
    utc_epochs: list[ttime.Time], station_coordinates: np.ndarray
) -> list[ttime.Time]:

    # Initialize time converter
    converter = ttime.default_time_scale_converter()

    # Perform transformation
    return [
        converter.convert_time_object(
            input_scale=ttime.TimeScales.utc_scale,
            output_scale=ttime.TimeScales.tdb_scale,
            input_value=utc_epoch,
            earth_fixed_position=station_coordinates,
        )
        for utc_epoch in utc_epochs
    ]


def get_station_reference_state_itrf2000(station_name: str) -> np.ndarray:
    """Get ITRF2000 state of station at 2000-01-01T00:00:00

    The state of the ground stations is read from the catalogs in
    `.tudat/resource/station_locations`

    :param station_name: Name of the ground station as defined in the catalog
    :return reference_state: Cartesian state of station at reference epoch in ITRF2000
    """

    evn_positions = tenvs.ground_station.get_vlbi_station_positions()
    evn_velocities = tenvs.ground_station.get_vlbi_station_velocities()

    if station_name not in evn_positions:
        raise KeyError(
            f"The ITRF2000 state of station {station_name} is not available"
        )

    # Transform position and velocity to SI units
    station_position = evn_positions[station_name]
    station_velocity = evn_velocities[station_name] / (1000 * 365 * 86400)

    return np.array([station_position, station_velocity]).flatten()
