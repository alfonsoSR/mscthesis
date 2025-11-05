import numpy as np
from tudatpy.dynamics.environment_setup import ground_station as tgss
from .logging import log


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


def get_id_from_station_name(station_name: str) -> int:

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


def get_ground_station_reference_state_itrf(
    station: str, source: str = "approximated_dsn"
) -> np.ndarray:

    # Get reference state for station
    match source:

        case "from_glo":

            log.debug("Station coordinates from glo")

            xpos = tgss.get_vlbi_station_positions()[station]
            xvel = tgss.get_vlbi_station_velocities()[station]
            xvel /= 1000.0 * 365.0 * 86400.0
            return np.array([*xpos, *xvel])

        case "approximated_dsn":

            log.debug("Approximated station coordinates")

            xpos = tgss.get_approximate_dsn_ground_station_positions()[station]
            xvel = np.zeros(3)
            return np.array([*xpos, *xvel])

        case _:
            raise NotImplementedError(
                f"Invalid type of available position: {source}"
            )
