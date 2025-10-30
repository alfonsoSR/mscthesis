from ..core import SetupBase
from .planets import PlanetSetup
from .stations import StationSetup
from .vehicles import VehicleSetup


class EnvironmentGeneralSetup(SetupBase):

    global_frame_origin: str
    global_frame_orientation: str
    spacecraft: str
    center: str


class EnvironmentSetup(SetupBase):

    general: EnvironmentGeneralSetup
    planets: dict[str, PlanetSetup]
    stations: dict[str, StationSetup]
    vehicles: dict[str, VehicleSetup]
