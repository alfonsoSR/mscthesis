from ..core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from .planets import PlanetSetup
from .stations import StationSetup
from .vehicles import VehicleSetup


class EnvironmentGeneralSetup(SetupBase):

    global_frame_origin: str = NotImplemented
    global_frame_orientation: str = NotImplemented
    spacecraft: str = NotImplemented
    center: str = NotImplemented


@dataclass
class EnvironmentSetup(SetupCollectionBase):

    general: EnvironmentGeneralSetup
    planets: dict[str, PlanetSetup]
    stations: dict[str, StationSetup]
    vehicles: dict[str, VehicleSetup]
