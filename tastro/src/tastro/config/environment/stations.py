from ..core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from tudatpy.astro import time_representation as ttime


class StationCoordinatesSetup(SetupBase):

    reference_epoch: ttime.Time = NotImplemented
    available_position: str = NotImplemented
    itrf_version: str = NotImplemented
    linear_motion: bool = NotImplemented
    body_deformation: bool = NotImplemented


class StationLightTimeCorrectionsSetup(SetupBase):

    tropospheric_correction: bool = NotImplemented
    ionospheric_correction: bool = NotImplemented
    relativistic_correction: bool = NotImplemented


@dataclass
class StationSetup(SetupCollectionBase):

    present: bool
    coordinates: StationCoordinatesSetup
    ltcorrections: StationLightTimeCorrectionsSetup
