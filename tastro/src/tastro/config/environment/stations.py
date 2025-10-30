from ..core import SetupBase
from tudatpy.astro import time_representation as ttime


class StationCoordinatesSetup(SetupBase):

    reference_epoch: ttime.Time
    available_position: str
    itrf_version: str
    linear_motion: bool
    body_deformation: bool


class StationLightTimeCorrectionsSetup(SetupBase):

    tropospheric_correction: bool
    ionospheric_correction: bool
    relativistic_correction: bool


class StationSetup(SetupBase):

    present: bool
    coordinates: StationCoordinatesSetup
    ltcorrections: StationLightTimeCorrectionsSetup
