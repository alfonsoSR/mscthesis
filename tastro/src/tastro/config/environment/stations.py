from ..core import SetupBase
from tudatpy.astro import time_representation as ttime


class StationCoordinatesSetup(SetupBase):

    reference_epoch: ttime.Time
    available_position: str
    itrf_version: str
    linear_motion: bool
    body_deformation: bool


class UplinkFrequencySetup(SetupBase):

    present: bool
    start: list[ttime.Time]
    end: list[ttime.Time]
    ref_freq: list[float]
    ramp_rate: list[float]


class StationSetup(SetupBase):

    present: bool
    coordinates: StationCoordinatesSetup
    uplink: UplinkFrequencySetup
