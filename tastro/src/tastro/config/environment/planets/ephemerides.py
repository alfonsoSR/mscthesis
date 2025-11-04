from ...core import SetupBase
from tudatpy.astro import time_representation as ttime


class PlanetEphemeridesSetup(SetupBase):

    present: bool
    model: str
    ephemeris_frame_origin: str
    ephemeris_frame_orientation: str
    interpolation_step: ttime.Time
    interpolation_buffer: ttime.Time
