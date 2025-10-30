from ..core import SetupBase

from tudatpy.astro import time_representation as ttime


class PlanetEphemeridesSetup(SetupBase):

    present: bool
    model: str
    ephemeris_frame_origin: str
    ephemeris_frame_orientation: str
    interpolation_step: ttime.Time
    interpolation_buffer: ttime.Time


class PlanetRotationSetup(SetupBase):

    present: bool
    model: str
    base_frame: str
    target_frame: str


class PlanetShapeSetup(SetupBase):

    present: bool
    model: str
    equatorial_radius: float
    flattening_factor: float


class PlanetGravitySetup(SetupBase):

    present: bool
    model: str
    spherical_harmonics_file: str
    spherical_harmonics_degree: int
    spherical_harmonics_order: int
    spherical_harmonics_frame: str


class DirectRadiationSourceSetup(SetupBase):

    luminosity: float


class PlanetRadiationSourceSetup(SetupBase):

    present: bool
    model: str
    direct_setup: DirectRadiationSourceSetup


class PlanetExponentialAtmosphereSetup(SetupBase):

    scale_height: float
    surface_density: float
    constant_temperature: float
    specific_gas_constant: float
    ratio_specific_heats: float


class PlanetAtmosphereSetup(SetupBase):

    present: bool
    model: str
    exponential_settings: PlanetExponentialAtmosphereSetup


class PlanetSetup(SetupBase):

    present: bool
    ephemerides: PlanetEphemeridesSetup
    rotation: PlanetRotationSetup
    shape: PlanetShapeSetup
    gravity: PlanetGravitySetup
    radiation: PlanetRadiationSourceSetup
    atmosphere: PlanetAtmosphereSetup
