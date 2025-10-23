from ..core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from tudatpy.astro import time_representation as ttime


class PlanetEphemeridesSetup(SetupBase):

    model: str = NotImplemented
    ephemeris_frame_origin: str = NotImplemented
    ephemeris_frame_orientation: str = NotImplemented
    interpolation_step: ttime.Time = NotImplemented
    interpolation_buffer: ttime.Time = NotImplemented


class PlanetRotationSetup(SetupBase):

    model: str = NotImplemented
    base_frame: str = NotImplemented
    target_frame: str | None = NotImplemented


class PlanetShapeSetup(SetupBase):

    model: str = NotImplemented
    equatorial_radius: float = NotImplemented
    flattening_factor: float = NotImplemented


class PlanetGravitySetup(SetupBase):

    model: str = NotImplemented
    spherical_harmonics_file: str | None = NotImplemented
    spherical_harmonics_degree: int | None = NotImplemented
    spherical_harmonics_order: int | None = NotImplemented
    spherical_harmonics_frame: str | None = NotImplemented


class DirectRadiationSourceSetup(SetupBase):

    luminosity: float = NotImplemented


@dataclass
class PlanetRadiationSourceSetup(SetupCollectionBase):

    present: bool
    model: str
    direct_setup: DirectRadiationSourceSetup


class PlanetExponentialAtmosphereSetup(SetupBase):

    scale_height: float | None = NotImplemented
    surface_density: float | None = NotImplemented
    constant_temperature: float | None = NotImplemented
    specific_gas_constant: float | None = NotImplemented
    ratio_specific_heats: float | None = NotImplemented


@dataclass
class PlanetAtmosphereSetup(SetupCollectionBase):

    present: bool
    model: str
    exponential_settings: PlanetExponentialAtmosphereSetup


@dataclass
class PlanetSetup(SetupCollectionBase):

    present: bool
    ephemerides: PlanetEphemeridesSetup
    rotation: PlanetRotationSetup
    shape: PlanetShapeSetup
    gravity: PlanetGravitySetup
    radiation: PlanetRadiationSourceSetup
    atmosphere: PlanetAtmosphereSetup
