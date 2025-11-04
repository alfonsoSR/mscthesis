from ...core import SetupBase


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
