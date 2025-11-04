from ...core import SetupBase


class DirectRadiationSourceSetup(SetupBase):

    luminosity: float


class PlanetRadiationSourceSetup(SetupBase):

    present: bool
    model: str
    direct_setup: DirectRadiationSourceSetup
