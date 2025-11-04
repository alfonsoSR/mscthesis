from ...core import SetupBase


class PlanetShapeSetup(SetupBase):

    present: bool
    model: str
    equatorial_radius: float
    flattening_factor: float
