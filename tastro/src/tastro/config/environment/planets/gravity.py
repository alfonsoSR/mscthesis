from ...core import SetupBase


class SolidK2TideSetup(SetupBase):

    tide_raising_bodies: list[str]
    k2: float
    present: bool = False


class SolidTideSetup(SetupBase):

    present: bool
    model: str
    solid_k2_setup: SolidK2TideSetup


class PlanetGravitySetup(SetupBase):

    present: bool
    model: str
    spherical_harmonics_file: str
    spherical_harmonics_degree: int
    spherical_harmonics_order: int
    spherical_harmonics_frame: str
