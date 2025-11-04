from ...core import SetupBase


class PlanetRotationSetup(SetupBase):

    present: bool
    model: str
    base_frame: str
    target_frame: str
