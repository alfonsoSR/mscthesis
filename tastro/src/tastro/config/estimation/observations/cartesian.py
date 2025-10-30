from ...core import SetupBase
from pathlib import Path

from tudatpy.astro import time_representation as ttime


class CartesianSourceSetup(SetupBase):

    path: Path
    link: str
    use_ephemerides: bool


class CartesianObservationsSetup(SetupBase):

    sources: list[CartesianSourceSetup]
