from ....core import SetupBase
from pathlib import Path
import sys
import os


class CartesianSourceSetup(SetupBase):

    path: Path
    link: str
    use_ephemerides: bool


class CartesianObservationsSetup(SetupBase):

    sources: list[CartesianSourceSetup]
