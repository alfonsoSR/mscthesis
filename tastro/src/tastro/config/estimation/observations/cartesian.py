from ...core import SetupBase, SetupCollectionBase
from pathlib import Path
from dataclasses import dataclass
from tudatpy.astro import time_representation as ttime


class CartesianSourceSetup(SetupBase):

    path: Path = NotImplemented
    link: str = NotImplemented


@dataclass
class CartesianObservationsSetup(SetupBase):

    sources: list[CartesianSourceSetup] = NotImplemented
