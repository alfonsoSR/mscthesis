from ...core import SetupCollectionBase
from dataclasses import dataclass
from .convergence import LightTimeConvergence
from .corrections import LightTimeCorrectionsSetup


@dataclass
class LightTimeSetup(SetupCollectionBase):

    corrections: LightTimeCorrectionsSetup
    convergence: LightTimeConvergence
