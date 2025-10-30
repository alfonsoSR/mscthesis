from ...core import SetupBase

from .convergence import LightTimeConvergence
from .corrections import LightTimeCorrectionsSetup


class LightTimeSetup(SetupBase):

    corrections: LightTimeCorrectionsSetup
    convergence: LightTimeConvergence
