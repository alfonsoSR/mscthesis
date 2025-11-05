from ...core import SetupBase
from pathlib import Path


class ObservationSourceSetup(SetupBase):

    path: Path
    station: str
