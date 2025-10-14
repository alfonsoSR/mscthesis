from ...core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from pathlib import Path


class TroposphericCorrectionSetup(SetupBase):

    model: str = NotImplemented
    sources: list[Path] = NotImplemented


class IonosphericCorrectionSetup(SetupBase):

    model: str = NotImplemented
    sources: list[Path] = NotImplemented


class RelativisticCorrectionSetup(SetupBase):

    model: str = NotImplemented
    bodies: list[str] = NotImplemented


@dataclass
class LightTimeCorrectionsSetup(SetupCollectionBase):

    tropospheric: TroposphericCorrectionSetup
    ionospheric: IonosphericCorrectionSetup
    relativistic: RelativisticCorrectionSetup
