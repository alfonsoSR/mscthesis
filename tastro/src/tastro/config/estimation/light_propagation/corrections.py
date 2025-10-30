from ...core import SetupBase

from pathlib import Path


class TroposphericCorrectionSetup(SetupBase):

    present: bool
    model: str
    sources: list[Path]


class IonosphericCorrectionSetup(SetupBase):

    present: bool
    model: str
    sources: list[Path]


class RelativisticCorrectionSetup(SetupBase):

    present: bool
    model: str
    bodies: list[str]


class LightTimeCorrectionsSetup(SetupBase):

    tropospheric: TroposphericCorrectionSetup
    ionospheric: IonosphericCorrectionSetup
    relativistic: RelativisticCorrectionSetup
