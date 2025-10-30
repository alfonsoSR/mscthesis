from ...core import SetupBase
from pathlib import Path

from tudatpy.astro import time_representation as ttime
from .filters import FiltersSetup


class ODFSourceSetup(SetupBase):

    path: Path
    station: str


class SourcesSetup(SetupBase):

    ifms: list[Path]
    odf: list[ODFSourceSetup]


class UplinkFrequencySetup(SetupBase):

    start: list[ttime.Time]
    end: list[ttime.Time]
    ref_freq: list[float]


class ObservationCompressionSetup(SetupBase):

    present: bool
    integration_time: ttime.Time


class ClosedLoopObservationsSetup(SetupBase):

    sources: SourcesSetup
    uplink: dict[str, UplinkFrequencySetup]
    compression: ObservationCompressionSetup
    filters: FiltersSetup
