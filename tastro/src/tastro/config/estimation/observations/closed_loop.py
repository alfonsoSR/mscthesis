from ...core import SetupBase, SetupCollectionBase
from pathlib import Path
from dataclasses import dataclass
from tudatpy.astro import time_representation as ttime


class ODFSourceSetup(SetupBase):

    path: Path = NotImplemented
    station: str = NotImplemented


class SourcesSetup(SetupBase):

    ifms: list[Path] = NotImplemented
    odf: list[ODFSourceSetup] = NotImplemented


class UplinkFrequencySetup(SetupBase):

    start: list[ttime.Time] = NotImplemented
    end: list[ttime.Time] = NotImplemented
    ref_freq: list[float] = NotImplemented


class ObservationCompressionSetup(SetupBase):

    integration_time: ttime.Time = NotImplemented


@dataclass
class ClosedLoopObservationsSetup(SetupCollectionBase):

    sources: SourcesSetup
    uplink: dict[str, UplinkFrequencySetup]
    compression: ObservationCompressionSetup
