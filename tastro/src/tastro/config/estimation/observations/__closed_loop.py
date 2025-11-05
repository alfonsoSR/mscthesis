raise DeprecationWarning("This module is deprecated")

from ...core import SetupBase
from pathlib import Path

from tudatpy.astro import time_representation as ttime
from .filters import FiltersSetup
from .common import ObservationSourceSetup


class SourcesSetup(SetupBase):

    ifms: list[Path]
    odf: list[ObservationSourceSetup]


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
