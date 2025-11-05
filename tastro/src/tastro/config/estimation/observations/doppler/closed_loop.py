from ....core import SetupBase
from pathlib import Path
from ..common import ObservationSourceSetup
from ..filters import FiltersSetup
from tudatpy.astro import time_representation as ttime


class ClosedLoopSourcesSetup(SetupBase):

    ifms: list[Path]
    odf: list[ObservationSourceSetup]


class UplinkFrequencySetup(SetupBase):

    start: list[ttime.Time]
    end: list[ttime.Time]
    ref_freq: list[float]
    present: bool = False


class ObservationCompressionSetup(SetupBase):

    present: bool
    integration_time: ttime.Time


class ClosedLoopObservationsSetup(SetupBase):

    sources: ClosedLoopSourcesSetup
    uplink: dict[str, UplinkFrequencySetup]
    compression: ObservationCompressionSetup
    filters: FiltersSetup
