from ....core import SetupBase
from tudatpy.estimation.observations_setup import ancillary_settings as tanc
from tudatpy.astro import time_representation as ttime


class DopplerAncillarySettingsSetup(SetupBase):

    uplink_band: tanc.FrequencyBands = NotImplemented
    downlink_band: tanc.FrequencyBands = NotImplemented
    reference_band: tanc.FrequencyBands = NotImplemented
    reference_frequency: float = NotImplemented
    integration_time: ttime.Time = NotImplemented
