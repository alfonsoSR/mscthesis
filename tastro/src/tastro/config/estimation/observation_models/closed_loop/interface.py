from ..links import LinkEndSetup

from ....core import SetupBase
from ...light_propagation import LightTimeSetup
from tudatpy.estimation.observations_setup import ancillary_settings as tanc
from .ancillary import DopplerAncillarySettingsSetup


class DopplerLinkDefinitionSetup(SetupBase):

    transmitter: LinkEndSetup
    retransmitter: LinkEndSetup
    receiver: LinkEndSetup


class ClosedLoopDopplerSetup(SetupBase):

    present: bool
    link_definitions: dict[str, DopplerLinkDefinitionSetup]
    # light_time: LightTimeSetup
    ancillary: DopplerAncillarySettingsSetup
