from ..links import LinkEndSetup
from dataclasses import dataclass
from ....core import SetupCollectionBase
from ...light_propagation import LightTimeSetup
from tudatpy.estimation.observations_setup import ancillary_settings as tanc
from .ancillary import DopplerAncillarySettingsSetup


@dataclass
class DopplerLinkDefinitionSetup(SetupCollectionBase):

    transmitter: LinkEndSetup
    retransmitter: LinkEndSetup
    receiver: LinkEndSetup


@dataclass
class ClosedLoopDopplerSetup(SetupCollectionBase):

    present: bool
    link_definitions: dict[str, DopplerLinkDefinitionSetup]
    # light_time: LightTimeSetup
    ancillary: DopplerAncillarySettingsSetup
