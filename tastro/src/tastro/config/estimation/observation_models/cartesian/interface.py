from ..links import LinkEndSetup
from dataclasses import dataclass
from ....core import SetupCollectionBase


@dataclass
class CartesianLinkDefinitionSetup(SetupCollectionBase):

    observer: LinkEndSetup
    observed_body: LinkEndSetup


@dataclass
class CartesianSetup(SetupCollectionBase):

    present: bool
    link_definitions: dict[str, CartesianLinkDefinitionSetup]
