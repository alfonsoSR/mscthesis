from ..links import LinkEndSetup

from ....core import SetupBase


class CartesianLinkDefinitionSetup(SetupBase):

    observer: LinkEndSetup
    observed_body: LinkEndSetup


class CartesianSetup(SetupBase):

    present: bool
    link_definitions: dict[str, CartesianLinkDefinitionSetup]
