from ...core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from tudatpy.estimation.observable_models_setup import links as tlinks


class LinkEndSetup(SetupBase):

    body: str = NotImplemented
    reference_point: str = NotImplemented
