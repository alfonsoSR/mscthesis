from ...core import SetupBase

from tudatpy.estimation.observable_models_setup import links as tlinks


class LinkEndSetup(SetupBase):

    body: str
    reference_point: str
