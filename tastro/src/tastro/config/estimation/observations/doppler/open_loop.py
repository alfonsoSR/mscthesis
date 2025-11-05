from ....core import SetupBase
from ..common import ObservationSourceSetup
from ..filters import FiltersSetup


class OpenLoopObservationsSetup(SetupBase):

    sources: list[ObservationSourceSetup]
    filters: FiltersSetup
    present: bool = False
