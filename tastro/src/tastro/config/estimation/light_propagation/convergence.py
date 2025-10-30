from ...core import SetupBase
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
)


class LightTimeConvergence(SetupBase):

    iterate_corrections: bool
    max_iterations: int
    on_failure: tlight.LightTimeFailureHandling
