from ...core import SetupBase
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
)


class LightTimeConvergence(SetupBase):

    iterate_corrections: bool = NotImplemented
    max_iterations: int = NotImplemented
    on_failure: tlight.LightTimeFailureHandling = NotImplemented
