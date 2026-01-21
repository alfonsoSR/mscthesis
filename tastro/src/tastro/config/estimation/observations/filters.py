from ...core import SetupBase

from tudatpy.astro import time_representation as ttime


class AbsoluteFilterSetup(SetupBase):

    value: float
    filter_out: bool
    use_opposite: bool
    station: str = "all"


class BetweenEpochsFilterSetup(SetupBase):

    first_epoch: ttime.Time
    second_epoch: ttime.Time
    filter_out: bool
    use_opposite: bool
    station: str = "all"


class ResidualBasedFilterSetup(SetupBase):

    sigmas: float
    positive_offset: float
    splitter_value: ttime.Time
    number_of_iterations: int
    station: str = "all"


class FiltersSetup(SetupBase):

    present: bool
    absolute_value: list[AbsoluteFilterSetup]
    between_epochs: list[BetweenEpochsFilterSetup]
    residual_based: list[ResidualBasedFilterSetup]
