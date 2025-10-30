from ...core import SetupBase

from tudatpy.astro import time_representation as ttime


class AbsoluteFilterSetup(SetupBase):

    value: float
    filter_out: bool
    use_opposite: bool


class BetweenEpochsFilterSetup(SetupBase):

    first_epoch: ttime.Time
    second_epoch: ttime.Time
    filter_out: bool
    use_opposite: bool


class FiltersSetup(SetupBase):

    present: bool
    absolute_value: list[AbsoluteFilterSetup]
    between_epochs: list[BetweenEpochsFilterSetup]
