from ...core import SetupBase, SetupCollectionBase
from dataclasses import dataclass
from tudatpy.astro import time_representation as ttime


class AbsoluteFilterSetup(SetupBase):

    value: float = NotImplemented
    filter_out: bool = NotImplemented
    use_opposite: bool = NotImplemented


class BetweenEpochsFilterSetup(SetupBase):

    first_epoch: ttime.Time = NotImplemented
    second_epoch: ttime.Time = NotImplemented
    filter_out: bool = NotImplemented
    use_opposite: bool = NotImplemented


@dataclass
class FiltersSetup(SetupCollectionBase):

    present: bool
    absolute_value: list[AbsoluteFilterSetup]
    between_epochs: list[BetweenEpochsFilterSetup]
