from ..core import SetupBase, SetupCollectionBase
from dataclasses import dataclass


class AccelerationGravitationalSetup(SetupBase):
    pass


class AccelerationAerodynamicSetup(SetupBase):
    pass


class AccelerationRadiationSetup(SetupBase):
    pass


class AccelerationRelativisticSetup(SetupBase):

    use_karl: bool = NotImplemented
    use_lense: bool = NotImplemented


class AccelerationThrustSetup(SetupBase):
    pass


class AccelerationImpulsiveSetup(SetupBase):
    pass


class AccelerationEmpiricalSetup(SetupBase):
    pass


@dataclass
class PlanetAccelerationsSetup(SetupCollectionBase):

    gravitational: AccelerationGravitationalSetup
    aerodynamic: AccelerationAerodynamicSetup
    radiation: AccelerationRadiationSetup
    relativistic: AccelerationRelativisticSetup


@dataclass
class VehicleAccelerationsSetup(SetupCollectionBase):

    thrust: AccelerationThrustSetup
    impulsive: AccelerationImpulsiveSetup
    empirical: AccelerationEmpiricalSetup


@dataclass
class TargetAccelerationSetup(SetupCollectionBase):

    external: dict[str, PlanetAccelerationsSetup]
    internal: VehicleAccelerationsSetup
