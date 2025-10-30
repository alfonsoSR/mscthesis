from ..core import SetupBase


class AccelerationGravitationalSetup(SetupBase):

    present: bool


class AccelerationAerodynamicSetup(SetupBase):

    present: bool


class AccelerationRadiationSetup(SetupBase):

    present: bool
    occulting_bodies: list[str]


class AccelerationRelativisticSetup(SetupBase):

    present: bool
    use_karl: bool
    use_lense: bool


class AccelerationThrustSetup(SetupBase):

    present: bool


class AccelerationImpulsiveSetup(SetupBase):

    present: bool


class AccelerationEmpiricalSetup(SetupBase):

    present: bool


class PlanetAccelerationsSetup(SetupBase):

    gravitational: AccelerationGravitationalSetup
    aerodynamic: AccelerationAerodynamicSetup
    radiation: AccelerationRadiationSetup
    relativistic: AccelerationRelativisticSetup


class VehicleAccelerationsSetup(SetupBase):

    thrust: AccelerationThrustSetup
    impulsive: AccelerationImpulsiveSetup
    empirical: AccelerationEmpiricalSetup


class TargetAccelerationSetup(SetupBase):

    external: dict[str, PlanetAccelerationsSetup]
    internal: VehicleAccelerationsSetup
