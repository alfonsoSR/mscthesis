from ..core import SetupBase, SetupCollectionBase
from tudatpy.astro import time_representation as ttime
from dataclasses import dataclass
from tudatpy.dynamics.environment_setup import aerodynamic_coefficients as taero


class VehicleEphemerisSetup(SetupBase):

    model: str = NotImplemented
    ephemeris_frame_origin: str = NotImplemented
    ephemeris_frame_orientation: str = NotImplemented
    interpolation_step: ttime.Time | None = NotImplemented
    interpolation_buffer: ttime.Time | None = NotImplemented


class VehicleRotationSetup(SetupBase):

    model: str = NotImplemented
    base_frame: str = NotImplemented
    target_frame: str = NotImplemented


class CannonballRadiationTargetSetup(SetupBase):

    reference_area: float | None = NotImplemented
    radiation_pressure_coefficient: float | None = NotImplemented


@dataclass
class VehicleRadiationTargetSetup(SetupCollectionBase):

    present: bool
    model: str
    cannonball_settings: CannonballRadiationTargetSetup


class VehicleSystemsSetup(SetupBase):

    turnaround_ratio: str = NotImplemented
    reference_point: str = NotImplemented
    mass: float | None = NotImplemented


class VehicleShapeSetup(SetupBase):

    model: str = NotImplemented


class CannonballAerodynamicsSetup(SetupBase):

    reference_area: float = NotImplemented
    force_coefficients: list[float] = NotImplemented
    coefficient_frame: taero.AerodynamicCoefficientFrames = NotImplemented


@dataclass
class VehicleAerodynamicSetup(SetupCollectionBase):

    present: bool
    model: str
    cannonball_settings: CannonballAerodynamicsSetup


@dataclass
class VehicleSetup(SetupCollectionBase):

    present: bool
    ephemerides: VehicleEphemerisSetup
    rotation: VehicleRotationSetup
    radiation: VehicleRadiationTargetSetup
    systems: VehicleSystemsSetup
    shape: VehicleShapeSetup
    aerodynamics: VehicleAerodynamicSetup
    # reference_point: VehicleReferencePointSetup
