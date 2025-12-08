from ..core import SetupBase
from tudatpy.astro import time_representation as ttime
import numpy as np

from tudatpy.dynamics.environment_setup import aerodynamic_coefficients as taero


class VehicleEphemerisSetup(SetupBase):

    present: bool
    model: str
    ephemeris_frame_origin: str
    ephemeris_frame_orientation: str
    interpolation_step: ttime.Time
    interpolation_buffer: ttime.Time


class VehicleRotationSetup(SetupBase):

    present: bool
    model: str
    base_frame: str
    target_frame: str


class CannonballRadiationTargetSetup(SetupBase):

    reference_area: float
    radiation_pressure_coefficient: float


class PanelReflectionLawSetup(SetupBase):

    panels: list[str]
    specular_reflectivity: float
    diffuse_reflectivity: float
    instantaneous_reradiation: bool


class PaneledRadiationTargetSetup(SetupBase):

    instantaneous_reradiation: bool = True
    present: bool = True


class VehicleRadiationTargetSetup(SetupBase):

    present: bool
    model: str
    cannonball_settings: CannonballRadiationTargetSetup
    paneled_settings: list[PanelReflectionLawSetup]

    # paneled_settings: PaneledRadiationTargetSetup


class VehicleSystemsSetup(SetupBase):

    present: bool
    turnaround_ratio: str
    reference_point: str
    mass: float


class SinglePanelShapeSettings(SetupBase):

    panel_id: str
    area: float
    frame: str
    normal_direction: np.ndarray


class VehicleShapeSetup(SetupBase):

    present: bool
    model: str
    paneled_settings: list[SinglePanelShapeSettings]


class CannonballAerodynamicsSetup(SetupBase):

    reference_area: float
    force_coefficients: list[float]
    coefficient_frame: taero.AerodynamicCoefficientFrames


class VehicleAerodynamicSetup(SetupBase):

    present: bool
    model: str
    cannonball_settings: CannonballAerodynamicsSetup


class VehicleSetup(SetupBase):

    present: bool
    ephemerides: VehicleEphemerisSetup
    rotation: VehicleRotationSetup
    radiation: VehicleRadiationTargetSetup
    systems: VehicleSystemsSetup
    shape: VehicleShapeSetup
    aerodynamics: VehicleAerodynamicSetup
    # reference_point: VehicleReferencePointSetup
