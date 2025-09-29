import yaml
from pathlib import Path
from dataclasses import dataclass
from astropy import time


class PropagationPeriod:

    def __init__(
        self,
        start: str,
        end: str,
        step: float,
        buffer: float,
        starting_point: str,
        terminate_exactly: bool,
        custom_start_epoch: str,
    ) -> None:

        _buffer = time.TimeDelta(buffer, format="sec", scale="tdb")
        _start = time.Time(start, scale="tdb")
        _end = time.Time(end, scale="tdb")

        self.step = step
        self.start = start
        self.end = end
        self.start_buffer = str((_start - _buffer).isot)
        self.end_buffer = str((_end + _buffer).isot)
        self.terminate_exactly = terminate_exactly
        self.starting_point = starting_point
        self.custom_start_epoch = custom_start_epoch

        return None


class BodySetup:

    def __init__(
        self,
        name: str,
        center: bool,
        fixed_frame: str,
        gravity_field: str,
        ephemerides: str,
        rotation: str,
        shape: str,
    ) -> None:

        self.name = name
        self.fixed_frame = fixed_frame
        self.is_center = center
        self.gravity_field = gravity_field
        self.ephemerides = ephemerides
        self.rotation = rotation
        self.shape = shape

        return None


@dataclass
class GravityFieldAccelerationSettings:

    use: bool
    sh_model: str
    sh_order: int
    sh_degree: int


@dataclass
class RelativisticAccelerationSettings:

    use: bool
    karl: bool
    lense: bool


class AccelerationSettings:

    def __init__(self, settings_per_type: dict[str, dict]) -> None:

        self.gravity = GravityFieldAccelerationSettings(
            **settings_per_type["gravity_field"]
        )

        self.relativity = RelativisticAccelerationSettings(
            **settings_per_type["relativity"]
        )

        return None


class EnvironmentSetup:

    def __init__(
        self,
        global_frame_origin: str,
        global_frame_orientation: str,
        spacecraft: str,
        interpolation_step: float,
    ) -> None:

        self.global_frame_origin = global_frame_origin
        self.global_frame_orientation = global_frame_orientation
        self.spacecraft = spacecraft
        self.interpolation_step = interpolation_step

        return None


@dataclass
class IntegrationSettings:

    propagator: str
    integrator: str
    rk_coefficients: str
    rk_order_to_integrate: str
    step_size: float
    rtol: float
    atols: list[float]
    min_step: float
    max_step: float


@dataclass
class EphemeridesSettings:

    spacecraft: str
    planets: str
    masses: str


@dataclass
class DependentVariableSettings:

    variables: dict[str, bool]


@dataclass
class PlotSettings:

    rsw_error: bool
    keplerian_elements: bool
    keplerian_difference: bool
    dependent_variables: list[str]
    orbit: bool
    show: bool
    save: bool
    cartesian_elements: bool = False


class PropSettings:

    def __init__(self, config_file: Path) -> None:

        config = yaml.safe_load(config_file.open())

        self.time = PropagationPeriod(**config["PropagationPeriod"])
        self.env = EnvironmentSetup(**config["EnvironmentSetup"])

        # Bodies in the environment
        self.bodies: dict[str, BodySetup] = {
            name: BodySetup(name=name, **body_settings)
            for name, body_settings in config["Bodies"].items()
        }

        # Find central body
        _center: BodySetup | None = None
        for item in self.bodies.values():
            if item.is_center:
                _center = item
                break
        if _center is None:
            raise ValueError("Missing central body in configuration")
        self.center = _center

        # Acceleration model setup
        self.accelerations = {
            str(target): {
                str(body): AccelerationSettings(btsettings)
                for body, btsettings in tconfig.items()
            }
            for target, tconfig in config["AccelerationSettings"].items()
        }

        # Dependent variable settings
        self.variables = DependentVariableSettings(
            variables=config["DependentVariables"]
        )

        # Propagator settings
        self.integration = IntegrationSettings(**config["IntegrationSettings"])

        # Plotting settings
        self.plots = PlotSettings(**config["Plotting"])

        return None
