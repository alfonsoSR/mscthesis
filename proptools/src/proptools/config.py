import yaml
from pathlib import Path
from dataclasses import dataclass
import argparse
from tudatpy.astro import time_representation as ttime


class PropagationPeriod:

    def __init__(
        self,
        start: str,
        end: str,
        step: float,
        buffer: float,
        terminate_exactly: bool,
    ) -> None:

        self.step = ttime.Time(step)
        _buffer = ttime.Time(buffer)
        self.start = ttime.DateTime.from_iso_string(
            start
        ).to_epoch_time_object()
        self.start_buffer = self.start - _buffer
        self.end = ttime.DateTime.from_iso_string(end).to_epoch_time_object()
        self.end_buffer = self.end + _buffer
        self.terminate_exactly = terminate_exactly

        return None


class BodySetup:

    def __init__(
        self, name: str, center: bool, gravity_field: str, ephemerides: str
    ) -> None:

        self.name = name
        self.is_center = center
        self.gravity_field = gravity_field
        self.ephemerides = ephemerides

        return None


@dataclass
class GravityFieldAccelerationSettings:

    use: bool
    sh_model: str
    sh_order: int
    sh_degree: int


class AccelerationSettings:

    def __init__(self, settings_per_type: dict[str, dict]) -> None:

        self.gravity = GravityFieldAccelerationSettings(
            **settings_per_type["gravity_field"]
        )

        return None


class EnvironmentSetup:

    def __init__(
        self,
        global_frame_origin: str,
        global_frame_orientation: str,
        central_body: str,
        spacecraft: str,
        interpolation_step: float,
    ) -> None:

        self.global_frame_origin = global_frame_origin
        self.global_frame_orientation = global_frame_orientation
        self.central_body = central_body
        self.spacecraft = spacecraft
        self.interpolation_step = ttime.Time(interpolation_step)

        return None


@dataclass
class IntegrationSettings:

    propagator: str
    integrator: str
    rk_coefficients: str
    rk_order_to_integrate: str


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

        # Propagator settings
        self.integration = IntegrationSettings(**config["IntegrationSettings"])

        return None
