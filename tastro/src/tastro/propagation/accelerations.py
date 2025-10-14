from tudatpy.dynamics.propagation_setup import acceleration as tacs
from ..config import CaseSetup
from dataclasses import dataclass
from ..logging import log


class ExternalAccelerationSettingsGenerator:

    def __init__(self, vehicle: str, source: str, config: CaseSetup) -> None:

        self.vehicle = vehicle
        self.source = source
        self.config = config
        self.environment = self.config.environment.planets[self.source]
        self.acceleration = self.config.propagation.accelerations[
            self.vehicle
        ].external[self.source]

        return None

    def gravity_settings(self) -> tacs.AccelerationSettings:

        log.debug(f"Gravitational attraction: {self.source} on {self.vehicle}")

        match self.environment.gravity.model:

            case "spice_point_mass":

                return tacs.point_mass_gravity()

            case "spherical_harmonics":

                if self.environment.gravity.spherical_harmonics_degree is None:
                    raise ValueError("Spherical harmonics max degree not set")

                if self.environment.gravity.spherical_harmonics_order is None:
                    raise ValueError("Spherical harmonics max order not set")

                return tacs.spherical_harmonic_gravity(
                    maximum_degree=self.environment.gravity.spherical_harmonics_degree,
                    maximum_order=self.environment.gravity.spherical_harmonics_order,
                )

            case _:
                raise NotImplementedError(
                    f"Invalid gravity field model: {self.vehicle} - {self.source}"
                )

    def relativistic_settings(self) -> tacs.AccelerationSettings:

        log.debug(f"Relativistic correction: {self.source} on {self.vehicle}")

        return tacs.relativistic_correction(
            use_schwarzschild=self.acceleration.relativistic.use_karl,
            use_lense_thirring=self.acceleration.relativistic.use_lense,
            use_de_sitter=False,
        )
