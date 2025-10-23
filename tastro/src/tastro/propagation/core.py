from ..config import CaseSetup
from tudatpy.dynamics.propagation_setup import (
    propagator as tprops,
    integrator as tigrs,
    acceleration as tacs,
)
from tudatpy.astro import time_representation as ttime
from .accelerations import ExternalAccelerationSettingsGenerator
from ..logging import log


class PropagationSettings:

    def __init__(self, name: str, config: CaseSetup) -> None:

        self.name = name
        self.config = config
        self.integrator = self.config.propagation.integrator
        self.accelerations = self.config.propagation.accelerations

        return None

    def integrator_settings(self) -> tigrs.IntegratorSettings:

        match self.integrator.general.integrator_type:

            case "rkf_fixed":

                log.debug("RKF fixed integrator")

                # Fail if not chosen in table
                if not self.integrator.rkf_fixed.present:
                    raise ValueError("Missing setup for fixed step integrator")

                # Define step size
                step = self.integrator.rkf_fixed.step_size
                if self.integrator.general.starting_point == "end":
                    step = step * -1

                # Return integrator settings
                return tigrs.runge_kutta_fixed_step(
                    time_step=step,
                    coefficient_set=self.integrator.rkf_fixed.coefficients,
                    order_to_use=tigrs.OrderToIntegrate.lower,
                    assess_termination_on_minor_steps=False,
                )

            case "rkf_variable":

                log.debug("RKF variable integrator")

                # Fail if not chosen in table
                if not self.integrator.rkf_variable.present:
                    raise ValueError(
                        "Missing setup for variable step integrator"
                    )

                # Define step-size control settings
                control = tigrs.step_size_control_elementwise_matrix_tolerance(
                    relative_error_tolerance=self.integrator.rkf_variable.rtols,
                    absolute_error_tolerance=self.integrator.rkf_variable.atols,
                    safety_factor=self.integrator.rkf_variable.safety_factor,
                    minimum_factor_increase=self.integrator.rkf_variable.min_increment,
                    maximum_factor_increase=self.integrator.rkf_variable.max_increment,
                )

                # Define step-size validation settings
                validation = tigrs.step_size_validation(
                    minimum_step=self.integrator.rkf_variable.min_step,
                    maximum_step=self.integrator.rkf_variable.max_step,
                )

                # Return integrator settings
                return tigrs.runge_kutta_variable_step(
                    initial_time_step=self.integrator.rkf_variable.initial_step,
                    coefficient_set=self.integrator.rkf_variable.coefficients,
                    step_size_control_settings=control,
                    step_size_validation_settings=validation,
                    assess_termination_on_minor_steps=False,
                )

            case _:
                raise NotImplementedError(
                    f"Invalid integrator type: {self.integrator.general.integrator_type}"
                )

    def starting_point_and_termination_settings(
        self,
    ) -> tuple[ttime.Time, tprops.PropagationTerminationSettings]:

        # Retrieve general limits for the simulation interval
        initial_epoch = self.config.time.initial_epoch
        final_epoch = self.config.time.final_epoch

        # Define starting point and termination condition
        match self.integrator.general.starting_point:

            case "start":

                log.debug("Propagation from start")

                # Set initial epoch as start
                propagation_start = initial_epoch

                # Set final epoch as termination condition
                termination_condition = tprops.time_termination(
                    termination_time=final_epoch,
                    terminate_exactly_on_final_condition=self.integrator.general.terminate_exactly,
                )

                return propagation_start, termination_condition

            case "end":

                log.debug("Propagation from end")

                # Set final epoch as start
                propagation_start = final_epoch

                # Set initial epoch as termination condition
                termination_condition = tprops.time_termination(
                    termination_time=initial_epoch,
                    terminate_exactly_on_final_condition=self.integrator.general.terminate_exactly,
                )

                return propagation_start, termination_condition

            case "middle":

                log.debug("Propagation from middle")

                # Set middle epoch as start
                propagation_start = (
                    initial_epoch + (final_epoch - initial_epoch) / 2.0
                )

            case "custom":

                log.debug("Propagation from custom epoch")

                # Set custom epoch as start
                propagation_start = self.integrator.general.custom_start_epoch
                if propagation_start is None:
                    raise ValueError("Custom propagation start not set")

            case _:
                raise NotImplementedError(
                    f"Invalid starting point: {self.integrator.general.starting_point}"
                )

        # Only reached if middle or custom: Set ends as termination
        forward_termination = tprops.time_termination(
            termination_time=final_epoch,
            terminate_exactly_on_final_condition=self.integrator.general.terminate_exactly,
        )
        backward_termination = tprops.time_termination(
            termination_time=initial_epoch,
            terminate_exactly_on_final_condition=self.integrator.general.terminate_exactly,
        )
        termination_condition = tprops.non_sequential_termination(
            forward_termination_settings=forward_termination,
            backward_termination_settings=backward_termination,
        )

        return propagation_start, termination_condition

    def acceleration_settings(
        self,
    ) -> dict[str, dict[str, list[tacs.AccelerationSettings]]]:

        # Initialize container for acceleration settings
        acceleration_settings: dict[
            str, dict[str, list[tacs.AccelerationSettings]]
        ] = {}

        # Iterate over all the vehicles
        for vehicle, vehicle_setup in self.accelerations.items():

            # Initialize container for vehicle accelerations
            vehicle_acceleration_settings: dict[
                str, list[tacs.AccelerationSettings]
            ] = {}

            # Loop over external accelerations for vehicle
            for body, body_acceleration_setup in vehicle_setup.external.items():

                if not self.config.environment.planets[body].present:
                    continue

                # Initialize container for body accelerations
                body_acceleration_settings: list[tacs.AccelerationSettings] = []

                # Initialize settings generator for body
                body_generator = ExternalAccelerationSettingsGenerator(
                    vehicle, body, self.config
                )

                # Gravitational acceleration
                if body_acceleration_setup.gravitational.present:
                    body_acceleration_settings.append(
                        body_generator.gravity_settings()
                    )

                # Relativistic correction
                if body_acceleration_setup.relativistic.present:
                    body_acceleration_settings.append(
                        body_generator.relativistic_settings()
                    )

                # Radiation pressure settings
                if body_acceleration_setup.radiation.present:
                    body_acceleration_settings.append(
                        body_generator.radiation_pressure_settings()
                    )

                # Aerodynamics settings
                if body_acceleration_setup.aerodynamic.present:
                    body_acceleration_settings.append(
                        body_generator.aerodynamic_settings()
                    )

                # Add to vehicle acceleration settings if not empty
                if len(body_acceleration_settings) > 0:
                    vehicle_acceleration_settings[body] = (
                        body_acceleration_settings
                    )

            # Add vehicle settings to container if not empty
            if len(vehicle_acceleration_settings.keys()) > 0:
                acceleration_settings[vehicle] = vehicle_acceleration_settings

        return acceleration_settings
