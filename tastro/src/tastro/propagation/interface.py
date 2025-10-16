from ..config import CaseSetup
from typing import TYPE_CHECKING
from .core import PropagationSettings
from tudatpy.dynamics.propagation_setup import (
    create_acceleration_models,
    propagator as tprops,
)
from tudatpy.interface import spice
import numpy as np
from ..logging import log

if TYPE_CHECKING:
    from tudatpy.dynamics.environment import SystemOfBodies


def translational_propagator_settings_from_config(
    config: "CaseSetup", bodies: "SystemOfBodies"
) -> tprops.TranslationalStatePropagatorSettings:

    log.info("Generating settings for translational propagator")

    # Initialize settings generator
    generator = PropagationSettings("", config)

    # Define integrator settings
    integrator_settings = generator.integrator_settings()

    # Define propagation start and termination conditions
    propagation_start, termination_condition = (
        generator.starting_point_and_termination_settings()
    )

    # Define acceleration settings
    acceleration_settings = generator.acceleration_settings()

    # Retrieve propagation centers and bodies to propagate from settings
    central_bodies = [config.environment.general.center]
    bodies_to_propagate = list(acceleration_settings.keys())

    # Retrieve initial state of bodies to propagate from SPICE
    initial_states_per_body = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                target_body_name=vehicle,
                observer_body_name=center,
                reference_frame_name=config.environment.general.global_frame_orientation,
                aberration_corrections="none",
                ephemeris_time=propagation_start,
            )
            for vehicle, center in zip(bodies_to_propagate, central_bodies)
        ]
    )

    # Define acceleration model
    acceleration_model = create_acceleration_models(
        body_system=bodies,
        selected_acceleration_per_body=acceleration_settings,
        bodies_to_propagate=bodies_to_propagate,
        central_bodies=central_bodies,
    )

    # Return propagator settings
    propagator_settings = tprops.translational(
        central_bodies=central_bodies,
        acceleration_models=acceleration_model,
        bodies_to_integrate=bodies_to_propagate,
        initial_states=initial_states_per_body.T,
        initial_time=propagation_start,
        integrator_settings=integrator_settings,
        termination_settings=termination_condition,
        propagator=config.propagation.integrator.general.state_representation,
    )

    console_print_settings = propagator_settings.print_settings
    # console_print_settings.enable_all_printing(0, 0)
    # console_print_settings.print_state_indices = True
    # console_print_settings.print_dependent_variable_indices = True
    # console_print_settings.print_propagation_clock_time = True
    # console_print_settings.print_termination_reason = True
    # console_print_settings.print_number_of_function_evaluations = True
    # console_print_settings.print_initial_and_final_conditions = True

    log.info("Generated settings for translational propagator")
    return propagator_settings
