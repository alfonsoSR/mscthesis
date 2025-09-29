# Goal: Define propagation setup for closed-loop estimation
# Subgoal: Find a propagator that is sufficiently good for MEX flyby: ORMM - ROB
# I want an error wrt ephemerides that is consistent with the difference between
# ROB and ORMM, which is under one km over the time-span of the observations
# I will be propagating between 2013-12-27T00:00:00 and 2014-01-01T00:00:00


from tudatpy.astro import (
    frame_conversion as tframe,
    time_representation as ttime,
)
from tudatpy.math import root_finders as troot
from proptools import config as pcon, io as pio
import argparse
from pathlib import Path
from tudatpy.dynamics import (
    environment as tenv,
    environment_setup as tenvs,
    propagation as tprop,
    propagation_setup as tprops,
    simulator as tsim,
)
from tudatpy import data as tdata
from tudatpy.interface import spice
import numpy as np

Parser = argparse.ArgumentParser()
Parser.add_argument("config_file", help="Path to configuration file")


def __generate_body_settings_from_config(
    name: str,
    config: pcon.PropSettings,
) -> tenvs.BodySettings:

    # Define empty settings for body
    settings = tenvs.BodySettings()

    # Get body-specific configuration
    body_config = config.bodies[name]

    # Ephemeris settings
    match body_config.ephemerides:
        case "none":
            pass
        case "interpolated":
            settings.ephemeris_settings = tenvs.ephemeris.interpolated_spice(
                initial_time=ttime.DateTime.from_iso_string(
                    config.time.start_buffer
                ).to_epoch_time_object(),
                final_time=ttime.DateTime.from_iso_string(
                    config.time.end_buffer
                ).to_epoch_time_object(),
                time_step=ttime.Time(config.time.step),
                frame_origin=config.env.global_frame_origin,
                frame_orientation=config.env.global_frame_orientation,
            )
        case "direct":
            settings.ephemeris_settings = tenvs.ephemeris.direct_spice(
                frame_orientation=config.env.global_frame_orientation,
                frame_origin=config.env.global_frame_origin,
                body_name_to_use=name,
            )
        case _:
            raise NotImplementedError(f"Invalid body ephemerides: {name}")

    # Orientation settings
    match body_config.rotation:

        case "none":
            pass
        case "spice":
            settings.rotation_model_settings = tenvs.rotation_model.spice(
                base_frame=config.env.global_frame_orientation,
                target_frame=f"iau_{name}".upper(),
            )
        case "precise":

            # Precise rotation model only available for Mars
            if name != "Mars":
                raise ValueError(
                    f"Precise rotation model not available for {name}"
                )

            settings.rotation_model_settings = (
                tenvs.rotation_model.mars_high_accuracy(
                    base_frame=config.env.global_frame_orientation
                )
            )

        case _:
            raise NotImplementedError(f"Invalid rotation model: {name}")

    # Shape settings
    match body_config.shape:

        case "none":
            pass
        case "spherical_spice":
            settings.shape_settings = tenvs.shape.spherical_spice()
        case _:
            raise ValueError(f"Invalid shape model: {name}")

    # Gravity field settings
    match body_config.gravity_field:
        case "none":
            pass
        case "point_mass":
            settings.gravity_field_settings = tenvs.gravity_field.central_spice(
                body_config.name
            )
        case "spherical_harmonics":

            # Load additional config from acceleration settings
            gravity_settings = config.accelerations[config.env.spacecraft][
                name
            ].gravity
            if gravity_settings.sh_model == "none":
                raise ValueError(
                    f"Requested spherical harmonics gravity field for {name}"
                    " with missing sh_model"
                )

            # Get path to directory with gravity field models
            gravity_path = Path(tdata.get_gravity_models_path()).resolve()
            gravity_path /= f"{name}/{gravity_settings.sh_model}"
            if not gravity_path.exists():
                raise ValueError(f"Gravity data not found in {gravity_path}")

            settings.gravity_field_settings = (
                tenvs.gravity_field.from_file_spherical_harmonic(
                    file=str(gravity_path),
                    maximum_degree=gravity_settings.sh_degree,
                    maximum_order=gravity_settings.sh_order,
                )
            )

        case _:
            raise ValueError(f"Invalid gravity field settings for {name}")

    return settings


def system_of_bodies_from_config(
    config: pcon.PropSettings,
) -> tenv.SystemOfBodies:

    # Initialize empty environment settings
    environment_settings = tenvs.BodyListSettings(
        frame_origin=config.env.global_frame_origin,
        frame_orientation=config.env.global_frame_orientation,
    )

    # Define settings for spacecraft
    sc_settings = tenvs.BodySettings()
    environment_settings.add_settings(sc_settings, config.env.spacecraft)

    # Define settings for massive bodies
    for body in config.bodies:

        environment_settings.add_settings(
            settings_to_add=__generate_body_settings_from_config(body, config),
            body_name=body,
        )

    # Create system of bodies
    bodies = tenvs.create_system_of_bodies(environment_settings)

    return bodies


def acceleration_settings_from_config(
    config: pcon.PropSettings,
) -> dict[str, dict[str, list[tprops.acceleration.AccelerationSettings]]]:

    # Initialize dictionary for acceleration settings
    acceleration_settings: dict[
        str, dict[str, list[tprops.acceleration.AccelerationSettings]]
    ] = {}

    for target, tconfig in config.accelerations.items():

        # Initialize settings for target
        target_settings: dict[
            str, list[tprops.acceleration.AccelerationSettings]
        ] = {}

        # Define target-settings for each body
        for body, btconfig in tconfig.items():

            # Initialize list of settings for body
            bsettings: list[tprops.acceleration.AccelerationSettings] = []

            # Gravity field settings
            if btconfig.gravity.use:

                # Get model type from environment configuration
                match config.bodies[body].gravity_field:
                    case "point_mass":
                        bsettings.append(
                            tprops.acceleration.point_mass_gravity()
                        )
                    case "spherical_harmonics":
                        bsettings.append(
                            tprops.acceleration.spherical_harmonic_gravity(
                                maximum_degree=btconfig.gravity.sh_degree,
                                maximum_order=btconfig.gravity.sh_order,
                            )
                        )
                    case _:
                        raise ValueError(
                            f"Requested gravity of {body} with invalid settings"
                        )

            # Relativistic correction
            if btconfig.relativity.use:

                bsettings.append(
                    tprops.acceleration.relativistic_correction(
                        use_schwarzschild=btconfig.relativity.karl,
                        use_lense_thirring=btconfig.relativity.lense,
                    )
                )

            # Add settings to dictionary if not empty
            if len(bsettings) > 0:
                target_settings[body] = bsettings

        # Add target settings to dictionary if not empty
        if len(target_settings.keys()) > 0:
            acceleration_settings[target] = target_settings

    return acceleration_settings


def integration_settings_from_config(
    config: pcon.PropSettings,
) -> tprops.integrator.IntegratorSettings:

    # Adjust the sign of the step based on starting point for propagation
    step = ttime.Time(config.integration.step_size)
    if config.time.starting_point == "end":
        step = ttime.Time(-config.integration.step_size)

    match config.integration.integrator:

        # Runge-Kutta fixed-step
        case "rk_fixed":

            # Get coefficient set from configuration
            coefficients = getattr(
                tprops.integrator.CoefficientSets,
                config.integration.rk_coefficients,
            )

            # Get order to integrate from configuration
            rk_use_order = getattr(
                tprops.integrator.OrderToIntegrate,
                config.integration.rk_order_to_integrate,
            )

            # Define integrator settings
            integrator = tprops.integrator.runge_kutta_fixed_step(
                time_step=step,
                coefficient_set=coefficients,
                order_to_use=rk_use_order,
            )

        # Runge-Kutta variable step
        case "rk_variable":

            # Get coefficient set from configuration
            coefficients = getattr(
                tprops.integrator.CoefficientSets,
                config.integration.rk_coefficients,
            )

            # Step-size control settings
            rtols = np.ones(6, dtype=float) * config.integration.rtol
            atols = np.array(config.integration.atols, dtype=float)

            step_control = tprops.integrator.step_size_control_elementwise_matrix_tolerance(
                relative_error_tolerance=rtols,
                absolute_error_tolerance=atols,
                safety_factor=0.8,
                minimum_factor_increase=0.1,
                maximum_factor_increase=4,
            )

            # Step-size validation settings
            step_validation = tprops.integrator.step_size_validation(
                minimum_step=ttime.Time(config.integration.min_step),
                maximum_step=ttime.Time(config.integration.max_step),
            )

            # Define integrator settings
            integrator = tprops.integrator.runge_kutta_variable_step(
                initial_time_step=step,
                coefficient_set=coefficients,
                step_size_control_settings=step_control,
                step_size_validation_settings=step_validation,
                assess_termination_on_minor_steps=False,
            )

        # Any other option not implemented yet
        case _:

            raise NotImplementedError(f"Integrator not implemented")

    # # Define integrator type
    # match config.integration.integrator:
    #     case "rk_fixed":
    #         integrator_type = (
    #             tprops.integrator.AvailableIntegrators.runge_kutta_fixed_step_size_type
    #         )
    #     case "rk_variable":
    #         integrator_type = (
    #             tprops.integrator.AvailableIntegrators.runge_kutta_variable_step_size_type
    #         )
    #     case "adams_bashforth":
    #         integrator_type = (
    #             tprops.integrator.AvailableIntegrators.adams_bashforth_moulton_type
    #         )
    #     case "bulirsch":
    #         integrator_type = (
    #             tprops.integrator.AvailableIntegrators.bulirsch_stoer_type
    #         )
    #     case _:
    #         raise ValueError("Invalid integrator")

    return integrator


def starting_point_and_termination_condition_from_config(
    config: pcon.PropSettings,
) -> tuple[ttime.Time, tprops.propagator.PropagationTerminationSettings]:

    # Create time objects for the limits of the propagation interval
    initial_epoch = ttime.DateTime.from_iso_string(
        config.time.start
    ).to_epoch_time_object()
    final_epoch = ttime.DateTime.from_iso_string(
        config.time.end
    ).to_epoch_time_object()

    # Define initial epoch for propagation based on settings
    match config.time.starting_point:
        case "start":

            # Define initial epoch
            propagation_start = initial_epoch

            # Define termination condition
            termination_condition = tprops.propagator.time_termination(
                termination_time=final_epoch,
                terminate_exactly_on_final_condition=config.time.terminate_exactly,
            )

        case "end":

            # Define initial epoch
            propagation_start = final_epoch

            # Define termination condition
            termination_condition = tprops.propagator.time_termination(
                termination_time=initial_epoch,
                terminate_exactly_on_final_condition=config.time.terminate_exactly,
            )

        case "middle" | "custom":

            # Start of propagation
            if config.time.starting_point == "middle":
                _midpoint = (final_epoch - initial_epoch) / 2
                propagation_start = initial_epoch + _midpoint
            else:
                propagation_start = ttime.DateTime.from_iso_string(
                    config.time.custom_start_epoch
                ).to_epoch_time_object()

            # Termination condition
            forward_termination = tprops.propagator.time_termination(
                termination_time=final_epoch,
                terminate_exactly_on_final_condition=config.time.terminate_exactly,
            )
            backward_termination = tprops.propagator.time_termination(
                termination_time=initial_epoch,
                terminate_exactly_on_final_condition=config.time.terminate_exactly,
            )
            termination_condition = (
                tprops.propagator.non_sequential_termination(
                    forward_termination_settings=forward_termination,
                    backward_termination_settings=backward_termination,
                )
            )

        case _:
            raise ValueError("Invalid starting point for propagation")

    return propagation_start, termination_condition


def dependent_variable_settings_from_config(
    config: pcon.PropSettings,
) -> list[tprops.dependent_variable.SingleDependentVariableSaveSettings]:

    dependent_variables: list[
        tprops.dependent_variable.SingleDependentVariableSaveSettings
    ] = []

    for variable, use in config.variables.variables.items():

        if not use:
            continue

        if "_altitude" in variable:
            dependent_variables.append(
                tprops.dependent_variable.altitude(
                    body=config.env.spacecraft,
                    central_body=variable.split("_")[0],
                )
            )

        elif "_distance" in variable:
            body, target, _ = variable.split("_")
            dependent_variables.append(
                tprops.dependent_variable.relative_distance(
                    body=target, relative_body=body
                )
            )

        else:
            raise ValueError(f"Invalid dependent variable: {variable}")

    return dependent_variables


def propagator_settings_from_config(
    config: pcon.PropSettings,
    system_of_bodies: tenv.SystemOfBodies,
) -> tprops.propagator.TranslationalStatePropagatorSettings:

    # Create acceleration model
    acceleration_settings = acceleration_settings_from_config(config)
    acceleration_model = tprops.create_acceleration_models(
        body_system=system_of_bodies,
        selected_acceleration_per_body=acceleration_settings,
        bodies_to_propagate=list(acceleration_settings.keys()),
        central_bodies=[config.center.name],
    )

    # Define propagation starting point and termination condition
    propagation_start, termination_condition = (
        starting_point_and_termination_condition_from_config(config)
    )

    # Get initial state of integrated bodies from spice
    initial_states_per_target = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                target_body_name=target,
                observer_body_name=config.center.name,
                reference_frame_name=config.env.global_frame_orientation,
                aberration_corrections="none",
                ephemeris_time=propagation_start,
            )
            for target in acceleration_settings
        ]
    )

    # Check if propagator type is valid
    if not hasattr(
        tprops.propagator.TranslationalPropagatorType,
        config.integration.propagator,
    ):
        raise ValueError(
            f"Invalid propagator type: {config.integration.propagator}"
        )

    # Dependent variable settings
    dependent_variables = dependent_variable_settings_from_config(config)

    # # Get propagator type from configuration
    # match config.integration.propagator:
    #     case "cowell":
    #         propagator_type = (
    #             tprops.propagator.TranslationalPropagatorType.cowell
    #         )
    #     case _:
    #         raise ValueError("Invalid propagator type")

    # Define intergration settings from configuration
    integrator = integration_settings_from_config(config)

    return tprops.propagator.TranslationalStatePropagatorSettings(
        central_bodies=[config.center.name],
        acceleration_models=acceleration_model,
        bodies_to_propagate=list(acceleration_settings.keys()),
        initial_body_states=initial_states_per_target[0],
        initial_time=propagation_start,
        integrator_settings=integrator,
        termination_settings=termination_condition,
        propagator=getattr(
            tprops.propagator.TranslationalPropagatorType,
            config.integration.propagator,
        ),
        dependent_variable_settings=dependent_variables,
    )

    # # # Define settings for translational propagator
    # # translational_settings = tprops.propagator.translational(
    # #     central_bodies=[config.center.name],
    # #     acceleration_models=acceleration_model,
    # #     bodies_to_integrate=list(acceleration_settings.keys()),
    # #     initial_states=initial_states_per_target[0],
    # #     initial_time=propagation_start,
    # #     integrator_settings=integrator,
    # #     termination_settings=termination_condition,
    # #     propagator=propagator_type,
    # # )

    # return translational_settings


def extract_simulation_output(
    simulation: tprop.SingleArcSimulationResults, config: pcon.PropSettings
) -> pio.PropagationOutput:

    # Get ephemerides of spacecraft at propagation epochs
    reference_state_j2000 = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                target_body_name=config.env.spacecraft,
                observer_body_name=config.center.name,
                reference_frame_name=config.env.global_frame_orientation,
                aberration_corrections="NONE",
                ephemeris_time=epoch,
            )
            for epoch in simulation.state_history_time_object.keys()
        ]
    )

    # Propagation epochs as float
    propagation_epochs = np.array(list(simulation.state_history.keys()))

    # Propagated state vector (Cartesian global frame)
    propagated_state_j2000 = np.array(list(simulation.state_history.values()))

    # Get states in RSW frame (ephemerides)
    reference_state_rsw = np.zeros_like(reference_state_j2000)
    propagated_state_rsw = np.zeros_like(propagated_state_j2000)
    for idx, rstate in enumerate(reference_state_j2000):

        # Calculate rotation matrix
        rotation_matrix = tframe.inertial_to_rsw_rotation_matrix(rstate)

        # Reference state in RSW
        _refpos_rsw = rotation_matrix @ rstate[:3]
        _refvel_rsw = rotation_matrix @ rstate[3:]
        reference_state_rsw[idx] = np.array(
            [_refpos_rsw, _refvel_rsw]
        ).flatten()

        # Propagated state in RSW
        _ppos_rsw = rotation_matrix @ propagated_state_j2000[idx][:3]
        _pvel_rsw = rotation_matrix @ propagated_state_j2000[idx][3:]
        propagated_state_rsw[idx] = np.array([_ppos_rsw, _pvel_rsw]).flatten()

    # Dependent variables
    dependent_variable_values = np.array(
        [x for x in simulation.dependent_variable_history.values()]
    ).T
    if len(dependent_variable_values.shape) == 1:
        dependent_variable_values = dependent_variable_values[None, :]

    dependent_variables: dict[str, np.ndarray] = {}
    cidx: int = 0
    for name, use in config.variables.variables.items():
        if use:
            dependent_variables[name] = dependent_variable_values[cidx]
            cidx += 1

    # Save output
    return pio.PropagationOutput(
        epochs=propagation_epochs,
        cstate_j2000=propagated_state_j2000,
        rstate_j2000=reference_state_j2000,
        cstate_rsw=propagated_state_rsw,
        rstate_rsw=reference_state_rsw,
        dvars=dependent_variables,
        number_of_function_evaluations=simulation.total_number_of_function_evaluations,
        ustate=np.array(list(simulation.unprocessed_state_history.values())).T,
    )


if __name__ == "__main__":

    # Load configuration
    config_path = Path(Parser.parse_args().config_file).resolve()
    if not config_path.exists():
        raise FileNotFoundError("Configuration file does not exist")
    config = pcon.PropSettings(config_path)

    # Load metakernel (Should be next to configuration)
    metakernel = Path(config_path.parent / "metak.tm").resolve()
    try:
        spice.load_kernel(str(metakernel))

        # Define system of bodies from configuration
        bodies = system_of_bodies_from_config(config)

        # Define propagation settings from configuration
        propagator = propagator_settings_from_config(config, bodies)

        # Propagate equations of motion
        simulator = tsim.create_dynamics_simulator(
            bodies=bodies,
            propagator_settings=propagator,
            simulate_dynamics_on_creation=True,
        )
        assert isinstance(simulator, tsim.SingleArcSimulator)
        simulation = simulator.propagation_results
        assert isinstance(simulation, tprop.SingleArcSimulationResults)

        # Extract and save output
        output = extract_simulation_output(simulation, config)
        output.save_to_file(config_path.parent / "results.pkl")

    finally:
        spice.clear_kernels()

    # # Define system of bodies from configuration
    # bodies = system_of_bodies_from_config(config)

    # # Define propagation settings from configuration
    # propagator_settings = propagator_settings_from_config(config)

    # # Propagate equations of motion
    # simulator = tsim.create_dynamics_simulator(
    #     bodies=bodies,
    #     propagator_settings=None,
    #     simulate_dynamics_on_creation=True,
    # )
