from tastro import (
    io as nio,
    config as ncon,
    environment as nenv,
    propagation as nprop,
    estimation as nest,
)
from pathlib import Path
from tudatpy.interface import spice
from tudatpy.dynamics import (
    parameters_setup as tpars,
    simulator as tsim,
    propagation as tprop,
)
from tudatpy.estimation import estimation_analysis as testa
import numpy as np
from nastro import graphics as ng, types as nt

if __name__ == "__main__":

    # Load configuration from command line input
    user_input = nio.UserInputParser().parse_args()
    config = ncon.CaseSetup.from_config_file(user_input.config_file)

    # Define path to metakernel
    metakernel = Path(user_input.config_file.parent / "metak.tm").resolve()
    try:
        # Load metakernel
        spice.load_kernel(str(metakernel))

        # Create system of bodies from configuration
        bodies = nenv.system_of_bodies_from_config(config)

        # Define propagator settings
        propagator = nprop.translational_propagator_settings_from_config(
            config=config,
            bodies=bodies,
        )

        # Run propagation if requested
        if config.propagation.present:

            # Create simulator for translational dynamics
            dynamics_simulator = tsim.create_dynamics_simulator(
                bodies=bodies,
                propagator_settings=propagator,
                simulate_dynamics_on_creation=True,
            )

            # Save simulation output to file
            assert isinstance(dynamics_simulator, tsim.SingleArcSimulator)
            simulation_results = dynamics_simulator.propagation_results
            assert isinstance(
                simulation_results, tprop.SingleArcSimulationResults
            )
            nio.PropagationOutput.from_simulation(
                simulation_results, config
            ).save_to_file(user_input.config_file.parent / "results.pkl")

        # Run estimation if requested
        if config.estimation.present:

            # Load observations and create observation models
            estimation_manager = nest.EstimationManager(config)
            observations = estimation_manager.observation_collection(bodies)
            observation_models = estimation_manager.observation_models(
                observations
            )

            # Load cartesian observations
            reference_solution = nio.PropagationOutput.from_config_file(
                user_input.config_file
            )
            epochs = reference_solution.epochs

            # Define parameters to estimate
            parameter_settings = tpars.initial_states(
                propagator_settings=propagator,
                bodies=bodies,
            )
            parameters_to_estimate = tpars.create_parameter_set(
                parameter_settings=parameter_settings,
                bodies=bodies,
                propagator_settings=propagator,
                consider_parameters_names=[],
            )

            # Create the estimator
            estimator = testa.Estimator(
                bodies=bodies,
                estimated_parameters=parameters_to_estimate,
                observation_settings=observation_models,
                propagator_settings=propagator,
            )

            # Define input for the estimator
            estimation_input = testa.EstimationInput(
                observations_and_times=observations,
                convergence_checker=testa.estimation_convergence_checker(
                    maximum_iterations=5
                ),
            )
            estimation_input.define_estimation_settings(
                reintegrate_equations_on_first_iteration=True
            )

            # Perform the estimation
            estimation_results = estimator.perform_estimation(estimation_input)  # type: ignore
            assert isinstance(estimation_results, testa.EstimationOutput)

            results = nio.EstimationResults.from_estimation_output(
                estimation_results, config
            )
            results.save_to_file(
                user_input.config_file.parent / "estimation.pkl"
            )

        # # Compare results to ephemerides
        # reference = nt.CartesianPosition(
        #     *np.array(
        #         [
        #             spice.get_body_cartesian_position_at_epoch(
        #                 "MEX", "SSB", "J2000", "none", ti
        #             )
        #             for ti in epochs
        #         ]
        #     ).T
        # )

        # # Residual history
        # residual_history = estimation_results.residual_history.T
        # epochs_float = [ti for ti in epochs]

        # with ng.Mosaic("ab;cd") as fig:

        #     x_subfig = fig.subplot()
        #     y_subfig = fig.subplot()
        #     z_subfig = fig.subplot()
        #     r_subfig = fig.subplot()

        #     for idx, residuals in enumerate(residual_history):

        #         # Separate residuals per component
        #         assert isinstance(residuals, np.ndarray)
        #         residual_components = residuals.reshape(len(epochs_float), 3).T
        #         residuals_mag = np.linalg.norm(residual_components, axis=0)

        #         x_subfig.line(
        #             epochs_float,
        #             residual_components[0],
        #             fmt=".",
        #             label=f"Iteration {idx}",
        #         )

        #         y_subfig.line(
        #             epochs_float,
        #             residual_components[1],
        #             fmt=".",
        #             label=f"Iteration {idx}",
        #         )

        #         z_subfig.line(
        #             epochs_float,
        #             residual_components[2],
        #             fmt=".",
        #             label=f"Iteration {idx}",
        #         )

        #     r_subfig.line(
        #         epochs_float,
        #         np.linalg.norm(
        #             estimation_results.final_residuals.reshape(
        #                 len(epochs_float), 3
        #             ),
        #             axis=-1,
        #         ),
        #         fmt=".",
        #         label=f"Best iteration",
        #     )

        # with ng.ParasiteAxis() as fig:

        #     fig.line(
        #         doppler_observations.concatenated_times,
        #         cartesian[:, 0],
        #         fmt=".",
        #         axis="left",
        #     )
        #     fig.line(
        #         doppler_observations.concatenated_times,
        #         cartesian[:, 1],
        #         fmt=".",
        #         axis="right",
        #     )
        #     fig.line(
        #         doppler_observations.concatenated_times,
        #         cartesian[:, 2],
        #         fmt=".",
        #         axis="parasite",
        #     )
        # fig.line(
        #     doppler_observations.concatenated_times,
        #     doppler_observations.concatenated_observations,
        #     fmt=".",
        #     axis="right",
        # )

        # check = np.array(
        #     [
        #         spice.get_body_cartesian_position_at_epoch(
        #             "MEX", "Earth", "J2000", "LT", ti
        #         )
        #         for ti in epochs
        #     ]
        # )
        # print(ref.shape)
        # print(eps.shape)
        # print(
        #     len(
        #         doppler_observations.get_concatenated_observation_times_objects()
        #     )
        # )
        # print(check.shape)

        # cstate = nt.CartesianState(
        #     *np.array(
        #         [
        #             bodies.get("MEX").ephemeris.cartesian_state(eti)
        #             for eti in epochs
        #         ]
        #     ).T
        # ).to_keplerian(bodies.get("Mars").gravitational_parameter)
        # epochs_float = np.array([ti.to_float() for ti in epochs])

        # with ng.PlotKeplerianState() as fig:

        #     fig.add_state(epochs_float, cstate)

        # # Add tabulated ephemerides to vehicle
        # source_file = user_input.config_file.parent / "state_history.pkl"
        # with source_file.open("rb") as buffer:
        #     state_history = pickle.load(buffer)
        # ephemeris_settings = tephs.tabulated(
        #     body_state_history=state_history,
        #     frame_origin=config.environment.general.global_frame_origin,
        #     frame_orientation=config.environment.general.global_frame_orientation,
        # )
        # bodies.get(config.environment.general.spacecraft).ephemeris = (
        #     tephs.create_ephemeris(
        #         ephemeris_settings, config.environment.general.spacecraft
        #     )
        # )

        # # Define link ends
        # link_ends = {
        #     tslinks.LinkEndType.observed_body: tslinks.body_origin_link_end_id(
        #         config.environment.general.spacecraft
        #     ),
        #     tslinks.LinkEndType.observer: tslinks.body_origin_link_end_id(
        #         "Earth"
        #     ),
        # }
        # link_definition = tslinks.LinkDefinition(link_ends)

        # observation_simulation_settings = tsosim.tabulated_simulation_settings(
        #     observable_type=tsmodel.ObservableType.position_observable_type,
        #     link_ends=link_definition,
        #     simulation_times=list(state_history.keys()),
        #     reference_link_end_type=tslinks.LinkEndType.observer,
        # )

        # tsobs.simulate_observations(
        #     simulation_settings=[observation_simulation_settings],
        #     observation_simulators=
        # )

        # # Bias settings
        # absolute_bias_settings = tomss.biases.absolute_bias(np.array([0.0]))
        # relative_bias_settings = tomss.biases.relative_bias(np.array([0.0]))
        # bias_settings = tomss.biases.combined_bias(
        #     [absolute_bias_settings, relative_bias_settings]
        # )

        # # Define observation model
        # observation_model = tomss.model_settings.cartesian_position(
        #     link_ends=link_definition,
        #     bias_settings=bias_settings,
        # )

    finally:
        spice.clear_kernels()
