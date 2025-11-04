from tastro import (
    io as nio,
    config as ncon,
    environment as nenv,
    propagation as nprop,
    estimation as nest,
)
from ...io.command_line.runner import (
    CommandLineInputRunner,
)
from tudatpy.interface import spice
from tudatpy.dynamics import (
    parameters_setup as tpars,
    simulator as tsim,
    propagation as tprop,
)
from tudatpy.dynamics.propagation_setup import acceleration as tacc
from tudatpy.astro import time_representation as ttime
from tudatpy.estimation import (
    estimation_analysis as testa,
    observations as tobs,
)
from tudatpy.estimation.observations_setup import (
    observations_simulation_settings as tosim,
)
from ...logging import log
from ..core import AnalysisManagerBase
from pathlib import Path

import typing

if typing.TYPE_CHECKING:
    from tudatpy.dynamics.environment import SystemOfBodies
    from tudatpy.dynamics.propagation_setup.propagator import (
        TranslationalStatePropagatorSettings,
    )


class SimulationManager(AnalysisManagerBase[CommandLineInputRunner]):

    def __init__(self, user_input: "CommandLineInputRunner") -> None:

        super().__init__(user_input)

        # Load configuration from command line input
        self.config = ncon.CaseSetup.from_config_file(user_input.config_file)
        self.config.perform_propagation = user_input.propagate
        self.config.perform_estimation = user_input.estimate

        return None

    def run_simulations(self) -> None:

        for source_dir in self.user_input.source_dirs:

            runner_single(
                source_dir, self.user_input.propagate, self.user_input.estimate
            )

        return None


def propagate_translational_dynamics(
    bodies: "SystemOfBodies",
    propagator: "TranslationalStatePropagatorSettings",
    config: "ncon.CaseSetup",
) -> "nio.PropagationOutput":

    log.info("Propagating dynamics")

    # Create simulator for translational dynamics
    dynamics_simulator = tsim.create_dynamics_simulator(
        bodies=bodies,
        propagator_settings=propagator,
        simulate_dynamics_on_creation=True,
    )
    # Save simulation output to file
    assert isinstance(dynamics_simulator, tsim.SingleArcSimulator)
    simulation_results = dynamics_simulator.propagation_results
    assert isinstance(simulation_results, tprop.SingleArcSimulationResults)
    results = nio.PropagationOutput.from_simulation(simulation_results, config)

    return results
    # .save_to_file(self.user_input.source_dir / "results.pkl")
    # log.info("Saved propagation results")

    # return None


def perform_estimation(
    bodies: "SystemOfBodies",
    propagator: "TranslationalStatePropagatorSettings",
    config: "ncon.CaseSetup",
) -> "nio.EstimationResults":

    log.info("Performing estimation")

    # Initialize estimation manager
    estimation_manager = nest.EstimationManager(config)

    # Generate observation collection
    observations = estimation_manager.observation_collection(bodies)

    # Create observation model
    observation_models = estimation_manager.observation_models(observations)

    # Define parameters to estimate
    parameters_to_estimate = estimation_manager.parameters_to_estimate(
        propagator=propagator,
        bodies=bodies,
    )

    # # Define parameters to estimate
    # parameter_settings = tpars.initial_states(
    #     propagator_settings=propagator,
    #     bodies=bodies,
    # )
    # parameters_to_estimate = tpars.create_parameter_set(
    #     parameter_settings=parameter_settings,
    #     bodies=bodies,
    #     propagator_settings=propagator,
    #     consider_parameters_names=[],
    # )

    # Initialize estimator
    log.info("Initializing estimator")
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
            maximum_iterations=5,
            minimum_residual_change=0.01,
            number_of_iterations_without_improvement=2,
        ),
    )
    estimation_input.define_estimation_settings(
        reintegrate_equations_on_first_iteration=True
    )

    # Perform the estimation
    log.info("Performing estimation")
    estimation_results = estimator.perform_estimation(estimation_input)  # type: ignore
    assert isinstance(estimation_results, testa.EstimationOutput)

    log.info("Saving estimation results")
    return nio.EstimationResults.from_estimation_output(
        estimation_results, observations
    )

    # results.save_to_file(self.source_dir / "estimation.pkl")
    # log.info("Saved estimation results")

    # return None


def runner_single(source_dir: Path, propagate: bool, estimate: bool) -> None:

    # Define paths to configuration and metakernel
    config_path = source_dir / "configuration.yaml"
    metakernel_path = source_dir / "metak.tm"

    # Load configuration
    config = ncon.CaseSetup.from_config_file(config_path)
    config.perform_estimation = estimate
    config.perform_propagation = propagate

    try:
        # Load metakernel
        spice.load_kernel(str(metakernel_path))

        # Create system of bodies from configuration
        bodies = nenv.system_of_bodies_from_config(config)

        # Define propagator settings
        propagator = nprop.translational_propagator_settings_from_config(
            config=config,
            bodies=bodies,
        )

        # Propagate translational dynamics
        if propagate:

            results = propagate_translational_dynamics(
                bodies=bodies, propagator=propagator, config=config
            )
            log.info("Saving propagation results")
            results.save_to_file(source_dir / "results.pkl")

        # Perform estimation
        if estimate:

            results = perform_estimation(
                bodies=bodies,
                propagator=propagator,
                config=config,
            )
            log.info("Saving estimation results")
            results.save_to_file(source_dir / "estimation.pkl")

    finally:
        # Clear kernel pool
        spice.clear_kernels()

    return None


# def runner(user_input: CommandLineInputRunner) -> None:

#     # Initialize manager
#     manager = SimulationManager(user_input)

#     try:
#         # Load metakernel
#         spice.load_kernel(str(manager.user_input.metakernel))

#         # Create system of bodies from configuration
#         bodies = nenv.system_of_bodies_from_config(manager.config)

#         # Define propagator settings
#         propagator = nprop.translational_propagator_settings_from_config(
#             config=manager.config,
#             bodies=bodies,
#         )

#         # Run propagation if requested
#         if manager.config.perform_propagation:

#             propagation_results = propagate_translational_dynamics(
#                 bodies, propagator, manager.config
#             )
#             propagation_results.save_to_file(manager.source_dir)

#         # Run estimation if requested
#         if manager.config.perform_estimation:

#             manager.perform_estimation(bodies, propagator)

#     finally:
#         spice.clear_kernels()

#     return None
