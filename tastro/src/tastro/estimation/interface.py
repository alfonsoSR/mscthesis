from tudatpy.dynamics import (
    environment as tenv,
    parameters_setup as tpars,
    parameters as tpar,
)
from tudatpy.dynamics.propagation import (
    SimulationResults,
    SingleArcSimulationResults,
    SingleArcVariationalSimulationResults,
)
from tudatpy.dynamics.environment_setup import ephemeris as tephs
from tudatpy.dynamics.propagation_setup import propagator as tprop
from tudatpy.interface import spice
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
    model_settings as toms,
)
from tudatpy.estimation.observations_setup import (
    observations_simulation_settings as tosim,
)
from tudatpy.estimation import (
    observations as tobs,
    estimation_analysis as testa,
)
from ..config import CaseSetup
from ..logging import log
from .cartesian import CartesianSettingsGenerator
from .closed_loop import ClosedLoopSettingsGenerator
from .open_loop import OpenLoopSettingsGenerator
from .common import ObservationModelSettingsGenerator
from typing import Type, TYPE_CHECKING
from tudatpy.astro import time_representation as ttime
from ..io import EstimationResults, PrefitResults

if TYPE_CHECKING:
    from tudatpy.dynamics.environment import SystemOfBodies
    from tudatpy.dynamics.propagation_setup.propagator import (
        TranslationalStatePropagatorSettings,
    )


class EstimationManager:

    available_generators: dict[str, Type[ObservationModelSettingsGenerator]] = {
        "closed_loop": ClosedLoopSettingsGenerator,
        "cartesian": CartesianSettingsGenerator,
        "open_loop": OpenLoopSettingsGenerator,
    }

    def __init__(self, config: CaseSetup) -> None:

        self.config = config

        # Initialize generators and flags
        self.flags: dict[str, bool] = {}
        self.generators: dict[str, ObservationModelSettingsGenerator] = {}

        for model, generator in self.available_generators.items():

            # Get setup for current observation model
            setup = getattr(self.config.estimation.observation_models, model)

            # Update flags
            self.flags[model] = getattr(setup, "present")

            # Update generators
            self.generators[model] = generator("", setup, config)

        # # Flags
        # self.closed_loop_present = (
        #     self.config.estimation.observation_models.closed_loop.present
        # )
        # self.cartesian_present = (
        #     self.config.estimation.observation_models.cartesian.present
        # )
        # self.open_loop_present = (
        #     self.config.estimation.observation_models.open_loop.present
        # )

        # # Initialize model settings generators
        # self.cartesian_generator = CartesianSettingsGenerator(
        #     "", config.estimation.observation_models.cartesian, config
        # )
        # self.closed_loop_generator = ClosedLoopSettingsGenerator(
        #     "", config.estimation.observation_models.closed_loop, config
        # )
        # self.open_loop_generator = OpenLoopSettingsGenerator(
        #     "",
        #     config.estimation.observation_models.open_loop,
        #     config,
        # )

        return None

    def observation_collection(
        self, bodies: tenv.SystemOfBodies
    ) -> tobs.ObservationCollection:

        # Initialize container for observation collections
        observation_collections: list[tobs.ObservationCollection] = []

        # Loop over generators and update
        for model, generator in self.generators.items():

            if self.flags[model]:
                observation_collections.append(
                    generator.observation_collection()
                )

        # Return merged observation collections
        return tobs.merge_observation_collections(observation_collections)

    def observation_models(
        self, observations: tobs.ObservationCollection
    ) -> list[toms.ObservationModelSettings]:

        # Initialize container for observation models
        observation_models: list[toms.ObservationModelSettings] = []

        # Loop over generators and update
        for model, generator in self.generators.items():

            if self.flags[model]:
                observation_models += generator.observation_model_settings(
                    observations
                )

        return observation_models

    def parameters_to_estimate(
        self,
        propagator: tprop.PropagatorSettings,
        bodies: tenv.SystemOfBodies,
    ) -> tpar.EstimatableParameterSet:

        log.info("Defining parameters to estimate")

        # Initialize list of parameters to estimate
        parameters: list[tpars.EstimatableParameterSettings] = []

        # Initial states of propagated bodies
        if self.config.estimation.parameters.initial_states:

            log.debug("Estimation of initial states")

            parameters += tpars.initial_states(
                propagator_settings=propagator,
                bodies=bodies,
            )

        # Constant drag coefficient
        if self.config.estimation.parameters.drag_coefficient:

            log.debug("Estimation of constant drag coefficient")

            parameters.append(
                tpars.constant_drag_coefficient(
                    body=self.config.environment.general.spacecraft
                )
            )

        # Drag coefficient
        if self.config.estimation.parameters.arcwise_drag_coefficient:

            log.debug("Estimation of arcwise drag coefficient")

            t0 = self.config.propagation.integrator.general.custom_start_epoch
            deltas = [-9.5, -2.5, 4.5, 11.5]
            times = [t0 + ttime.Time(dt * 3600.0) for dt in deltas]

            parameters.append(
                tpars.arcwise_constant_drag_coefficient(
                    body=self.config.environment.general.spacecraft,
                    arc_initial_times=times,
                )
            )

            # parameters.append(
            #     tpars.constant_drag_coefficient(
            #         body=self.config.environment.general.spacecraft
            #     )
            # )

        # Radiation pressure coefficient
        if self.config.estimation.parameters.radiation_pressure_coefficient:

            log.debug("Estimation of radiation pressure coefficients")

            # parameters.append(
            #     tpars.radiation_pressure_coefficient(
            #         body=self.config.environment.general.spacecraft
            #     )
            # )

            parameters.append(
                tpars.radiation_pressure_target_direction_scaling(
                    target_body=self.config.environment.general.spacecraft,
                    source_body="Sun",
                )
            )

        if len(parameters) == 0:
            raise ValueError(
                "Requested estimation without parameters to estimate"
            )

        # Create parameter set
        parameter_set = tpars.create_parameter_set(
            parameter_settings=parameters,
            bodies=bodies,
            propagator_settings=propagator,
            # consider_parameters
        )

        log.info("Created parameter set")
        return parameter_set


def perform_estimation(
    bodies: "SystemOfBodies",
    propagator: "TranslationalStatePropagatorSettings",
    config: "CaseSetup",
) -> "EstimationResults":

    log.info("Performing estimation")

    # Initialize estimation manager
    estimation_manager = EstimationManager(config)

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
            minimum_residual_change=0.001,
            number_of_iterations_without_improvement=2,
        ),
    )
    estimation_input.define_estimation_settings(
        reintegrate_equations_on_first_iteration=True,
        save_state_history_per_iteration=True,
    )

    # Perform the estimation
    log.info("Performing estimation")
    estimation_results = estimator.perform_estimation(estimation_input)  # type: ignore
    assert isinstance(estimation_results, testa.EstimationOutput)

    log.info("Saving estimation results")
    return EstimationResults.from_estimator_output_and_observations(
        estimator=estimator,
        estimation_output=estimation_results,
        observations=observations,
        config=config,
    )


def perform_prefit_analysis(
    bodies: "SystemOfBodies",
    config: "CaseSetup",
) -> PrefitResults:

    log.info("Running pre-fit analysis")

    # Initialize estimation manager
    estimation_manager = EstimationManager(config)

    # Load observations and create observation models
    observations = estimation_manager.observation_collection(bodies)
    observation_models = estimation_manager.observation_models(observations)

    # Create observation simulator
    observation_simulator = tosim.create_observation_simulators(
        observation_settings=observation_models,
        bodies=bodies,
    )

    # Calculate pre-fit residuals
    log.info(f"Calculating pre-fit residuals")
    tobs.compute_residuals_and_dependent_variables(
        observation_collection=observations,
        observation_simulators=observation_simulator,
        bodies=bodies,
    )

    # Pack results in data structure and return
    return PrefitResults.from_observation_collection(observations)
