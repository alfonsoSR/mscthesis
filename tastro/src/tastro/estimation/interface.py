from tudatpy.dynamics import (
    environment as tenv,
    parameters_setup as tpars,
    parameters as tpar,
)
from nastro import graphics as ng
import copy
from tudatpy.dynamics.propagation import (
    SimulationResults,
    SingleArcSimulationResults,
    SingleArcVariationalSimulationResults,
)
from pathlib import Path
from tudatpy.dynamics.environment_setup import ephemeris as tephs
from tudatpy.dynamics.propagation_setup import propagator as tprop
from tudatpy.interface import spice
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
    model_settings as toms,
    links as tlinks,
)
from tudatpy.estimation.observations import observations_processing as tobsp
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
from .open_loop import OpenLoopSettingsGenerator, RESIDUAL_FILTERING_OFFSET
from .common import ObservationModelSettingsGenerator
from typing import Type, TYPE_CHECKING
from tudatpy.astro import time_representation as ttime
from ..io import EstimationResults, PrefitResults
import traceback
from scipy.signal import find_peaks
import numpy as np
from nastro import graphics as ng

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

        # Initialize estimation time range from propagation limits
        self.__time_boundaries: tuple[ttime.Time, ttime.Time] = (
            self.config.time.initial_epoch,
            self.config.time.final_epoch,
        )

        # Initialize count of estimated biases
        self.__estimated_biases: int = 0

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

        # Merge observation collections
        observations = tobs.merge_observation_collections(
            observation_collections
        )

        # Remove empty observation subsets
        observations.remove_empty_observation_sets()

        # Update time boundaries of the estimation
        self.__time_boundaries = observations.time_bounds_time_object
        self.__links_per_observable = observations.link_ends_per_observable_type

        # Return merged observation collections
        return observations

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
        observations: tobs.ObservationCollection,
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

            # Get distance to central body over estimation
            center = self.config.environment.general.center
            target = self.config.environment.general.spacecraft
            __step = ttime.Time(60.0)  # One minute
            epochs: list[ttime.Time] = np.arange(
                self.__time_boundaries[0],
                self.__time_boundaries[1] + __step,
                __step,
            ).tolist()
            x_pos = np.array(
                [
                    spice.get_body_cartesian_position_at_epoch(
                        target_body_name=target,
                        observer_body_name=center,
                        reference_frame_name="J2000",
                        aberration_corrections="None",
                        ephemeris_time=epoch,
                    )
                    for epoch in epochs
                ]
            )
            center_distance = np.linalg.norm(x_pos, axis=-1)

            # Find peaks of center distance and associated epochs
            peak_indices, _ = find_peaks(center_distance)
            cd_epochs: list[ttime.Time] = [epochs[idx] for idx in peak_indices]

            # # Add initial or final epochs based on direction of propagation
            # match self.config.propagation.integrator.general.starting_point:

            #     # Add arc between t0 and first apoapsis
            #     case "start":
            #         cd_epochs = [self.__time_boundaries[0]] + cd_epochs

            #     # Add arc between last apoapsis and tend
            #     case "end":
            #         cd_epochs = cd_epochs + [self.__time_boundaries[-1]]

            #     # Special cases when propagating backwards and forwards
            #     # TODO: Ensure arc between start and closest limit
            #     case _:
            #         log.warning(
            #             "Potentially incorrect handling of CD per orbit"
            #         )

            # Log estimated parameters
            for epoch in cd_epochs:
                epoch_str = ttime.DateTime.from_epoch_time_object(
                    epoch
                ).to_iso_string(add_T=True, number_of_digits_seconds=0)
                log.debug(f"Estimation of drag coefficient at {epoch_str}")

            # Define estimation of C_D at epochs
            parameters.append(
                tpars.arcwise_constant_drag_coefficient(
                    body=self.config.environment.general.spacecraft,
                    arc_initial_times=cd_epochs,
                )
            )

            # parameters.append(
            #     tpars.constant_drag_coefficient(
            #         body=self.config.environment.general.spacecraft
            #     )
            # )

        # Radiation pressure coefficient (Sun direction)
        if self.config.estimation.parameters.radiation_pressure_coefficient:

            spacecraft = self.config.environment.general.spacecraft
            match self.config.environment.vehicles[spacecraft].radiation.model:
                case "cannonball":

                    log.debug("Radiation pressure coefficient: cannonball")
                    parameters.append(
                        tpars.radiation_pressure_coefficient(
                            body=self.config.environment.general.spacecraft
                        )
                    )

                case "paneled":

                    log.debug("Sun-pointing radiation pressure coefficient")
                    parameters.append(
                        tpars.radiation_pressure_target_direction_scaling(
                            target_body=self.config.environment.general.spacecraft,
                            source_body="Sun",
                        )
                    )

                case _:
                    log.fatal(traceback.extract_stack()[-1])
                    log.fatal("Invalid radiation target settings")
                    exit(1)

        # Radiation pressure coefficient (Normal to Sun)
        if self.config.estimation.parameters.k2_radiation_coefficient:

            log.debug("Estimation of k2 radiation pressure coefficient")

            parameters.append(
                tpars.radiation_pressure_target_perpendicular_direction_scaling(
                    target_body=self.config.environment.general.spacecraft,
                    source_body="Sun",
                )
            )

        # GM of Phobos
        if self.config.estimation.parameters.gm_phobos:

            log.debug("Estimation of GM of Phobos")

            parameters.append(tpars.gravitational_parameter("Phobos"))

        # Absolute bias for open-loop data
        if self.config.estimation.parameters.open_loop_biases:

            # Ensure open-loop observations are used
            if "open_loop" not in self.generators:
                log.fatal(
                    "Attempted to set up estimation of bias "
                    "without open-loop observations"
                )
                log.fatal(traceback.extract_stack()[-2])
                exit(1)

            # Retrieve link ends involved in open-loop observations
            open_loop_observable_type = self.generators[
                "open_loop"
            ].observable_type
            link_definitions = (
                observations.get_link_definitions_for_observables(
                    open_loop_observable_type
                )
            )

            # Define biases for open-loop links
            for link_definition in link_definitions:

                # Update count of estimated biases
                self.__estimated_biases += 1

                # Log information about link
                _transmitter = link_definition.link_end_id(
                    tlinks.LinkEndType.transmitter
                ).reference_point
                _receiver = link_definition.link_end_id(
                    tlinks.LinkEndType.receiver
                ).reference_point
                log.debug(
                    f"Absolute bias for link: " f"{_transmitter} :: {_receiver}"
                )

                # Update list of parameters
                parameters.append(
                    tpars.absolute_observation_bias(
                        observable_type=open_loop_observable_type,
                        link_ends=link_definition,
                    )
                )

            log.debug(
                f"Estimation {self.__estimated_biases} absolute "
                "observation biases"
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

    def prefit_based_filtering_and_weighting(
        self, collection: tobs.ObservationCollection
    ) -> None:

        # Ensure configuration specifies base directory
        if self.config.base_directory is None:
            log.fatal(
                "Requested residual-based filtering, but configuration "
                "does not specify base directory"
            )
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Update observation collection with offset pre-fit residuals
        prefits = PrefitResults.from_file(
            self.config.base_directory / "prefit_results.pkl"
        )
        collection.set_residuals(prefits.residual)

        # Apply pre-fit-based filtering
        for model, generator in self.generators.items():

            # Skip if not used
            if not self.flags[model]:
                continue

            # Perform filtering
            collection = generator.apply_residual_based_filter(collection)

        return None

    def apply_residual_based_filtering(
        self, collection: tobs.ObservationCollection
    ) -> tobs.ObservationCollection:

        # Ensure configuration specifies base directory
        if self.config.base_directory is None:
            log.fatal(
                "Requested residual-based filtering, but configuration "
                "does not specify base directory"
            )
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Load pre-fit residuals into observation collection
        prefits = PrefitResults.from_file(
            self.config.base_directory / "prefit_results.pkl"
        )
        collection.set_residuals(prefits.residual + RESIDUAL_FILTERING_OFFSET)

        # Loop over observable types
        for model, generator in self.generators.items():

            # Skip if not used
            if not self.flags[model]:
                continue

            # Perform filtering
            collection = generator.apply_residual_based_filter(collection)

        # Remove offset from residuals
        collection.set_residuals(collection.get_concatenated_residuals() * 0.0)

        return collection


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

    # Apply pre-fit-based filtering (Does nothing if not requested)
    observations = estimation_manager.apply_residual_based_filtering(
        observations
    )

    # Create observation model
    observation_models = estimation_manager.observation_models(observations)

    # Define parameters to estimate
    parameters_to_estimate = estimation_manager.parameters_to_estimate(
        propagator=propagator,
        bodies=bodies,
        observations=observations,
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

    # # Apply residual-based filtering
    # prefits = observations.get_concatenated_residuals()
    # observations.set_residuals(prefits + RESIDUAL_FILTERING_OFFSET)

    # # Loop over observable types
    # for model, generator in estimation_manager.generators.items():

    #     # Skip if not used
    #     if not estimation_manager.flags[model]:
    #         continue

    #     # Perform filtering
    #     observations = generator.apply_residual_based_filter(observations)

    # observations.set_residuals(
    #     observations.get_concatenated_residuals() - RESIDUAL_FILTERING_OFFSET
    # )

    # Pack results in data structure and return
    return PrefitResults.from_observation_collection(observations)
