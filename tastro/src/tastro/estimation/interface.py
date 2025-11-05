from tudatpy.dynamics import (
    environment as tenv,
    parameters_setup as tpars,
    parameters as tpar,
)
from tudatpy.dynamics.environment_setup import ephemeris as tephs
from tudatpy.dynamics.propagation_setup import propagator as tprop
from tudatpy.interface import spice
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
    model_settings as toms,
)
from tudatpy.estimation import (
    observations as tobs,
)
from ..config import CaseSetup
from ..logging import log
from .cartesian import CartesianSettingsGenerator
from .closed_loop import ClosedLoopSettingsGenerator
from .open_loop import OpenLoopSettingsGenerator
from .common import ObservationModelSettingsGenerator
from typing import Type


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

        # # Closed-loop observations
        # if self.closed_loop_present:
        #     observation_collections.append(
        #         self.closed_loop_generator.observation_collection(bodies)
        #     )

        # # Cartesian observations
        # if self.cartesian_present:
        #     observation_collections.append(
        #         self.cartesian_generator.observation_collection()
        #     )

        # # Open-loop observations
        # if self.open_loop_present:
        #     observation_collections.append(
        #         self.open_loop_generator.observation_collection()
        #     )

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

        # # Closed-loop Doppler
        # if self.closed_loop_present:

        #     # Generate settings for closed-loop observations
        #     observation_models += (
        #         self.closed_loop_generator.observation_model_settings(
        #             observations
        #         )
        #     )

        # # Cartesian observations
        # if self.cartesian_present:

        #     # Generate settings for cartesian observations
        #     observation_models += (
        #         self.cartesian_generator.observation_model_settings(
        #             observations
        #         )
        #     )

        # # Open-loop Doppler
        # if self.open_loop_present:

        #     # Generate settings for open-loop observations
        #     observation_models += (
        #         self.open_loop_generator.observation_model_settings(
        #             observations
        #         )
        #     )

        return observation_models

    def parameters_to_estimate(
        self, propagator: tprop.PropagatorSettings, bodies: tenv.SystemOfBodies
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

        # Drag coefficient
        if self.config.estimation.parameters.drag_coefficient:

            log.debug("Estimation of drag coefficient")

            parameters.append(
                tpars.constant_drag_coefficient(
                    body=self.config.environment.general.spacecraft
                )
            )

        # Radiation pressure coefficient
        if self.config.estimation.parameters.radiation_pressure_coefficient:

            log.debug("Estimation of radiation pressure coefficients")

            parameters.append(
                tpars.radiation_pressure_coefficient(
                    body=self.config.environment.general.spacecraft
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
