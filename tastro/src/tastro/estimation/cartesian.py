from .common import ObservationModelSettingsGenerator
from ..config.estimation.observation_models.cartesian import CartesianSetup
from tudatpy.estimation.observable_models_setup import (
    model_settings as toms,
    links as tlinks,
    biases as tbias,
)
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.estimation import observations as tobs
from ..logging import log
from ..io import PropagationOutput
import numpy as np


class CartesianSettingsGenerator(
    ObservationModelSettingsGenerator[CartesianSetup]
):

    def observation_collection(self) -> tobs.ObservationCollection:

        log.info(
            "Generating observation collection from cartesian observations"
        )

        # Cartesian link definitions
        link_definitions = self.link_definitions()

        # Define single observation sets
        observation_collection_contents: list[tobs.SingleObservationSet] = []
        for source in self.config.estimation.observations.cartesian.sources:

            log.debug(f"Cartesian observations: {source.path.parent}")

            # Load propagation results from file
            raw_observations = PropagationOutput.from_file(source.path)
            cartesian_positions = [
                item for item in raw_observations.cstate_j2000[:, :3]
            ]

            # Generate single observation set
            observation_collection_contents.append(
                tobs.single_observation_set(
                    observable_type=toms.ObservableType.relative_position_observable_type,
                    link_definition=link_definitions[source.link],
                    observations=cartesian_positions,
                    observation_times=raw_observations.epochs.tolist(),
                    reference_link_end=tlinks.LinkEndType.observed_body,
                )
            )

        # Define observation collection
        observation_collection = tobs.ObservationCollection(
            observation_collection_contents
        )

        log.info("Generated observation collection from cartesian observations")

        return observation_collection

    def bias_settings(self) -> tbias.ObservationBiasSettings:

        raise NotImplementedError("Bias settings not implemented")

    def observation_model_settings(
        self, observations: tobs.ObservationCollection
    ) -> list[toms.ObservationModelSettings]:

        log.info("Generating settings for cartesian observation model")

        # Get link definitions from observation collection
        link_definitions = observations.get_link_definitions_for_observables(
            observable_type=toms.ObservableType.relative_position_observable_type
        )

        # Define cartesian observation model for each link
        observation_models: list[toms.ObservationModelSettings] = [
            toms.relative_cartesian_position(
                link_ends=link_definition,
            )
            for link_definition in link_definitions
        ]
        return observation_models
