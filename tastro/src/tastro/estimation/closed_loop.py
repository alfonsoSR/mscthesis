from ..config.estimation.observation_models.closed_loop import (
    ClosedLoopDopplerSetup,
)
from tudatpy.estimation import observations as tobs
from tudatpy.estimation.observable_models_setup import (
    model_settings as toms,
    links as tlinks,
)
from tudatpy.estimation.observations_setup import (
    ancillary_settings as tancs,
    observations_wrapper as towpr,
)
from tudatpy.dynamics.environment import SystemOfBodies
from .common import ObservationModelSettingsGenerator
from ..logging import log
from ..io.observations import load_doppler_observations_from_config
import numpy as np


class ClosedLoopSettingsGenerator(
    ObservationModelSettingsGenerator[ClosedLoopDopplerSetup]
):

    def observation_collection(
        self,
        bodies: SystemOfBodies,
    ) -> tobs.ObservationCollection:

        log.info("Generating closed-loop observation collection")

        # Ancillary settings for observation collection
        ancillary_settings = self.ancillary_settings()

        # Doppler link definitions
        link_definitions = self.link_definitions()

        # Load raw observation data per station
        data_per_station = load_doppler_observations_from_config(self.config)

        # Define observation collection
        observation_collection_contents: list[tobs.SingleObservationSet] = []
        for station, station_data in data_per_station.items():

            log.debug(f"Single observation set: {station}")

            # Assign frequency interpolator to the ground station object
            bodies.get("Earth").get_ground_station(
                station
            ).set_transmitting_frequency_calculator(
                station_data.uplink_interpolator
            )

            # Get link definition for current station
            link_definition = link_definitions[station]

            # # Link definition for current station
            # link_definition = doppler_link_definition_from_config(
            #     closed_loop_setup.link_definitions[station]
            # )

            # Create observation set
            values = [np.array([x]) for x in station_data.observations]
            observation_collection_contents.append(
                tobs.single_observation_set(
                    observable_type=toms.ObservableType.dsn_n_way_averaged_doppler_type,
                    link_definition=link_definition,
                    observations=values,
                    observation_times=station_data.epochs.tolist(),
                    reference_link_end=tlinks.LinkEndType.receiver,
                    ancilliary_settings=ancillary_settings,
                )
            )

        # Define observation collection
        observation_collection = tobs.ObservationCollection(
            observation_collection_contents
        )

        # Compress observations if requested
        if self.config.estimation.observations.closed_loop.compression.present:

            log.info("Compressing closed-loop observations")

            compression_ratio = int(
                self.config.estimation.observations.closed_loop.compression.integration_time.to_float()
                / self.local.ancillary.integration_time.to_float()
            )

            observation_collection = towpr.create_compressed_doppler_collection(
                original_observation_collection=observation_collection,
                compression_ratio=compression_ratio,
            )

        log.info(
            "Generated observation collection from closed-loop observations"
        )

        return observation_collection

    def ancillary_settings(
        self,
    ) -> tancs.ObservationAncilliarySimulationSettings:

        log.debug("Ancillary settings")

        return tancs.dsn_n_way_doppler_ancilliary_settings(
            frequency_bands=[
                self.local.ancillary.uplink_band,
                self.local.ancillary.downlink_band,
            ],
            reference_frequency_band=self.local.ancillary.reference_band,
            reference_frequency=self.local.ancillary.reference_frequency,
            integration_time=self.local.ancillary.integration_time,
        )

    def observation_model_settings(
        self, observations: tobs.ObservationCollection
    ) -> list[toms.ObservationModelSettings]:

        log.info("Generating settings for closed-loop observation model")

        # Get link definitions for closed-loop observable
        link_definitions = observations.get_link_definitions_for_observables(
            observable_type=toms.ObservableType.dsn_n_way_averaged_doppler_type
        )

        # Get settings for light-time corrections and convergence
        light_time_corrections = self.light_time_correction_settings()
        light_time_convergence = self.light_time_convergence_settings()

        # Define closed-loop observation model for each link
        observation_models: list[toms.ObservationModelSettings] = [
            toms.dsn_n_way_doppler_averaged(
                link_ends=link_definition,
                light_time_correction_settings=light_time_corrections,
                light_time_convergence_settings=light_time_convergence,
                subtract_doppler_signature=False,
            )
            for link_definition in link_definitions
        ]
        return observation_models
