from tudatpy.estimation.observations import observations_processing as tobsp
from tudatpy.estimation import observations as tobs
from tudatpy.astro import time_representation as ttime
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
from ..io.observations import load_open_loop_doppler_observations_from_config
import numpy as np

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..config.estimation.observation_models.doppler import (
        OpenLoopDopplerSetup,
    )


class OpenLoopSettingsGenerator(
    ObservationModelSettingsGenerator["OpenLoopDopplerSetup"]
):

    def filter_settings(self) -> dict[str, list[tobsp.ObservationFilterBase]]:

        # Initialize list of filters
        filters_per_station: dict[str, list[tobsp.ObservationFilterBase]] = {}

        # Return if filters are not set
        if not self.config.estimation.observations.open_loop.filters.present:
            return filters_per_station

        # Alias for filter setup
        filter_setup = self.config.estimation.observations.open_loop.filters

        # Absolute filters
        for absolute_setup in filter_setup.absolute_value:

            # Update container with station if not present
            if absolute_setup.station not in filters_per_station:
                filters_per_station[absolute_setup.station] = []

            # Update filters for current station
            filters_per_station[absolute_setup.station].append(
                tobsp.observation_filter(
                    filter_type=tobsp.ObservationFilterType.absolute_value_filtering,
                    filter_value=absolute_setup.value,
                    filter_out=absolute_setup.filter_out,
                    use_opposite_condition=absolute_setup.use_opposite,
                )
            )

        # Between-epochs filters
        for epochs_setup in filter_setup.between_epochs:

            # Update container with station if not present
            if epochs_setup.station not in filters_per_station:
                filters_per_station[epochs_setup.station] = []

            # Update filters for current station
            filters_per_station[epochs_setup.station].append(
                tobsp.observation_filter(
                    filter_type=tobsp.ObservationFilterType.time_bounds_filtering,
                    first_filter_value=epochs_setup.first_epoch,
                    second_filter_value=epochs_setup.second_epoch,
                    filter_out=epochs_setup.filter_out,
                    use_opposite_condition=epochs_setup.use_opposite,
                )
            )

        return filters_per_station

    def observation_collection(
        self,
    ) -> tobs.ObservationCollection:

        log.info("Generating open-loop observation collection")

        # Ancillary settings for observation collection
        ancillary_settings = self.ancillary_settings()

        # Doppler link definitions
        link_definitions = self.link_definitions()

        # Load raw observation data per station
        data_per_station = load_open_loop_doppler_observations_from_config(
            self.config
        )

        # Define settings for observation filters
        filters = self.filter_settings()

        # Define observation collection
        log.info("Filtering open-loop observations")
        observation_collection_contents: list[tobs.SingleObservationSet] = []
        for station, station_data in data_per_station.items():

            # Get link definition for current station
            link_definition = link_definitions[station]

            # Create observation set
            values = [np.array([x]) for x in station_data.observations]
            observation_set = tobs.single_observation_set(
                observable_type=toms.ObservableType.doppler_measured_frequency_type,
                link_definition=link_definition,
                observations=values,
                observation_times=station_data.epochs.tolist(),
                reference_link_end=tlinks.LinkEndType.receiver,
                ancilliary_settings=ancillary_settings,
            )

            # Define list of filters for current station
            station_filters: list[tobsp.ObservationFilterBase] = []
            if "all" in filters:
                station_filters += filters["all"]
            if station in filters:
                station_filters += filters[station]

            # Apply filters to observations
            for _filter in station_filters:
                observation_set.filter_observations(_filter)

            # Display debug information for station
            nobs_raw = len(station_data.observations)
            nobs_filtered = len(observation_set.concatenated_observations)
            nobs_delta = nobs_raw - nobs_filtered
            log.debug(
                f"Station {station} :: Raw {nobs_raw} :: "
                f"Filtered {nobs_filtered} :: Delta {nobs_delta}"
            )

            # Add observation set to collection contents
            observation_collection_contents.append(observation_set)

        # Define observation collection
        observation_collection = tobs.ObservationCollection(
            observation_collection_contents
        )

        log.info("Generated observation collection from open-loop observations")

        return observation_collection

    def ancillary_settings(
        self,
    ) -> tancs.ObservationAncilliarySimulationSettings:

        log.debug("Ancillary settings")

        return tancs.doppler_measured_frequency_ancillary_settings(
            frequency_bands=[
                self.local.ancillary.uplink_band,
                self.local.ancillary.downlink_band,
            ]
        )

        # return tancs.dsn_n_way_doppler_ancilliary_settings(
        #     frequency_bands=[
        #         self.local.ancillary.uplink_band,
        #         self.local.ancillary.downlink_band,
        #     ],
        #     reference_frequency_band=self.local.ancillary.reference_band,
        #     reference_frequency=self.local.ancillary.reference_frequency,
        #     integration_time=self.local.ancillary.integration_time,
        # )

    def observation_model_settings(
        self, observations: tobs.ObservationCollection
    ) -> list[toms.ObservationModelSettings]:

        log.info("Generating settings for open-loop observation model")

        # Get link definitions for closed-loop observable
        link_definitions = observations.get_link_definitions_for_observables(
            observable_type=toms.ObservableType.doppler_measured_frequency_type
        )

        # Get settings for light-time corrections and convergence
        light_time_corrections = self.light_time_correction_settings()
        light_time_convergence = self.light_time_convergence_settings()

        # Define closed-loop observation model for each link
        observation_models: list[toms.ObservationModelSettings] = [
            toms.doppler_measured_frequency(
                link_ends=link_definition,
                light_time_correction_settings=light_time_corrections,
                light_time_convergence_settings=light_time_convergence,
            )
            for link_definition in link_definitions
        ]
        return observation_models
