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
from .common import (
    ObservationModelSettingsGenerator,
    link_end_from_reference_point,
)
from ..logging import log
from ..io.observations import load_open_loop_doppler_observations_from_config
import numpy as np
import spiceypy

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..config.estimation.observation_models.doppler import (
        OpenLoopDopplerSetup,
    )
    from ..config.environment.stations import StationSetup


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

        # Get set of stations with uplink information
        uplink_stations: dict[tuple[ttime.Time, ttime.Time], str] = {}
        for station, station_setup in self.config.environment.stations.items():

            if station_setup.uplink.present:

                # Get uplink intervals for current station
                for start, end in zip(
                    station_setup.uplink.start, station_setup.uplink.end
                ):

                    uplink_stations[(start, end)] = station

        # Define link end for spacecraft
        log.warning("TEMPORARY FIX FOR SPACECRAFT LINK END")
        spacecraft_link_end = link_end_from_reference_point("HGA", "MEX", "")

        # Define observation collection
        log.info("Filtering open-loop observations")
        observation_collection_contents: list[tobs.SingleObservationSet] = []
        for station, station_data in data_per_station.items():

            # Estimate transmission epochs from observation timestamps
            one_way_light_time = np.array(
                spiceypy.spkezr(
                    targ=self.config.environment.general.spacecraft,
                    et=station_data.epochs,
                    ref=self.config.environment.general.global_frame_orientation,
                    abcorr="CN",
                    obs="Earth",
                )[-1]
            )
            estimated_tx_epochs = np.array(
                [
                    ttime.Time(epoch)
                    for epoch in station_data.epochs - 2.0 * one_way_light_time
                ]
            )

            # Map observation epochs to active uplink stations
            epochs_per_uplink: dict[str, list[ttime.Time]] = {}
            observations_per_uplink: dict[str, list[float]] = {}
            for coverage, uplink_station in uplink_stations.items():

                # Get epochs and observations using current uplink
                mask = (estimated_tx_epochs >= coverage[0]) * (
                    estimated_tx_epochs <= coverage[1]
                )
                current_epochs = station_data.epochs[mask]
                current_observations = station_data.observations[mask]

                # Skip if no observations in current section
                if len(current_epochs) == 0:
                    continue

                # Update containers
                if uplink_station not in epochs_per_uplink:
                    epochs_per_uplink[uplink_station] = current_epochs.tolist()
                    observations_per_uplink[uplink_station] = (
                        current_observations.tolist()
                    )
                else:
                    epochs_per_uplink[uplink_station] += current_epochs.tolist()
                    observations_per_uplink[
                        uplink_station
                    ] += current_observations.tolist()

            # Define link end for current station
            station_link_end = link_end_from_reference_point(
                station, "Earth", station
            )

            # Loop over active uplink stations and create observation sets
            for uplink_station in epochs_per_uplink:

                # Create link definition
                uplink = link_end_from_reference_point(
                    uplink_station, "Earth", uplink_station
                )
                link_definition = tlinks.LinkDefinition(
                    {
                        tlinks.LinkEndType.transmitter: uplink,
                        tlinks.LinkEndType.retransmitter: spacecraft_link_end,
                        tlinks.LinkEndType.receiver: station_link_end,
                    }
                )

                # Create observation subset
                values = [
                    np.array([x])
                    for x in observations_per_uplink[uplink_station]
                ]
                observation_subset = tobs.single_observation_set(
                    observable_type=toms.ObservableType.doppler_measured_frequency_type,
                    link_definition=link_definition,
                    observations=values,
                    observation_times=epochs_per_uplink[uplink_station],
                    reference_link_end=tlinks.LinkEndType.receiver,
                    ancilliary_settings=ancillary_settings,
                )

                # # The observations are sorted, so if one is in a section, the next
                # # one can only be in that section or in one after

                # # Get link definition for current station
                # link_definition = link_definitions[station]

                # # Create observation set
                # values = [np.array([x]) for x in station_data.observations]
                # observation_set = tobs.single_observation_set(
                #     observable_type=toms.ObservableType.doppler_measured_frequency_type,
                #     link_definition=link_definition,
                #     observations=values,
                #     observation_times=station_data.epochs.tolist(),
                #     reference_link_end=tlinks.LinkEndType.receiver,
                #     ancilliary_settings=ancillary_settings,
                # )

                # Define list of filters for current station
                station_filters: list[tobsp.ObservationFilterBase] = []
                if "all" in filters:
                    station_filters += filters["all"]
                if station in filters:
                    station_filters += filters[station]

                # Apply filters to observations
                for _filter in station_filters:
                    observation_subset.filter_observations(_filter)

                # Display debug information for station
                nobs_raw = len(epochs_per_uplink[uplink_station])
                nobs_filtered = len(
                    observation_subset.concatenated_observations
                )
                nobs_delta = nobs_raw - nobs_filtered
                log.debug(
                    f"Station {station} :: Raw {nobs_raw} :: "
                    f"Filtered {nobs_filtered} :: Delta {nobs_delta}"
                )

                # Add observation set to collection contents
                observation_collection_contents.append(observation_subset)

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
