from tudatpy.estimation.observations import observations_processing as tobsp
from tudatpy.estimation import observations as tobs
from tudatpy.astro import time_representation as ttime
from tudatpy.estimation.observable_models_setup import (
    model_settings as toms,
    links as tlinks,
    biases as tbias,
)
from tudatpy.estimation.observations_setup import (
    ancillary_settings as tancs,
)
from .common import (
    ObservationModelSettingsGenerator,
    link_end_from_reference_point,
)
from ..logging import log
from ..io.observations import load_open_loop_doppler_observations_from_config
import numpy as np
import spiceypy
from nastro import graphics as ng

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..config.estimation.observation_models.doppler import (
        OpenLoopDopplerSetup,
    )


RESIDUAL_FILTERING_OFFSET: float = 1e4


class OpenLoopSettingsGenerator(
    ObservationModelSettingsGenerator["OpenLoopDopplerSetup"]
):

    @property
    def observable_type_id(self) -> str:
        return "open_loop"

    @property
    def observable_type(self) -> toms.ObservableType:
        return toms.ObservableType.doppler_measured_frequency_type

    def observation_collection(self) -> tobs.ObservationCollection:

        log.info("Generating open-loop observation collection")

        # Ancillary settings for observation collection
        ancillary_settings = self.ancillary_settings()

        # Load raw observation data per station
        data_per_station = load_open_loop_doppler_observations_from_config(
            self.config
        )

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
            noise_per_uplink: dict[str, list[float]] = {}
            for coverage, uplink_station in uplink_stations.items():

                # Get epochs and observations using current uplink
                mask = (estimated_tx_epochs >= coverage[0]) * (
                    estimated_tx_epochs <= coverage[1]
                )
                current_epochs = station_data.epochs[mask]
                current_observations = station_data.observations[mask]
                current_noise = station_data.noise[mask]

                # Skip if no observations in current section
                if len(current_epochs) == 0:
                    continue

                # Update containers
                if uplink_station not in epochs_per_uplink:
                    epochs_per_uplink[uplink_station] = current_epochs.tolist()
                    observations_per_uplink[uplink_station] = (
                        current_observations.tolist()
                    )
                    noise_per_uplink[uplink_station] = current_noise.tolist()

                else:
                    epochs_per_uplink[uplink_station] += current_epochs.tolist()
                    observations_per_uplink[
                        uplink_station
                    ] += current_observations.tolist()
                    noise_per_uplink[uplink_station] += current_noise.tolist()

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
                    observable_type=self.observable_type,
                    link_definition=link_definition,
                    observations=values,
                    observation_times=epochs_per_uplink[uplink_station],
                    reference_link_end=tlinks.LinkEndType.receiver,
                    ancilliary_settings=ancillary_settings,
                )

                # Add noise from the Fdets as raw weights
                observation_subset.set_tabulated_weights(
                    np.array(noise_per_uplink[uplink_station])
                )

                # Split observation subset into scans
                splitter_time = tobsp.observation_set_splitter(
                    splitter_type=tobsp.ObservationSetSplitterType.time_interval_splitter,
                    splitter_value=ttime.Time(5.0),
                    min_number_observations=0,
                )
                # splitter_nobs = tobsp.observation_set_splitter(
                #     splitter_type=tobsp.ObservationSetSplitterType.nb_observations_splitter,
                #     splitter_value=10,
                #     min_number_observations=1,
                # )
                __collection = tobs.ObservationCollection([observation_subset])
                __collection.split_observation_sets(splitter_time)
                # __collection.split_observation_sets(splitter_nobs)

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

                # Apply absolute value filters
                for __absolute_filter in self.get_absolute_filters(station):
                    __collection.filter_observations(
                        __absolute_filter, save_filtered_observations=False
                    )
                    # observation_subset.filter_observations(__absolute_filter)

                # Apply filters between epochs
                for __epoch_filter in self.get_between_epochs_filters(station):
                    __collection.filter_observations(
                        __epoch_filter, save_filtered_observations=False
                    )
                    # observation_subset.filter_observations(__epoch_filter)

                # # Define list of filters for current station
                # station_filters: list[tobsp.ObservationFilterBase] = []
                # if "all" in filters:
                #     station_filters += filters["all"]
                # if station in filters:
                #     station_filters += filters[station]

                # station_filters: list[tobsp.ObservationFilterBase] = (
                #     self.get_absolute_filters(station)
                #     + self.get_between_epochs_filters(station)
                # )

                # # Apply filters to observations
                # for _filter in station_filters:
                #     observation_subset.filter_observations(_filter)

                # Display debug information for station
                nobs_raw = len(epochs_per_uplink[uplink_station])
                nobs_filtered = len(__collection.concatenated_observations)
                nobs_delta = nobs_raw - nobs_filtered
                log.debug(
                    f"Station {station} :: Raw {nobs_raw} :: "
                    f"Filtered {nobs_filtered} :: Delta {nobs_delta}"
                )

                # # If observation subset is empty, skip it
                # if nobs_filtered == 0:
                #     continue

                # # Redefine weights after filtering
                # raw_weights = observation_subset.weights_vector
                # __mean = np.average(raw_weights)
                # __std = np.std(raw_weights)
                # print(raw_weights.shape)
                # print(f"{uplink_station} :: {station} :: {__mean} :: {__std}")

                # Add observation set to collection contents
                observation_collection_contents += (
                    __collection.get_single_observation_sets()
                )
                # observation_collection_contents.append(observation_subset)

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

        # # Define a time bias for each link
        # bias_settings = tbias.combined_bias(
        #     [
        #         tbias.clock_induced_bias(
        #             time_bias=0.0,
        #             associated_link_end=tlinks.LinkEndType.receiver,
        #         )
        #         for link_definition in link_definitions
        #     ]
        # )

        # Define open-loop observation model for each link
        observation_models: list[toms.ObservationModelSettings] = []
        for link_definition in link_definitions:

            # Get ID of link
            _transmitter = link_definition.link_end_id(
                tlinks.LinkEndType.transmitter
            ).reference_point
            _receiver = link_definition.link_end_id(
                tlinks.LinkEndType.receiver
            ).reference_point
            log.debug(
                f"Defining open-loop model settings: "
                f"{_transmitter} :: {_receiver}"
            )

            observation_models.append(
                toms.doppler_measured_frequency(
                    link_ends=link_definition,
                    light_time_correction_settings=light_time_corrections,
                    light_time_convergence_settings=light_time_convergence,
                    bias_settings=tbias.absolute_bias(np.array([0.0])),
                )
            )
        return observation_models

    def apply_residual_based_filter(
        self, collection: tobs.ObservationCollection
    ) -> tobs.ObservationCollection:

        log.info("Filtering open-loop observations based on pre-fit residuals")

        # Get filter settings and skip if filters are not active
        filter_config = self.local_observation_config.filters
        if not self.local_observation_config.filters.present:
            return collection

        # Define parser for this type of observation
        type_parser = tobsp.observation_parser(self.observable_type)

        # Get original number of open-loop observations for logging
        nobs_original = len(collection.get_concatenated_residuals(type_parser))

        # Loop over links with open-loop observations
        used_receivers: set[str] = set()
        for link in collection.get_link_definitions_for_observables(
            self.observable_type
        ):

            # Get link ID and reference point
            receiver_id = link.link_end_id(tlinks.LinkEndType.receiver)
            reference_point = receiver_id.reference_point

            # Skip if reference point already considered
            if reference_point in used_receivers:
                continue

            # Define observation parser for current receiver
            receiver_parser = tobsp.observation_parser(
                link_end_id=(
                    receiver_id.body_name,
                    receiver_id.reference_point,
                )
            )
            combined_parser = tobsp.observation_parser(
                observation_parsers=[type_parser, receiver_parser],
                combine_conditions=True,
            )

            # Loop over residual-based filter configurations
            for residual_config in filter_config.residual_based:

                # Skip if configuration does not apply for current station
                if residual_config.station not in ("all", reference_point):
                    continue

                # Loop over observation subsets
                for subset in collection.get_single_observation_sets(
                    combined_parser
                ):

                    # Define parser for current subset
                    bounds_parser = tobsp.observation_parser(
                        time_bounds=subset.time_bounds,
                        use_opposite_condition=False,
                    )
                    subset_parser = tobsp.observation_parser(
                        observation_parsers=[
                            type_parser,
                            receiver_parser,
                            bounds_parser,
                        ],
                        combine_conditions=True,
                    )

                    # Apply residual-based filtering to subset
                    for _ in range(residual_config.number_of_iterations):

                        # Load residuals for current subset
                        subset_residuals = (
                            np.array(
                                collection.get_concatenated_residuals(
                                    subset_parser
                                )
                            ).flatten()
                            - RESIDUAL_FILTERING_OFFSET
                        )

                        # Calculate average and standard deviation
                        subset_average = np.average(subset_residuals)
                        subset_threshold = residual_config.sigmas * np.std(
                            subset_residuals
                        )

                        # Define n-sigma filter around average
                        upper_residual_filter = tobsp.observation_filter(
                            filter_type=tobsp.ObservationFilterType.residual_filtering,
                            filter_value=(
                                RESIDUAL_FILTERING_OFFSET
                                + subset_average
                                + subset_threshold
                            ),
                            filter_out=True,
                            use_opposite_condition=False,
                        )
                        lower_residual_filter = tobsp.observation_filter(
                            filter_type=tobsp.ObservationFilterType.residual_filtering,
                            filter_value=(
                                RESIDUAL_FILTERING_OFFSET
                                + subset_average
                                - subset_threshold
                            ),
                            filter_out=True,
                            use_opposite_condition=True,
                        )

                        # Apply filters
                        collection.filter_observations(
                            upper_residual_filter,
                            subset_parser,
                            save_filtered_observations=True,
                        )
                        collection.filter_observations(
                            lower_residual_filter,
                            subset_parser,
                            save_filtered_observations=True,
                        )

                        # __dates = np.array(
                        #     collection.get_concatenated_observation_times(
                        #         subset_parser
                        #     )
                        # ).flatten()
                        # subset_noise = np.array(
                        #     collection.get_concatenated_weights(subset_parser)
                        # ).flatten()

                        # __res = (
                        #     np.array(
                        #         collection.get_concatenated_residuals(
                        #             subset_parser
                        #         )
                        #     ).flatten()
                        #     - RESIDUAL_FILTERING_OFFSET
                        # )
                        # with ng.SingleAxis() as figure:

                        #     label = f"{np.average(subset_noise):.2e} {np.std(subset_noise):.2e}"

                        #     figure.line(__dates, __res, fmt=".", label=label)

                        # subset_std = np.std(subset_noise)
                        # # print(
                        # #     len(subset_noise),
                        # #     np.average(subset_noise),
                        # #     subset_std,
                        # # )

                        # collection.set_constant_weight(
                        #     weight=(1.0 / (subset_std * subset_std)),
                        #     observation_parser=subset_parser,
                        # )

                        # __dt = (
                        #     np.array(
                        #         collection.get_concatenated_observation_times(
                        #             subset_parser
                        #         )
                        #     ).flatten()
                        #     - tref
                        # ) / 3600.0
                        # fig.line(
                        #     __dt,
                        #     np.array(
                        #         collection.get_concatenated_weights(
                        #             subset_parser
                        #         )
                        #     ).flatten(),
                        #     fmt=".",
                        # )

            # Loop over subsets and apply weights
            # with ng.SingleAxis(
            #     ng.PlotSetup(canvas_title=reference_point)
            # ) as fig:
            for subset in collection.get_single_observation_sets(
                combined_parser
            ):

                # Define parser for current subset
                bounds_parser = tobsp.observation_parser(
                    time_bounds=subset.time_bounds,
                    use_opposite_condition=False,
                )
                subset_parser = tobsp.observation_parser(
                    observation_parsers=[
                        type_parser,
                        receiver_parser,
                        bounds_parser,
                    ],
                    combine_conditions=True,
                )

                noise = collection.get_concatenated_weights(subset_parser)
                noise_std = np.std(noise)
                weight = float(1.0 / (noise_std * noise_std))
                collection.set_constant_weight_per_observation_parser(
                    weights_per_observation_parser={subset_parser: weight},
                )

                # fig.line(
                #     collection.get_concatenated_observation_times(
                #         subset_parser
                #     ),
                #     noise,
                #     fmt=".",
                # )

                if np.log10(noise_std * 1e3) > 1:

                    tref = ttime.DateTime.from_iso_string(
                        "2013-12-29 07:10:09"
                    ).to_epoch_time_object()
                    t0, tend = subset.time_bounds
                    dt0 = (t0 - tref).to_float() / 3600.0
                    dtend = (tend - tref).to_float() / 3600.0

                    log.warning(
                        f"Noise above 10mHz :: {(noise_std * 1e3):.2e} "
                        f":: {reference_point} :: {dt0:.2f} :: {dtend:.2f}"
                    )

            # # Define weights from standard deviation after filtering
            # filtered_residuals = (
            #     np.array(
            #         collection.get_concatenated_residuals(combined_parser)
            #     ).flatten()
            #     - RESIDUAL_FILTERING_OFFSET
            # )

            # filtered_noise = np.array(
            #     collection.get_concatenated_weights(combined_parser)
            # ).flatten()
            # __filtered_epochs = np.array(
            #     collection.get_concatenated_observation_times(combined_parser)
            # ).flatten()

            # tref = ttime.DateTime.from_iso_string(
            #     "2013-12-29 07:10:09"
            # ).to_epoch()
            # __dt = (__filtered_epochs - tref) / 3600.0
            # with ng.DoubleAxis(figure_setup) as fig:

            #     fig.line(
            #         __dt,
            #         filtered_noise,
            #         fmt=".",
            #         label=reference_point,
            #     )
            #     fig.line(__dt, filtered_residuals, fmt=".", axis="right")

            # Log information about weights
            current_weight = np.average(
                collection.get_concatenated_weights(combined_parser)
            )
            current_std = 1.0 / np.sqrt(current_weight)
            log.debug(
                f"{reference_point} :: "
                f"STD {current_std:.2e} :: WEIGHT {current_weight:.2e}"
            )

            # Update set of considered stations
            used_receivers.add(reference_point)

        # # Loop over link-ends with observations of this type
        # for _, residual_config in enumerate(filter_config.residual_based):

        #     # Initialize register of processed receivers
        #     used_receivers: set[str] = set()

        #     # Apply filters per link
        #     for link in collection.get_link_definitions_for_observables(
        #         self.observable_type
        #     ):

        #         # Get link ID for receiver
        #         receiver_id = link.link_end_id(tlinks.LinkEndType.receiver)

        #         # Skip if filter does not apply to current station
        #         if residual_config.station not in (
        #             "all",
        #             receiver_id.reference_point,
        #         ):
        #             continue

        #         # Skip if already considered and add to set otherwise
        #         if receiver_id.reference_point in used_receivers:
        #             continue
        #         used_receivers.add(receiver_id.reference_point)

        #         # Define combined observation parser for current receiver
        #         receiver_parser = tobsp.observation_parser(
        #             link_end_id=(
        #                 receiver_id.body_name,
        #                 receiver_id.reference_point,
        #             )
        #         )
        #         combined_parser = tobsp.observation_parser(
        #             observation_parsers=[type_parser, receiver_parser],
        #             combine_conditions=True,
        #         )

        #         # Apply filter to each observation subset
        #         for subset in collection.get_single_observation_sets(
        #             combined_parser
        #         ):

        #             # Define parser for current subset
        #             bounds_parser = tobsp.observation_parser(
        #                 time_bounds=subset.time_bounds,
        #                 use_opposite_condition=False,
        #             )
        #             subset_parser = tobsp.observation_parser(
        #                 observation_parsers=[
        #                     type_parser,
        #                     receiver_parser,
        #                     bounds_parser,
        #                 ],
        #                 combine_conditions=True,
        #             )

        #             # Apply residual-based filtering to subset
        #             for _ in range(residual_config.number_of_iterations):

        #                 # Load residuals for current subset
        #                 subset_residuals = (
        #                     collection.get_concatenated_residuals(subset_parser)
        #                 )

        #                 # Average and standard deviation of subset
        #                 _residuals = (
        #                     np.array(subset_residuals).flatten()
        #                     - RESIDUAL_FILTERING_OFFSET
        #                 )
        #                 _average = np.average(_residuals)
        #                 _sigma_factor = residual_config.sigmas * np.std(
        #                     _residuals
        #                 )

        #                 # Define n-sigma filter around average
        #                 upper_residual_filter = tobsp.observation_filter(
        #                     filter_type=tobsp.ObservationFilterType.residual_filtering,
        #                     filter_value=(
        #                         RESIDUAL_FILTERING_OFFSET
        #                         + _average
        #                         + _sigma_factor
        #                     ),
        #                     filter_out=True,
        #                     use_opposite_condition=False,
        #                 )
        #                 lower_residual_filter = tobsp.observation_filter(
        #                     filter_type=tobsp.ObservationFilterType.residual_filtering,
        #                     filter_value=(
        #                         RESIDUAL_FILTERING_OFFSET
        #                         + _average
        #                         - _sigma_factor
        #                     ),
        #                     filter_out=True,
        #                     use_opposite_condition=True,
        #                 )

        #                 # Apply filters
        #                 collection.filter_observations(
        #                     upper_residual_filter,
        #                     subset_parser,
        #                     save_filtered_observations=True,
        #                 )
        #                 collection.filter_observations(
        #                     lower_residual_filter,
        #                     subset_parser,
        #                     save_filtered_observations=True,
        #                 )

        #         # Define weights from standard deviation after filtering
        #         __residuals = (
        #             collection.get_concatenated_residuals(combined_parser)
        #             - RESIDUAL_FILTERING_OFFSET
        #         )
        #         __std = np.std(__residuals)
        #         __weight = 1.0 / (__std * __std)
        #         collection.set_constant_weight(
        #             weight=__weight,
        #             observation_parser=combined_parser,
        #         )

        #         print(f"LINK ID :: {receiver_id.reference_point}")
        #         print(f"STD :: {__std:.2e}")
        #         print(f"WEIGHT :: {__weight:.2e}")

        # Remove empty observation sets from observation collection
        collection.remove_empty_observation_sets()

        # Log change in number of observations
        nobs_after = len(collection.get_concatenated_residuals(type_parser))
        difference = nobs_original - nobs_after
        log.debug(
            f"Open-loop residual filter :: Before {nobs_original} :: "
            f"After {nobs_after} :: Removed {difference}"
        )

        return collection
