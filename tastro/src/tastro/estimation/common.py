from ..core import SettingsGenerator
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
    links as tlinks,
    model_settings as toms,
)
from tudatpy.estimation import observations as tobs
from tudatpy.estimation.observations import observations_processing as tobsp

# from .links import link_end_from_config
from ..config.estimation.observation_models.cartesian import CartesianSetup
from ..config.estimation.observation_models.doppler import (
    ClosedLoopDopplerSetup,
)
from ..config.estimation.observation_models.links import LinkEndSetup
from ..logging import log
import traceback

from typing import TYPE_CHECKING, Union
import abc

if TYPE_CHECKING:

    from ..config.estimation.observation_models.cartesian import CartesianSetup
    from ..config.estimation.observation_models.doppler import (
        ClosedLoopDopplerSetup,
        OpenLoopDopplerSetup,
    )
    from ..config.estimation.observations.interface import (
        ClosedLoopObservationsSetup,
        CartesianObservationsSetup,
        OpenLoopObservationsSetup,
    )


class ObservationModelSettingsGenerator[
    T: ("CartesianSetup", "ClosedLoopDopplerSetup", "OpenLoopDopplerSetup")
](SettingsGenerator[T], metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def observable_type_id(self) -> str: ...

    @property
    @abc.abstractmethod
    def observable_type(self) -> toms.ObservableType: ...

    @property
    def local_observation_config(self) -> Union[
        "ClosedLoopObservationsSetup",
        "CartesianObservationsSetup",
        "OpenLoopObservationsSetup",
    ]:

        observation_config = self.config.estimation.observations
        return getattr(observation_config, self.observable_type_id)

    def get_absolute_filters(
        self, station: str
    ) -> list[tobsp.ObservationFilterBase]:

        filter_setup = self.local_observation_config.filters
        absolute_filters: list[tobsp.ObservationFilterBase] = []
        for absolute_setup in filter_setup.absolute_value:

            if absolute_setup.station not in (station, "all"):
                continue

            absolute_filters.append(
                tobsp.observation_filter(
                    filter_type=tobsp.ObservationFilterType.absolute_value_filtering,
                    filter_value=absolute_setup.value,
                    filter_out=absolute_setup.filter_out,
                    use_opposite_condition=absolute_setup.use_opposite,
                )
            )

        return absolute_filters

    def get_between_epochs_filters(
        self, station: str
    ) -> list[tobsp.ObservationFilterBase]:

        filter_setup = self.local_observation_config.filters
        epoch_filters: list[tobsp.ObservationFilterBase] = []
        for epochs_setup in filter_setup.between_epochs:

            if epochs_setup.station not in (station, "all"):
                continue

            epoch_filters.append(
                tobsp.observation_filter(
                    filter_type=tobsp.ObservationFilterType.time_bounds_filtering,
                    first_filter_value=epochs_setup.first_epoch,
                    second_filter_value=epochs_setup.second_epoch,
                    filter_out=epochs_setup.filter_out,
                    use_opposite_condition=epochs_setup.use_opposite,
                )
            )

        return epoch_filters

    def get_residual_based_filters(
        self, station: str
    ) -> list[tobsp.ObservationFilterBase]:

        log.warning("Residual-based filtering only defined for open-loop")

        filter_setup = self.local_observation_config.filters
        residual_filters: list[tobsp.ObservationFilterBase] = []
        return residual_filters
        # for residual_setup in filter_setup.residual_based:

        #     if residual_setup.station not in (station, "all"):
        #         continue

        #     residual_filters.append(
        #         tobsp.observation_filter(
        #             filter_type=tobsp.ObservationFilterType.residual_filtering,
        #             filter_value=residual_setup.value,
        #             filter_out=residual_setup.filter_out,
        #             use_opposite_condition=residual_setup.use_opposite,
        #         )
        #     )

        # return residual_filters

    @abc.abstractmethod
    def apply_residual_based_filter(
        self, collection: tobs.ObservationCollection
    ) -> tobs.ObservationCollection: ...

    # def apply_residual_based_filter(
    #     self, collection: tobs.ObservationCollection
    # ) -> tobs.ObservationCollection:

    #     # Show logging information
    #     log.info(
    #         f"Filtering {self.observable_type_id} observations "
    #         "based on pre-fit residuals"
    #     )

    #     # Skip if filters not active
    #     if not self.local_observation_config.filters.present:
    #         return collection

    #     # Define parser for this type of observation
    #     type_parser = tobsp.observation_parser(self.observable_type)

    #     # Loop over link-ends with observations of this type
    #     used_receivers: set[str] = set()
    #     for link in collection.get_link_definitions_for_observables(
    #         self.observable_type
    #     ):

    #         # Get link ID for receiver
    #         receiver_id = link.link_end_id(tlinks.LinkEndType.receiver)

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

    #         # Define collection of residual filters settings
    #         residual_filter_settings = self.get_residual_based_filters(
    #             receiver_id.reference_point
    #         )
    #         if len(residual_filter_settings) == 0:
    #             continue

    #         # Apply filters to observation collection
    #         residual_filters: dict[
    #             tobsp.ObservationCollectionParser,
    #             tobsp.ObservationFilterBase,
    #         ] = {
    #             combined_parser: residual_filter
    #             for residual_filter in residual_filter_settings
    #         }
    #         collection.filter_observations(residual_filters, False)

    #     return collection

    def link_definitions(self) -> dict[str, tlinks.LinkDefinition]:

        link_definitions: dict[str, tlinks.LinkDefinition] = {}
        for link_id, link_config in self.local.link_definitions.items():

            log.debug(f"Link definition :: {link_id}")

            link_ends = {
                getattr(
                    tlinks.LinkEndType, link_end_type
                ): link_end_from_config(link_end_config, link_id)
                for link_end_type, link_end_config in link_config.__dict__.items()
            }
            link_definitions[link_id] = tlinks.LinkDefinition(link_ends)

        return link_definitions

    def light_time_correction_settings(
        self,
    ) -> list[tlight.LightTimeCorrectionSettings]:

        # Initialize container for light-time corrections
        light_time_corrections: list[tlight.LightTimeCorrectionSettings] = []
        light_time_setup = self.config.estimation.light_propagation

        # Data for tropospheric correction
        if light_time_setup.corrections.tropospheric.present:

            match light_time_setup.corrections.tropospheric.model:

                case "vmf3":

                    log.debug("VMF3 tropospheric correction")

                    light_time_corrections.append(
                        tlight.vmf3_tropospheric_light_time_correction(
                            body_with_atmosphere_name="Earth",
                            use_gradient_correction=True,
                        )
                    )

                case _:
                    raise NotImplementedError(
                        "Invalid tropospheric"
                        f" model: {light_time_setup.corrections.tropospheric.model}"
                    )

        # Data for ionospheric correction
        if light_time_setup.corrections.ionospheric.present:

            match light_time_setup.corrections.ionospheric.model:

                case "ionex":

                    log.debug("IONEX ionospheric correction")

                    light_time_corrections.append(
                        tlight.ionex_ionospheric_light_time_correction(
                            body_with_ionosphere_name="",
                            ionosphere_height=450.0,
                            first_order_delay_coefficient=40.3,
                        )
                    )

                case _:
                    raise NotImplementedError("Invalid ionospheric model")

        # Data for relativistic correction
        if light_time_setup.corrections.relativistic.present:

            match light_time_setup.corrections.relativistic.model:

                case "first_order":

                    log.debug("First order relativistic correction")

                    light_time_corrections.append(
                        tlight.first_order_relativistic_light_time_correction(
                            light_time_setup.corrections.relativistic.bodies
                        )
                    )

                case _:
                    raise NotImplementedError("Invalid relativistic model")

        return light_time_corrections

    def light_time_convergence_settings(
        self,
    ) -> tlight.LightTimeConvergenceCriteria:

        convergence_setup = self.config.estimation.light_propagation.convergence
        return tlight.light_time_convergence_settings(
            iterate_corrections=convergence_setup.iterate_corrections,
            maximum_number_of_iterations=convergence_setup.max_iterations,
            failure_handling=convergence_setup.on_failure,
        )

    @abc.abstractmethod
    def observation_collection(self) -> "tobs.ObservationCollection": ...

    @abc.abstractmethod
    def observation_model_settings(
        self, observations: "tobs.ObservationCollection"
    ) -> "list[toms.ObservationModelSettings]": ...


def link_end_from_config(
    link_end_config: LinkEndSetup, link_id: str
) -> tlinks.LinkEndId:

    return link_end_from_reference_point(
        reference_point=link_end_config.reference_point,
        body=link_end_config.body,
        link_id=link_id,
    )


def link_end_from_reference_point(
    reference_point: str, body: str, link_id: str
) -> tlinks.LinkEndId:

    # Process reference point
    match reference_point:

        case "origin":
            reference_point = body
        case "__id":
            reference_point = link_id
        case _:
            pass

    # Show debug information
    log.debug(f"Link end: {reference_point} on {body}")

    # Create link end
    if reference_point == body:
        return tlinks.body_origin_link_end_id(reference_point)
    else:
        return tlinks.body_reference_point_link_end_id(body, reference_point)
