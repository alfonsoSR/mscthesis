from ..core import SettingsGenerator
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
    links as tlinks,
)

# from .links import link_end_from_config
from ..config.estimation.observation_models.cartesian import CartesianSetup
from ..config.estimation.observation_models.closed_loop import (
    ClosedLoopDopplerSetup,
)
from ..config.estimation.observation_models.links import LinkEndSetup
from ..logging import log


class ObservationModelSettingsGenerator[
    T: (CartesianSetup, ClosedLoopDopplerSetup)
](SettingsGenerator[T]):

    def link_definitions(self) -> dict[str, tlinks.LinkDefinition]:

        link_definitions: dict[str, tlinks.LinkDefinition] = {}
        for link_id, link_config in self.local.link_definitions.items():

            log.debug(f"Link definition :: {link_id}")

            link_ends = {
                getattr(
                    tlinks.LinkEndType, link_end_type
                ): link_end_from_config(link_end_config)
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


def link_end_from_config(link_end_config: LinkEndSetup) -> tlinks.LinkEndId:

    log.debug(
        f"Link end: {link_end_config.reference_point} in {link_end_config.body}"
    )

    if link_end_config.reference_point == "origin":
        return tlinks.body_origin_link_end_id(link_end_config.body)
    else:
        return tlinks.body_reference_point_link_end_id(
            body_name=link_end_config.body,
            reference_point_id=link_end_config.reference_point,
        )
