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


def update_system_of_bodies(
    config: CaseSetup, bodies: tenv.SystemOfBodies
) -> tenv.SystemOfBodies:

    log.info("Updating system of bodies for estimation")

    # Environment updates for ionospheric correction
    light_corrections = config.estimation.light_propagation.corrections
    if light_corrections.ionospheric.present:

        match light_corrections.ionospheric.model:

            case "ionex":

                log.debug("IONEX ionospheric correction")

                # Add IONEX information to ground station objects
                ionex_files = [
                    str(file) for file in light_corrections.ionospheric.sources
                ]
                tlight.set_ionosphere_model_from_ionex(
                    data_files=ionex_files,
                    bodies=bodies,
                )

            case _:
                raise NotImplementedError(
                    "Invalid ionospheric model"
                    f": {light_corrections.ionospheric.model}"
                )

    # Environment updates for tropospheric correction
    if light_corrections.tropospheric.present:

        match light_corrections.tropospheric.model:

            case "vmf3":

                log.debug("VMF3 tropospheric correction")

                # Add VMF3 information to ground stations
                vmf3_files = [
                    str(file) for file in light_corrections.tropospheric.sources
                ]
                tlight.set_vmf_troposphere_data(
                    data_files=vmf3_files,
                    bodies=bodies,
                    file_has_meteo=True,
                    file_has_gradient=True,
                    set_meteo_data=True,
                    set_troposphere_data=True,
                )

            case _:
                raise NotImplementedError(
                    "Invalid tropospheric model: "
                    + light_corrections.tropospheric.model
                )

    # Environment updates for closed-loop observations
    if config.estimation.observation_models.closed_loop.present:

        for (
            vehicle,
            vehicle_setup,
        ) in config.environment.vehicles.items():

            log.info(f"Closed-loop estimation settings for {vehicle}")

            if vehicle != "MEX":
                raise NotImplementedError("Only MEX supported")

            # Set turnaround ratio
            match vehicle_setup.systems.turnaround_ratio:

                case "default":
                    log.debug("Default turnaround ratios")
                    bodies.get(
                        vehicle
                    ).system_models.set_default_transponder_turnaround_ratio_function()
                case _:
                    raise NotImplementedError("Invalid turnaround ratio setup")

            # Set reference point for tracking
            match vehicle_setup.systems.reference_point:

                case "origin":
                    log.debug("Reference point is origin")
                    pass
                case "HGA":
                    log.debug("Reference point is HGA")
                    # Get state of HGA wrt LVI in fixed frame
                    cstate_hga_lvi_fixed = (
                        spice.get_body_cartesian_state_at_epoch(
                            target_body_name="MEX_HGA",
                            observer_body_name="MEX_SPACECRAFT",
                            reference_frame_name="MEX_SPACECRAFT",
                            aberration_corrections="none",
                            ephemeris_time=config.time.initial_epoch,
                        )
                    )

                    # Define constant ephemerides
                    hga_ephemeris_settings = tephs.constant(
                        constant_state=cstate_hga_lvi_fixed,
                        frame_origin="MEX_SPACECRAFT",
                        frame_orientation="MEX_SPACECRAFT",
                    )
                    hga_ephemerides = tephs.create_ephemeris(
                        hga_ephemeris_settings, "HGA"
                    )

                    # Update MEX with reference point
                    bodies.get(vehicle).system_models.set_reference_point(
                        reference_point="HGA",
                        ephemeris=hga_ephemerides,
                    )

                case _:
                    raise NotImplementedError(
                        "Invalid reference point settings"
                    )

    return bodies


class EstimationManager:

    def __init__(self, config: CaseSetup) -> None:

        self.config = config

        # Flags
        self.closed_loop_present = (
            self.config.estimation.observation_models.closed_loop.present
        )
        self.cartesian_present = (
            self.config.estimation.observation_models.cartesian.present
        )

        # Initialize model settings generators
        self.cartesian_generator = CartesianSettingsGenerator(
            "", config.estimation.observation_models.cartesian, config
        )
        self.closed_loop_generator = ClosedLoopSettingsGenerator(
            "", config.estimation.observation_models.closed_loop, config
        )

        return None

    def observation_collection(
        self, bodies: tenv.SystemOfBodies
    ) -> tobs.ObservationCollection:

        # Initialize container for observation collections
        observation_collections: list[tobs.ObservationCollection] = []

        # Closed-loop observations
        if self.closed_loop_present:
            observation_collections.append(
                self.closed_loop_generator.observation_collection(bodies)
            )

        # Cartesian observations
        if self.cartesian_present:
            observation_collections.append(
                self.cartesian_generator.observation_collection()
            )

        # Return merged observation collections
        return tobs.merge_observation_collections(observation_collections)

    def observation_models(
        self, observations: tobs.ObservationCollection
    ) -> list[toms.ObservationModelSettings]:

        # Initialize container for observation models
        observation_models: list[toms.ObservationModelSettings] = []

        # Closed-loop Doppler
        if self.config.estimation.observation_models.closed_loop.present:

            # Generate settings for closed-loop observations
            observation_models += (
                self.closed_loop_generator.observation_model_settings(
                    observations
                )
            )

        # Cartesian observations
        if self.config.estimation.observation_models.cartesian.present:

            # Generate settings for cartesian observations
            observation_models += (
                self.cartesian_generator.observation_model_settings(
                    observations
                )
            )

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
        return tpars.create_parameter_set(
            parameter_settings=parameters,
            bodies=bodies,
            propagator_settings=propagator,
            # consider_parameters
        )
