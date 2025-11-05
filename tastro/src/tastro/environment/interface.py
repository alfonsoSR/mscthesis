from tudatpy.dynamics import environment_setup as tenvs, environment as tenv
from ..config import CaseSetup
from .vehicle import VehicleSettings
from .planet import PlanetSettings
from .station import StationSettings
from tudatpy.estimation.observable_models_setup import (
    light_time_corrections as tlight,
)
from pathlib import Path
from ..logging import log

# from ..estimation import update_system_of_bodies
import numpy as np
import typing
from tudatpy.dynamics.environment_setup import ephemeris as tephs
from tudatpy.interface import spice
import traceback

if typing.TYPE_CHECKING:

    from tudatpy.dynamics.environment_setup import ground_station as tgs


class EnvironmentGenerator:

    def __init__(self, configuration: "CaseSetup") -> None:

        self.config = configuration

        return None

    def update_vehicle_settings(
        self, vehicle_name: str, bodies: "tenv.SystemOfBodies"
    ) -> "tenv.SystemOfBodies":

        log.info(f"Updating vehicle settings for {vehicle_name}")

        # Get vehicle setup and vehicle object
        setup = self.config.environment.vehicles[vehicle_name]
        vehicle = bodies.get(vehicle_name)

        # Fail if vehicle is not MEX
        if vehicle_name != "MEX":
            log.fatal("Only MEX is supported")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Set vehicle mass
        log.debug(f"Vehicle mass :: {setup.systems.mass}")
        vehicle.mass = setup.systems.mass

        # Set turnaround ratio
        match setup.systems.turnaround_ratio:

            case "default":
                log.debug("Default turnaround ratios")
                vehicle.system_models.set_default_transponder_turnaround_ratio_function()

            case _:
                log.fatal("Invalid turnaround ratio settings")
                log.fatal(traceback.extract_stack()[-2])
                exit(1)

        # Set reference point for tracking
        match setup.systems.reference_point:

            case "origin":
                log.debug("Reference point at origin")

            case "HGA":
                log.debug("Reference point is HGA")

                # Get state of HGA wrt LVI in fixed frame
                cstate_hga_lvi_fixed = spice.get_body_cartesian_state_at_epoch(
                    target_body_name="MEX_HGA",
                    observer_body_name="MEX_SPACECRAFT",
                    reference_frame_name="MEX_SPACECRAFT",
                    aberration_corrections="none",
                    ephemeris_time=self.config.time.initial_epoch,
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
                vehicle.system_models.set_reference_point(
                    reference_point="HGA",
                    ephemeris=hga_ephemerides,
                )

            case _:
                log.fatal("Invalid reference point settings")
                log.fatal(traceback.extract_stack()[-2])
                exit(1)

        return bodies

    def update_station_settings(
        self, station_name: str, bodies: "tenv.SystemOfBodies"
    ) -> "tenv.SystemOfBodies":

        log.info(f"Updating station settings for {station_name}")

        # Get station setup and station object
        setup = self.config.environment.stations[station_name]
        station = bodies.get("Earth").get_ground_station(station_name)

        # Update with uplink frequency calculator
        if setup.uplink.present:

            # Display information
            log.debug(f"Loading uplink information")

            # Define uplink frequency interpolator
            uplink_interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
                start_times=setup.uplink.start,
                end_times=setup.uplink.end,
                ramp_rates=setup.uplink.ramp_rate,
                start_frequency=setup.uplink.ref_freq,
            )

            # Assign uplink interpolator to station
            station.set_transmitting_frequency_calculator(uplink_interpolator)

        # Backwards compatibility check
        elif self.config.estimation.observation_models.closed_loop.present:

            try:

                # Display information
                log.debug(f"Loading uplink information [DEPRECATED]")

                # Define uplink frequency calculator
                uplink = self.config.estimation.observations.closed_loop.uplink[
                    station_name
                ]
                uplink_interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
                    start_times=uplink.start,
                    end_times=uplink.end,
                    ramp_rates=np.zeros(len(uplink.ref_freq)).tolist(),
                    start_frequency=uplink.ref_freq,
                )

                # Assign uplink interpolator to station
                station.set_transmitting_frequency_calculator(
                    uplink_interpolator
                )

            except Exception as err:
                raise err

        # No uplink frequency
        else:
            pass

        return bodies

    def update_light_propagation_settings(
        self, bodies: "tenv.SystemOfBodies"
    ) -> "tenv.SystemOfBodies":

        log.info("Updating environment with light propagation settings")

        # Get light propagation setup
        setup = self.config.estimation.light_propagation.corrections

        # Load information for tropospheric correction
        if setup.tropospheric.present:

            match setup.tropospheric.model:

                case "vmf3":

                    log.debug("Loading VMF3 tropospheric data")

                    # Get list of source files
                    source_files: list[str] = [
                        str(file) for file in setup.tropospheric.sources
                    ]

                    # Update environment
                    tlight.set_vmf_troposphere_data(
                        data_files=source_files,
                        bodies=bodies,
                        file_has_meteo=True,
                        file_has_gradient=True,
                        set_meteo_data=True,
                        set_troposphere_data=True,
                    )

                case _:
                    log.fatal(
                        f"Invalid tropospheric model :: "
                        f"{setup.tropospheric.model}"
                    )
                    log.fatal(traceback.extract_stack()[-2])
                    exit(1)

        # Load information for ionospheric correction
        if setup.ionospheric.present:

            match setup.ionospheric.model:

                case "ionex":

                    log.debug("Loading IONEX ionospheric data")

                    # Get list of source files
                    source_files: list[str] = [
                        str(file) for file in setup.ionospheric.sources
                    ]

                    # Update environment
                    tlight.set_ionosphere_model_from_ionex(
                        data_files=source_files,
                        bodies=bodies,
                    )

                case _:
                    log.fatal(
                        "Invalid ionospheric model :: "
                        f"{setup.ionospheric.model}"
                    )
                    log.fatal(traceback.extract_stack()[-2])
                    exit(1)

        return bodies

    def generate_system_of_bodies(self) -> "tenv.SystemOfBodies":

        # Initialize empty environment settings
        environment_settings = tenvs.BodyListSettings(
            frame_origin=self.config.environment.general.global_frame_origin,
            frame_orientation=self.config.environment.general.global_frame_orientation,
        )

        # Define settings for vehicles
        for (
            spacecraft,
            spacecraft_setup,
        ) in self.config.environment.vehicles.items():

            # Skip spacecraft if not selected
            if not spacecraft_setup.present:
                continue

            # Display info
            log.info(f"Generating body settings for vehicle: {spacecraft}")

            # Define empty settings for vehicle
            environment_settings.add_empty_settings(spacecraft)
            spacecraft_settings = environment_settings.get(spacecraft)

            # Initialize object to generate settings from configuration
            generator = VehicleSettings(
                spacecraft, spacecraft_setup, self.config
            )

            # Ephemeris settings
            if spacecraft_setup.ephemerides.present:
                spacecraft_settings.ephemeris_settings = (
                    generator.ephemeris_settings()
                )

            # Rotation settings
            if spacecraft_setup.rotation.present:
                spacecraft_settings.rotation_model_settings = (
                    generator.rotation_settings()
                )

            # Shape settings
            if spacecraft_setup.shape.present:
                spacecraft_settings.vehicle_shape_settings = (
                    generator.shape_settings()
                )

            # Radiation target settings
            if spacecraft_setup.radiation.present:
                spacecraft_settings.radiation_pressure_target_settings = (
                    generator.radiation_target_settings()
                )

            # Aerodynamic interface
            if spacecraft_setup.aerodynamics.present:
                spacecraft_settings.aerodynamic_coefficient_settings = (
                    generator.aerodynamic_settings()
                )

        # Define settings for planets
        for planet, planet_setup in self.config.environment.planets.items():

            if not planet_setup.present:
                continue

            log.info(f"Generating body settings for planet: {planet}")

            # Define empty settings for planet
            environment_settings.add_empty_settings(planet)
            planet_settings = environment_settings.get(planet)

            # Initialize object to generate settings from configuration
            generator = PlanetSettings(planet, planet_setup, self.config)

            # Ephemeris settings
            if planet_setup.ephemerides.present:
                planet_settings.ephemeris_settings = (
                    generator.ephemeris_settings()
                )

            # Rotation settings
            if planet_setup.rotation.present:
                planet_settings.rotation_model_settings = (
                    generator.rotation_settings()
                )

            # Shape settings
            if planet_setup.shape.present:
                planet_settings.shape_settings = generator.shape_settings()

            # Gravity settings
            if planet_setup.gravity.present:
                planet_settings.gravity_field_settings = (
                    generator.gravity_settings()
                )

            # Gravity field variation settings
            if planet_setup.tides.present:
                planet_settings.gravity_field_variation_settings = (
                    generator.gravity_field_variation_settings()
                )

            # Radiation source settings
            if planet_setup.radiation.present:
                planet_settings.radiation_source_settings = (
                    generator.radiation_source_settings()
                )

            # Atmosphere settings
            if planet_setup.atmosphere.present:
                planet_settings.atmosphere_settings = (
                    generator.atmosphere_settings()
                )

        # Define settings for ground stations
        if "Earth" in self.config.environment.planets:

            # Define list of ground station settings
            ground_station_settings: list[tgs.GroundStationSettings] = []
            for (
                station,
                station_setup,
            ) in self.config.environment.stations.items():

                # Skip station if not selected
                if not station_setup.present:
                    continue

                # Display information and initialize settings generator
                log.info(f"Generating ground station settings: {station}")
                generator = StationSettings(station, station_setup, self.config)

                # Define settings
                station_settings = generator.station_settings()
                ground_station_settings.append(station_settings)

            # Add settings to object
            environment_settings.get("Earth").ground_station_settings = (
                ground_station_settings
            )

        # Create system of bodies
        log.info("Creating system of bodies")
        return tenvs.create_system_of_bodies(environment_settings)

    def update_system_of_bodies(
        self, bodies: "tenv.SystemOfBodies"
    ) -> "tenv.SystemOfBodies":

        # Update settings for vehicles
        for vehicle in self.config.environment.vehicles:

            if self.config.environment.vehicles[vehicle].systems.present:

                bodies = self.update_vehicle_settings(vehicle, bodies)

        # Update settings for ground stations
        if "Earth" in bodies.list_of_bodies():

            for station in self.config.environment.stations:

                if self.config.environment.stations[station].present:

                    bodies = self.update_station_settings(station, bodies)

        # Update environment with light propagation settings
        if self.config.perform_estimation:

            bodies = self.update_light_propagation_settings(bodies)

        return bodies


def system_of_bodies_from_config(config: "CaseSetup") -> "tenv.SystemOfBodies":

    # Initialize generator
    generator = EnvironmentGenerator(config)

    # Create system of bodies
    bodies = generator.generate_system_of_bodies()

    # Update system of bodies
    bodies = generator.update_system_of_bodies(bodies)

    return bodies


#     # Initialize empty environment settings
#     environment_settings = tenvs.BodyListSettings(
#         frame_origin=config.environment.general.global_frame_origin,
#         frame_orientation=config.environment.general.global_frame_orientation,
#     )

#     # Define settings for vehicles
#     for spacecraft, spacecraft_setup in config.environment.vehicles.items():

#         # Skip spacecraft if not selected
#         if not spacecraft_setup.present:
#             continue

#         # Display info
#         log.info(f"Generating body settings for vehicle: {spacecraft}")

#         # Define empty settings for vehicle
#         environment_settings.add_empty_settings(spacecraft)
#         spacecraft_settings = environment_settings.get(spacecraft)

#         # Initialize object to generate settings from configuration
#         generator = VehicleSettings(spacecraft, spacecraft_setup, config)

#         # Ephemeris settings
#         if spacecraft_setup.ephemerides.present:
#             spacecraft_settings.ephemeris_settings = (
#                 generator.ephemeris_settings()
#             )

#         # Rotation settings
#         if spacecraft_setup.rotation.present:
#             spacecraft_settings.rotation_model_settings = (
#                 generator.rotation_settings()
#             )

#         # Shape settings
#         if spacecraft_setup.shape.present:
#             spacecraft_settings.vehicle_shape_settings = (
#                 generator.shape_settings()
#             )

#         # Radiation target settings
#         if spacecraft_setup.radiation.present:
#             spacecraft_settings.radiation_pressure_target_settings = (
#                 generator.radiation_target_settings()
#             )

#         # Aerodynamic interface
#         if spacecraft_setup.aerodynamics.present:
#             spacecraft_settings.aerodynamic_coefficient_settings = (
#                 generator.aerodynamic_settings()
#             )

#     # Define settings for planets
#     for planet, planet_setup in config.environment.planets.items():

#         if not planet_setup.present:
#             continue

#         log.info(f"Generating body settings for planet: {planet}")

#         # Define empty settings for planet
#         environment_settings.add_empty_settings(planet)
#         planet_settings = environment_settings.get(planet)

#         # Initialize object to generate settings from configuration
#         generator = PlanetSettings(planet, planet_setup, config)

#         # Ephemeris settings
#         if planet_setup.ephemerides.present:
#             planet_settings.ephemeris_settings = generator.ephemeris_settings()

#         # Rotation settings
#         if planet_setup.rotation.present:
#             planet_settings.rotation_model_settings = (
#                 generator.rotation_settings()
#             )

#         # Shape settings
#         if planet_setup.shape.present:
#             planet_settings.shape_settings = generator.shape_settings()

#         # Gravity settings
#         if planet_setup.gravity.present:
#             planet_settings.gravity_field_settings = (
#                 generator.gravity_settings()
#             )

#         # Gravity field variation settings
#         if planet_setup.tides.present:
#             planet_settings.gravity_field_variation_settings = (
#                 generator.gravity_field_variation_settings()
#             )

#         # Radiation source settings
#         if planet_setup.radiation.present:
#             planet_settings.radiation_source_settings = (
#                 generator.radiation_source_settings()
#             )

#         # Atmosphere settings
#         if planet_setup.atmosphere.present:
#             planet_settings.atmosphere_settings = (
#                 generator.atmosphere_settings()
#             )

#     # Define settings for ground stations
#     if "Earth" in config.environment.planets:

#         # Define list of ground station settings
#         ground_station_settings: list[tgs.GroundStationSettings] = []
#         for station, station_setup in config.environment.stations.items():

#             # Skip station if not selected
#             if not station_setup.present:
#                 continue

#             # Display information and initialize settings generator
#             log.info(f"Generating ground station settings: {station}")
#             generator = StationSettings(station, station_setup, config)

#             # Define settings
#             station_settings = generator.station_settings()
#             ground_station_settings.append(station_settings)

#         # Add settings to object
#         environment_settings.get("Earth").ground_station_settings = (
#             ground_station_settings
#         )

#     # Create system of bodies
#     log.info("Creating system of bodies")
#     bodies = tenvs.create_system_of_bodies(environment_settings)

#     # # Update vehicles with mass if defined
#     # for vehicle, vehicle_setup in config.environment.vehicles.items():
#     #     if vehicle_setup.systems.mass is not None:
#     #         log.debug(
#     #             f"Setting mass of {vehicle} to {vehicle_setup.systems.mass}"
#     #         )
#     #         bodies.get(vehicle).mass = vehicle_setup.systems.mass

#     # If estimation is present, update system of bodies
#     # if config.perform_estimation:
#     bodies = update_system_of_bodies(config, bodies)

#     return bodies


# def update_system_of_bodies(
#     config: CaseSetup, bodies: tenv.SystemOfBodies
# ) -> tenv.SystemOfBodies:

#     log.info("Updating system of bodies for estimation")

#     # Add uplink frequency calculators to ground stations
#     for station_id, station_setup in config.environment.stations.items():

#         # Skip if Earth not in system of bodies
#         if "Earth" not in bodies.list_of_bodies():
#             continue

#         # Skip if station is not selected
#         if not station_setup.present:
#             continue

#         # Skip if uplink not selected
#         if not station_setup.uplink.present:
#             continue

#         # Display information
#         log.debug(f"Loading uplink frequencies for {station_id}")

#         # Define uplink frequency interpolator
#         uplink_interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
#             start_times=station_setup.uplink.start,
#             end_times=station_setup.uplink.end,
#             ramp_rates=np.zeros(len(station_setup.uplink.ref_freq)).tolist(),
#             start_frequency=station_setup.uplink.ref_freq,
#         )

#         # Assign uplink interpolator to station
#         station = bodies.get("Earth").get_ground_station(station_id)
#         station.set_transmitting_frequency_calculator(uplink_interpolator)

#     # Environment updates for ionospheric correction
#     light_corrections = config.estimation.light_propagation.corrections
#     if light_corrections.ionospheric.present:

#         match light_corrections.ionospheric.model:

#             case "ionex":

#                 log.debug("IONEX ionospheric correction")

#                 # Add IONEX information to ground station objects
#                 ionex_files = [
#                     str(file) for file in light_corrections.ionospheric.sources
#                 ]
#                 tlight.set_ionosphere_model_from_ionex(
#                     data_files=ionex_files,
#                     bodies=bodies,
#                 )

#             case _:
#                 raise NotImplementedError(
#                     "Invalid ionospheric model"
#                     f": {light_corrections.ionospheric.model}"
#                 )

#     # Environment updates for tropospheric correction
#     if light_corrections.tropospheric.present:

#         match light_corrections.tropospheric.model:

#             case "vmf3":

#                 log.debug("VMF3 tropospheric correction")

#                 # Add VMF3 information to ground stations
#                 vmf3_files = [
#                     str(file) for file in light_corrections.tropospheric.sources
#                 ]
#                 tlight.set_vmf_troposphere_data(
#                     data_files=vmf3_files,
#                     bodies=bodies,
#                     file_has_meteo=True,
#                     file_has_gradient=True,
#                     set_meteo_data=True,
#                     set_troposphere_data=True,
#                 )

#             case _:
#                 raise NotImplementedError(
#                     "Invalid tropospheric model: "
#                     + light_corrections.tropospheric.model
#                 )

#     # Environment updates for closed-loop observations
#     if config.estimation.observation_models.closed_loop.present:

#         for (
#             vehicle,
#             vehicle_setup,
#         ) in config.environment.vehicles.items():

#             log.info(f"Closed-loop estimation settings for {vehicle}")

#             if vehicle != "MEX":
#                 raise NotImplementedError("Only MEX supported")

#             # Set turnaround ratio
#             match vehicle_setup.systems.turnaround_ratio:

#                 case "default":
#                     log.debug("Default turnaround ratios")
#                     bodies.get(
#                         vehicle
#                     ).system_models.set_default_transponder_turnaround_ratio_function()
#                 case _:
#                     raise NotImplementedError("Invalid turnaround ratio setup")

#             # Set reference point for tracking
#             match vehicle_setup.systems.reference_point:

#                 case "origin":
#                     log.debug("Reference point is origin")
#                     pass
#                 case "HGA":
#                     log.debug("Reference point is HGA")
#                     # Get state of HGA wrt LVI in fixed frame
#                     cstate_hga_lvi_fixed = (
#                         spice.get_body_cartesian_state_at_epoch(
#                             target_body_name="MEX_HGA",
#                             observer_body_name="MEX_SPACECRAFT",
#                             reference_frame_name="MEX_SPACECRAFT",
#                             aberration_corrections="none",
#                             ephemeris_time=config.time.initial_epoch,
#                         )
#                     )

#                     # Define constant ephemerides
#                     hga_ephemeris_settings = tephs.constant(
#                         constant_state=cstate_hga_lvi_fixed,
#                         frame_origin="MEX_SPACECRAFT",
#                         frame_orientation="MEX_SPACECRAFT",
#                     )
#                     hga_ephemerides = tephs.create_ephemeris(
#                         hga_ephemeris_settings, "HGA"
#                     )

#                     # Update MEX with reference point
#                     bodies.get(vehicle).system_models.set_reference_point(
#                         reference_point="HGA",
#                         ephemeris=hga_ephemerides,
#                     )

#                 case _:
#                     raise NotImplementedError(
#                         "Invalid reference point settings"
#                     )

#     return bodies
