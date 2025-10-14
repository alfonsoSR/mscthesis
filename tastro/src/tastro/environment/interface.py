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
from ..estimation import update_system_of_bodies


def system_of_bodies_from_config(config: "CaseSetup") -> "tenv.SystemOfBodies":

    # Initialize empty environment settings
    environment_settings = tenvs.BodyListSettings(
        frame_origin=config.environment.general.global_frame_origin,
        frame_orientation=config.environment.general.global_frame_orientation,
    )

    # Define settings for vehicles
    for spacecraft, spacecraft_setup in config.environment.vehicles.items():

        if not spacecraft_setup.present:
            continue

        log.info(f"Generating body settings for vehicle: {spacecraft}")

        # Define empty settings for vehicle
        environment_settings.add_empty_settings(spacecraft)
        spacecraft_settings = environment_settings.get(spacecraft)

        # Initialize object to generate settings from configuration
        generator = VehicleSettings(spacecraft, spacecraft_setup, config)

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

        # Radiation target settings
        if spacecraft_setup.radiation.present:
            spacecraft_settings.radiation_pressure_target_settings = (
                generator.radiation_target_settings()
            )

    # Define settings for planets
    for planet, planet_setup in config.environment.planets.items():

        if not planet_setup.present:
            continue

        log.info(f"Generating body settings for planet: {planet}")

        # Define empty settings for planet
        environment_settings.add_empty_settings(planet)
        planet_settings = environment_settings.get(planet)

        # Initialize object to generate settings from configuration
        generator = PlanetSettings(planet, planet_setup, config)

        # Ephemeris settings
        if planet_setup.ephemerides.present:
            planet_settings.ephemeris_settings = generator.ephemeris_settings()

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

    # Define settings for ground stations
    if "Earth" in config.environment.planets and config.estimation.present:

        # Get settings for Earth
        ground_station_settings = []

        # Define settings for ground stations
        for station, station_setup in config.environment.stations.items():

            if not station_setup.present:
                continue

            log.info(f"Generating body settings for ground station: {station}")

            # Initialize object to generate settings from configuration
            generator = StationSettings(station, station_setup, config)

            # Define settings
            ground_station_settings.append(generator.station_settings())

        # Add settings to object
        environment_settings.get("Earth").ground_station_settings = (
            ground_station_settings
        )

    # Create system of bodies
    log.info("Creating system of bodies")
    bodies = tenvs.create_system_of_bodies(environment_settings)

    # If estimation is present, update system of bodies
    if config.estimation.present:
        bodies = update_system_of_bodies(config, bodies)

    return bodies
