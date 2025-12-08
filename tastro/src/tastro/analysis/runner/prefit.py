from ..core import AnalysisManagerBase
from ...io.command_line.core import CommandLineInput
from ...io import PrefitResults
from ...config import CaseSetup
from ...environment import system_of_bodies_from_config
from ...environment.vehicle import VehicleSettings
from ...estimation import EstimationManager
from tudatpy.estimation.observations_setup import (
    observations_simulation_settings as tosim,
)
from tudatpy.estimation import observations as tobs
from tudatpy.dynamics.environment_setup import ephemeris as teph
from ...logging import log
from tudatpy.interface import spice
from pathlib import Path
import multiprocessing as mp


class PrefitManager(AnalysisManagerBase[CommandLineInput]):

    def __init__(self, user_input: CommandLineInput) -> None:

        super().__init__(user_input)

        self.configs: dict[Path, "CaseSetup"] = {
            source: CaseSetup.from_config_file(source / "configuration.yaml")
            for source in self.user_input.source_dirs
        }

        # Update configurations to perform estimation
        for source in self.configs:
            self.configs[source].perform_estimation = True

        # # Load configuration from command line input
        # self.config = CaseSetup.from_config_file(
        #     self.source_dir / "configuration.yaml"
        # )
        # self.config.perform_estimation = True

        return None

    @staticmethod
    def calculate_prefit_residuals_single(
        source: "Path", config: "CaseSetup"
    ) -> None:

        log.info(f"Running pre-fit analysis for {source}")

        try:
            # Load metakernel
            spice.load_kernel(str(source / "metak.tm"))

            # Ensure vehicle has ephemerides
            vehicle_name = config.environment.general.spacecraft
            vehicle_config = config.environment.vehicles[vehicle_name]
            if not vehicle_config.ephemerides.present:
                log.warning(f"Forcing {vehicle_name} to have ephemerides")
                vehicle_config.ephemerides.present = True

            # Create system of bodies and estimation manager
            bodies = system_of_bodies_from_config(config)
            estimation_manager = EstimationManager(config)

            # Load observations and create observation models
            observations = estimation_manager.observation_collection(bodies)
            observation_models = estimation_manager.observation_models(
                observations
            )

            # Create observation simulator
            observation_simulator = tosim.create_observation_simulators(
                observation_settings=observation_models,
                bodies=bodies,
            )

            # Calculate pre-fit residuals
            log.info(f"Calculating pre-fit residuals")
            tobs.compute_residuals_and_dependent_variables(
                observation_collection=observations,
                observation_simulators=observation_simulator,
                bodies=bodies,
            )

            # Save output
            log.info(f"Saving results of pre-fit analysis: {source}")
            output = PrefitResults.from_observation_collection(observations)
            output.save_to_file(source / "prefit_results.pkl")

        finally:
            spice.clear_kernels()

        return None

    def calculate_prefit_residuals(self) -> None:

        with mp.Pool(4) as pool:

            pool.starmap(
                self.calculate_prefit_residuals_single,
                self.configs.items(),
            )

        return None


# def calculate_prefit_residuals(user_input: CommandLineInput) -> None:

#     # Initialize manager
#     manager = PrefitManager(user_input)

#     try:
#         # Load metakernel
#         spice.load_kernel(str(manager.source_dir / "metak.tm"))

#         # Ensure vehicle has ephemerides
#         vehicle_name = manager.config.environment.general.spacecraft
#         vehicle_config = manager.config.environment.vehicles[vehicle_name]
#         if not vehicle_config.ephemerides.present:
#             log.warning(f"Forcing {vehicle_name} to have ephemerides")
#             vehicle_config.ephemerides.present = True

#         # Create system of bodies from user input
#         bodies = system_of_bodies_from_config(manager.config)

#         # # Ensure vehicle has ephemerides
#         # vehicle_name = manager.config.environment.general.spacecraft
#         # vehicle_config = manager.config.environment.vehicles[vehicle_name]
#         # if not vehicle_config.ephemerides.present:

#         #     log.info(f"Manually setting ephemerides for {vehicle_name}")

#         #     # Initialize generator for vehicle settings
#         #     generator = VehicleSettings(
#         #         vehicle_name, vehicle_config, manager.config
#         #     )

#         #     # Create ephemerides from settings
#         #     bodies.get(vehicle_name).ephemeris = teph.create_ephemeris(
#         #         ephemeris_settings=generator.ephemeris_settings(),
#         #         body_name=vehicle_name,
#         #     )

#         # Create estimation manager
#         estimation_manager = EstimationManager(manager.config)

#         # Load observations and create observation models
#         observations = estimation_manager.observation_collection(bodies)
#         observation_models = estimation_manager.observation_models(observations)

#         # Create observation simulator
#         observation_simulators = tosim.create_observation_simulators(
#             observation_settings=observation_models,
#             bodies=bodies,
#         )

#         # Calculate pre-fit residuals
#         log.info("Calculating pre-fit residuals")
#         tobs.compute_residuals_and_dependent_variables(
#             observation_collection=observations,
#             observation_simulators=observation_simulators,
#             bodies=bodies,
#         )

#         # Save output
#         log.info("Saving results of pre-fit analysis")
#         output = PrefitResults.from_observation_collection(observations)
#         output.save_to_file(manager.source_dir / "prefit_results.pkl")

#     finally:
#         spice.clear_kernels()

#     return None
