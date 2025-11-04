raise DeprecationWarning("This submodule is deprecated")

from .. import (
    environment as nenv,
    io as nio,
    config as ncon,
    estimation as nest,
)
import argparse
from dataclasses import dataclass
from pathlib import Path
from tudatpy.estimation.observations_setup import (
    observations_simulation_settings as tosim,
)
from tudatpy.dynamics.environment_setup import ephemeris as teph
from tudatpy.estimation import observations as tobs
from tudatpy.estimation.observations import observations_processing as tobsp
import numpy as np
import pickle
from tudatpy.interface import spice
from tudatpy.astro import time_representation as ttime
from nastro import graphics as ng, types as nt
from .utils import get_propagation_start_epoch_from_config


@dataclass
class PrefitResults:

    epochs: np.ndarray
    simulated: np.ndarray
    measured: np.ndarray
    residual: np.ndarray

    @classmethod
    def from_prefit_results(
        cls, results: "tobs.ObservationCollection"
    ) -> "PrefitResults":

        epochs = np.array(results.get_concatenated_observation_times())
        simulated = np.array(results.get_concatenated_computed_observations())
        measured = np.array(results.get_concatenated_observations())
        residuals = np.array(results.get_concatenated_residuals())

        return PrefitResults(epochs, simulated, measured, residuals)

    def save_to_file(self, file_name: Path) -> None:

        with file_name.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @classmethod
    def from_file(cls, file_name: Path) -> "PrefitResults":

        with file_name.open("rb") as buffer:
            output = pickle.load(buffer)

        assert isinstance(output, cls)

        return output


@dataclass
class PrefitInput:

    source_dir: Path


class PrefitParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__(prog="tprefit")

        self.add_argument("source_dir", help="Source directory")

        return None

    def parse_args(self) -> "PrefitInput":

        default = super().parse_args()

        # Define path to source directory
        source_dir = Path(default.source_dir).resolve()
        if not source_dir.exists():
            raise FileNotFoundError(f"Invalid source directory: {source_dir}")

        return PrefitInput(source_dir=source_dir)


def prefit_closed_loop(setup: "PrefitInput") -> None:

    # Load configuration
    config = ncon.CaseSetup.from_config_file(
        setup.source_dir / "configuration.yaml"
    )
    config.perform_estimation = True
    config.perform_propagation = True

    # Load spice kernels
    metakernel = str(setup.source_dir / "metak.tm")
    try:
        spice.load_kernel(metakernel)

        # Create system of bodies
        bodies = nenv.system_of_bodies_from_config(config)

        # Load observations
        manager = nest.EstimationManager(config)
        observations = manager.observation_collection(bodies)
        observation_models = manager.observation_models(observations)

        # Create observation simulator
        simulator = tosim.create_observation_simulators(
            observation_settings=observation_models,
            bodies=bodies,
        )

        # Compute pre-fit residuals
        tobs.compute_residuals_and_dependent_variables(
            observation_collection=observations,
            observation_simulators=simulator,
            bodies=bodies,
        )

        with (setup.source_dir / "epochs.pkl").open("wb") as buffer:
            pickle.dump(observations.concatenated_times, buffer)

        # Save output
        PrefitResults.from_prefit_results(observations).save_to_file(
            setup.source_dir / "prefit_results.pkl"
        )

    finally:
        spice.clear_kernels()

    return None


def prefit_closed_loop_cli() -> None:

    # Read user input from command line
    user_input = PrefitParser().parse_args()

    # Calculate pre-fit residuals and save output
    prefit_closed_loop(user_input)

    return None


def plot_closed_loop_prefits(setup: "PrefitInput") -> None:

    # Load results
    results = PrefitResults.from_file(setup.source_dir / "prefit_results.pkl")

    # Get time of closest approach from configuration
    config = ncon.CaseSetup.from_config_file(
        setup.source_dir / "configuration.yaml"
    )
    t0 = get_propagation_start_epoch_from_config(config)
    t0_isot = ttime.DateTime.from_epoch_time_object(t0).to_iso_string(
        add_T=True, number_of_digits_seconds=0
    )

    # Load propagation results
    propagation = nio.PropagationOutput.from_file(
        setup.source_dir / "results.pkl"
    )

    # Turn epochs into hours past closest approach
    dt = (results.epochs - t0.to_float()) / 3600.0
    dt_propagation = (propagation.epochs - t0.to_float()) / 3600.0

    # Cartesian state with respect to Mars
    cstate = nt.CartesianState[nt.Vector](*propagation.rstate_j2000.T)

    with ng.Mosaic("a;b") as canvas:

        # Plot simulated observations
        with canvas.subplot(generator=ng.DoubleAxis) as subfig:

            subfig.line(dt, results.simulated, fmt=".", label="Simulated")
            subfig.line(dt, results.measured, fmt=".", label="Observed")

            subfig.line(dt_propagation, cstate.r_mag, alpha=0.4, axis="right")

        # Plot pre-fit residuals
        with canvas.subplot() as subfig:

            subfig.line(dt, results.residual, fmt=".", markersize=2)

    return None


def plot_closed_loop_prefits_cli() -> None:

    # Read user input from command line
    user_input = PrefitParser().parse_args()

    # Generate figure
    plot_closed_loop_prefits(user_input)

    return None
