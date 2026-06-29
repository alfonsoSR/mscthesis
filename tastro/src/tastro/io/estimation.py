from dataclasses import dataclass
import numpy as np
from .propagation import PropagationOutput
from tudatpy.dynamics import propagation as tprop
from tudatpy.estimation import estimation_analysis as testa
from tudatpy.estimation.observations import ObservationCollection
from ..config import CaseSetup
from pathlib import Path
import pickle
from ..logging import log
import traceback

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    from tudatpy.estimation import observations as tobs


class EstimationResults:

    def __init__(
        self,
        best_iteration_index: int,
        correlation_matrix: np.ndarray,
        covariance_matrix: np.ndarray,
        parameter_history: np.ndarray,
        residual_history: np.ndarray,
        estimation_epochs: np.ndarray,
        simulation_history: list[PropagationOutput] | None = None,
        covariance_history: dict[float, np.ndarray] | None = None,
    ) -> None:

        self.best_iteration_index = best_iteration_index
        self.correlation_matrix = correlation_matrix
        self.covariance_matrix = covariance_matrix
        self.parameter_history = parameter_history
        self.residual_history = residual_history
        self.simulation_history = simulation_history
        self.epochs = estimation_epochs
        self.covariance_history = covariance_history

        return None

    @staticmethod
    def from_estimator_output_and_observations(
        estimator: "testa.Estimator",
        estimation_output: "testa.EstimationOutput",
        observations: "ObservationCollection",
        config: "CaseSetup",
    ) -> "EstimationResults":

        # Extract epochs from observations
        estimation_epochs = np.array(
            observations.get_concatenated_observation_times()
        )

        # Load simulation results per iteration into data structures
        simulation_results: "list[PropagationOutput]" = []
        for idx, result in enumerate(
            estimation_output.simulation_results_per_iteration
        ):

            # Sanity check for type information
            if not isinstance(
                result, tprop.SingleArcVariationalSimulationResults
            ):
                continue

            # Extract translational results and verify type
            simulation = result.dynamics_results
            if not isinstance(simulation, tprop.SingleArcSimulationResults):
                continue

            # Update list
            simulation_results.append(
                PropagationOutput.from_simulation(simulation, config)
            )

        # Propagate covariance
        covariance_epochs = simulation_results[-1].epochs.tolist()
        print(f"Number of epochs: {len(covariance_epochs)}")
        covariance_history = testa.propagate_covariance(
            initial_covariance=estimation_output.covariance,
            state_transition_interface=estimator.state_transition_interface,
            output_times=covariance_epochs,
        )

        # Return data structure
        return EstimationResults(
            best_iteration_index=estimation_output.best_iteration,
            correlation_matrix=estimation_output.correlations,
            covariance_matrix=estimation_output.covariance,
            parameter_history=estimation_output.parameter_history,
            residual_history=estimation_output.residual_history,
            estimation_epochs=estimation_epochs,
            simulation_history=simulation_results,
            covariance_history=covariance_history,
        )

    @classmethod
    def from_estimation_output(
        cls,
        results: testa.EstimationOutput,
        observations: ObservationCollection,
    ) -> "EstimationResults":

        # Extract epochs from observations
        epochs = np.array(observations.get_concatenated_observation_times())

        return EstimationResults(
            best_iteration_index=results.best_iteration,
            correlation_matrix=results.correlations,
            covariance_matrix=results.covariance,
            parameter_history=results.parameter_history,
            residual_history=results.residual_history,
            estimation_epochs=epochs,
            simulation_history=None,
        )

    @property
    def final_residuals(self) -> np.ndarray:

        return self.residual_history[:, self.best_iteration_index]

    @property
    def final_parameters(self) -> np.ndarray:

        return self.parameter_history[:, self.best_iteration_index]

    @property
    def final_propagation_results(self) -> "PropagationOutput":

        if self.simulation_history is None:
            log.fatal("Simulation history not available")
            exit(1)

        return self.simulation_history[self.best_iteration_index]

    def save_to_file(self, file: Path) -> None:

        log.info(f"Saving estimation results in {file}")

        with file.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @classmethod
    def from_file(cls, file: Path) -> "EstimationResults":

        if not file.exists():
            log.error(
                f"Failed to load estimation results from {file} :"
                "File not found"
            )
            log.error(traceback.extract_stack()[-2])
            exit(1)

        with file.open("rb") as buffer:
            _self = pickle.load(buffer)

        if not isinstance(_self, cls):
            raise TypeError("Invalid input file for estimation results")

        return _self


@dataclass
class PrefitResults:

    epochs: np.ndarray
    simulated: np.ndarray
    measured: np.ndarray
    residual: np.ndarray

    @classmethod
    def from_observation_collection(
        cls, observations: "tobs.ObservationCollection"
    ) -> "PrefitResults":

        epochs = np.array(observations.get_concatenated_observation_times())
        simulated = np.array(
            observations.get_concatenated_computed_observations()
        )
        measured = np.array(observations.get_concatenated_observations())
        residuals = np.array(observations.get_concatenated_residuals())

        return PrefitResults(epochs, simulated, measured, residuals)

    @classmethod
    def from_prefit_results(
        cls, results: "tobs.ObservationCollection"
    ) -> "PrefitResults":
        """DEPRECATED: Use from_observation_collection instead"""

        log.warning(traceback.extract_stack()[-2])
        log.warning(
            "PrefitResults.from_prefit_results is deprecated. "
            "Use from_observation_collection instead"
        )

        return cls.from_observation_collection(results)

    def save_to_file(self, file_name: Path) -> None:

        with file_name.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @classmethod
    def from_file(cls, file_name: Path) -> "PrefitResults":

        if not file_name.exists():
            log.fatal(f"Failed to load pre-fit results: {file_name} not found")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        with file_name.open("rb") as buffer:
            output = pickle.load(buffer)

        assert isinstance(output, cls)

        return output
