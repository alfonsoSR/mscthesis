from dataclasses import dataclass
import numpy as np
from .propagation import PropagationOutput
from tudatpy.dynamics import propagation as tprop
from tudatpy.estimation import estimation_analysis as testa
from ..config import CaseSetup
from pathlib import Path
import pickle
from ..logging import log


class EstimationResults:

    def __init__(
        self,
        best_iteration_index: int,
        correlation_matrix: np.ndarray,
        covariance_matrix: np.ndarray,
        parameter_history: np.ndarray,
        residual_history: np.ndarray,
        simulation_history: list[PropagationOutput],
    ) -> None:

        self.best_iteration_index = best_iteration_index
        self.correlation_matrix = correlation_matrix
        self.covariance_matrix = covariance_matrix
        self.parameter_history = parameter_history
        self.residual_history = residual_history
        self.simulation_history = simulation_history

        return None

    @classmethod
    def from_estimation_output(
        cls, results: testa.EstimationOutput, config: CaseSetup
    ) -> "EstimationResults":

        best_iteration_index = results.best_iteration
        correlation_matrix = results.correlations
        covariance_matrix = results.covariance
        parameter_history = results.parameter_history
        residual_history = results.residual_history

        # Generate propagation output objects from simulation history
        simulation_history: list[PropagationOutput] = []
        for simulation in results.simulation_results_per_iteration:
            assert isinstance(simulation, tprop.SingleArcSimulationResults)
            simulation_history.append(
                PropagationOutput.from_simulation(simulation, config)
            )

        return EstimationResults(
            best_iteration_index=best_iteration_index,
            correlation_matrix=correlation_matrix,
            covariance_matrix=covariance_matrix,
            parameter_history=parameter_history,
            residual_history=residual_history,
            simulation_history=simulation_history,
        )

    @property
    def final_residuals(self) -> np.ndarray:

        return self.residual_history[:, self.best_iteration_index]

    @property
    def final_parameters(self) -> np.ndarray:

        return self.parameter_history[:, self.best_iteration_index]

    def save_to_file(self, file: Path) -> None:

        log.info(f"Saving estimation results in {file}")

        with file.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @classmethod
    def from_file(cls, file: Path) -> "EstimationResults":

        with file.open("rb") as buffer:
            _self = pickle.load(buffer)

        if not isinstance(_self, cls):
            raise TypeError("Invalid input file for estimation results")

        return _self
