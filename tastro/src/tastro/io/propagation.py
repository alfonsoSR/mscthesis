from ..config import CaseSetup
from tudatpy.astro import (
    time_representation as ttime,
    frame_conversion as tframe,
)
from tudatpy.dynamics import propagation as tprop
from tudatpy.interface import spice
import numpy as np
from dataclasses import dataclass
import pickle
from pathlib import Path


@dataclass
class PropagationOutput:

    epochs: np.ndarray
    cstate_j2000: np.ndarray
    rstate_j2000: np.ndarray
    cstate_rsw: np.ndarray
    rstate_rsw: np.ndarray
    number_of_function_evaluations: float
    ustate: np.ndarray

    def save_to_file(self, file: Path) -> None:

        with file.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @property
    def state_history(self) -> dict[ttime.Time, np.ndarray]:

        return {
            epoch: cstate for (epoch, cstate) in zip(self.epochs, self.cstate_j2000)
        }

    @staticmethod
    def from_file(file: Path) -> "PropagationOutput":

        with file.open("rb") as buffer:
            _self = pickle.load(buffer)

        return _self

    @staticmethod
    def from_config_file(config_file: Path) -> "PropagationOutput":

        # Get path to results from path to config file
        results_file = config_file.parent / "results.pkl"
        return PropagationOutput.from_file(results_file)

    @staticmethod
    def from_simulation(
        simulation: tprop.SingleArcSimulationResults, config: CaseSetup
    ) -> "PropagationOutput":

        # Get ephemerides of spacecraft at propagation epochs
        reference_state_j2000 = np.array(
            [
                spice.get_body_cartesian_state_at_epoch(
                    target_body_name=config.environment.general.spacecraft,
                    observer_body_name=config.environment.general.center,
                    reference_frame_name=config.environment.general.global_frame_orientation,
                    aberration_corrections="NONE",
                    ephemeris_time=epoch,
                )
                for epoch in simulation.state_history_time_object.keys()
            ]
        )

        # Propagation epochs as float
        propagation_epochs = np.array(list(simulation.state_history.keys()))

        # Propagated state vector (Cartesian global frame)
        propagated_state_j2000 = np.array(list(simulation.state_history.values()))

        # Get states in RSW frame (ephemerides)
        reference_state_rsw = np.zeros_like(reference_state_j2000)
        propagated_state_rsw = np.zeros_like(propagated_state_j2000)
        for idx, rstate in enumerate(reference_state_j2000):

            # Calculate rotation matrix
            rotation_matrix = tframe.inertial_to_rsw_rotation_matrix(rstate)

            # Reference state in RSW
            _refpos_rsw = rotation_matrix @ rstate[:3]
            _refvel_rsw = rotation_matrix @ rstate[3:]
            reference_state_rsw[idx] = np.array([_refpos_rsw, _refvel_rsw]).flatten()

            # Propagated state in RSW
            _ppos_rsw = rotation_matrix @ propagated_state_j2000[idx][:3]
            _pvel_rsw = rotation_matrix @ propagated_state_j2000[idx][3:]
            propagated_state_rsw[idx] = np.array([_ppos_rsw, _pvel_rsw]).flatten()

        return PropagationOutput(
            epochs=propagation_epochs,
            cstate_j2000=propagated_state_j2000,
            rstate_j2000=reference_state_j2000,
            cstate_rsw=propagated_state_rsw,
            rstate_rsw=reference_state_rsw,
            number_of_function_evaluations=simulation.total_number_of_function_evaluations,
            ustate=np.array(list(simulation.unprocessed_state_history.values())).T,
        )
