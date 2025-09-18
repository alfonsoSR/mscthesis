import numpy as np
from dataclasses import dataclass
from pathlib import Path
import pickle


@dataclass
class PropagationOutput:

    epochs: np.ndarray
    cstate_j2000: np.ndarray
    rstate_j2000: np.ndarray
    cstate_rsw: np.ndarray
    rstate_rsw: np.ndarray

    def save_to_file(self, file: Path) -> None:

        with file.open("wb") as buffer:
            pickle.dump(self, buffer)

        return None

    @staticmethod
    def from_file(file: Path) -> "PropagationOutput":

        with file.open("rb") as buffer:
            _self = pickle.load(buffer)

        return _self
