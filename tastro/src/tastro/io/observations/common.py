import numpy as np
from tudatpy.astro import time_representation as ttime
from ...logging import log
import traceback
import abc

from typing import TYPE_CHECKING, Self

if TYPE_CHECKING:
    from ...config import CaseSetup


class RawObservationRecord(metaclass=abc.ABCMeta):

    def __init__(
        self,
        epochs: np.ndarray,
        observations: np.ndarray,
        noise: np.ndarray | None = None,
    ) -> None:

        # Ensure that epochs are Time objects
        if not isinstance(epochs[0], ttime.Time):
            log.fatal("Observation epochs must be Time objects")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Sort observations and epochs chronologicaly
        idx_order = np.argsort([epoch.to_float() for epoch in epochs])
        self.epochs: np.ndarray = epochs[idx_order]
        self.observations: np.ndarray = observations[idx_order]

        # Load and sort noise
        if noise is None:
            self.noise: np.ndarray = np.zeros_like(self.observations)
        else:
            self.noise: np.ndarray = noise[idx_order]

        return None
