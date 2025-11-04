from ...io.command_line.core import CLInputFigure, CLParserFigure
from ..core import AnalysisFigureManagerBase
from ...io import EstimationResults, PropagationOutput
from nastro import types as nt, graphics as ng
import numpy as np


class CartesianResidualHistoryManager(AnalysisFigureManagerBase[CLInputFigure]):

    def __init__(self, user_input: CLInputFigure, residual_size: int) -> None:

        super().__init__(user_input)

        # Save residual size
        self.residual_size = residual_size

        # Load estimation results
        estimation = EstimationResults.from_file(
            self.source_dir / "estimation.pkl"
        )

        # Reshape residuals
        self.residual_history = self.reshape_residual_history(estimation)

        # Load propagation restuls
        propagation = PropagationOutput.from_file(
            self.source_dir / "results.pkl"
        )

        # Define propagation epochs and distance to Mars
        self.epochs = estimation.epochs
        self.propagation_epochs = propagation.epochs
        self.dmars = nt.CartesianState(*propagation.rstate_j2000.T).r_mag

        return None

    def reshape_residual_history(
        self, estimation: EstimationResults
    ) -> list[np.ndarray]:

        n_epochs = int(len(estimation.final_residuals) // self.residual_size)
        residual_history = [
            np.array(residual_set).reshape(n_epochs, self.residual_size).T
            for residual_set in estimation.residual_history.T
        ]
        return residual_history


def show_cartesian_residuals(user_input: CLInputFigure) -> None:

    # Define manager
    manager = CartesianResidualHistoryManager(user_input, 3)

    # Get source directory of observations from configuration
    observations_config = (
        manager.config.estimation.observations.cartesian.sources[0]
    )
    ref_id = manager.get_directory_id(observations_config.path.parent)
    if observations_config.use_ephemerides:
        ref_id += " [Eph]"

    # Calculate quanitities to plot
    dt = (manager.epochs - manager.ref_epoch) / 3600.0
    dt_prop = (manager.propagation_epochs - manager.ref_epoch) / 3600.0
    dmars = manager.dmars * 1e-7

    # Define setup for canvas
    source_id = manager.get_directory_id(manager.source_dir)
    canvas_setup = ng.PlotSetup(
        canvas_size=manager.default_canvas_size,
        canvas_title=f"Cartesian residuals :: {ref_id} :: {source_id}",
        show=manager.user_input.show,
        save=manager.user_input.save,
        dir=manager.user_input.source_dir,
        name="cartesian-residual-history.png",
    )

    # Define common setup for subfigures
    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {manager.ref_epoch_isot}",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
        scilimits=(-2, 3),
    )

    # Define figure
    with ng.Mosaic("ab;cd", canvas_setup) as canvas:

        for idx, label in enumerate(("x", "y", "z")):

            # Define setup for current subfigure
            current_setup = subfig_setup.version(
                ylabel=rf"$\Delta {label}$ [m]"
            )

            # Create subfigure
            with canvas.subplot(
                current_setup, generator=ng.DoubleAxis
            ) as subfig:

                subfig.line(
                    dt_prop, dmars, color="black", alpha=0.2, axis="right"
                )

                for jdx, residual_set in enumerate(manager.residual_history):

                    subfig.line(dt, residual_set[idx])

        # Add subfigure with residual magnitude
        rsetup = subfig_setup.version(ylabel=r"$||\Delta \mathbf{r}||$ [m]")
        with canvas.subplot(rsetup, generator=ng.DoubleAxis) as subfig:

            subfig.line(dt_prop, dmars, color="black", alpha=0.2, axis="right")
            for jdx, residual_set in enumerate(manager.residual_history):

                subfig.line(
                    dt,
                    np.linalg.norm(residual_set, axis=0),
                    label=f"Iteration {jdx}",
                )

    return None
