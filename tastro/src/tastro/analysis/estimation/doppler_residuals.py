from ...io.command_line.core import CLInputFigure, CLParserFigure
from ..core import AnalysisFigureManagerBase
from ...io import EstimationResults, PropagationOutput
from nastro import types as nt, graphics as ng
import numpy as np


class ResidualHistoryManager(AnalysisFigureManagerBase[CLInputFigure]):

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

        # Get index of best iteration
        self.best_iteration = estimation.best_iteration_index

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

    def get_1d_residual_set(self, residual_set: np.ndarray) -> np.ndarray:

        match self.residual_size:

            case 1:
                return residual_set[0]
            case 3:
                return np.linalg.norm(residual_set, axis=0)
            case _:
                raise NotImplementedError(
                    f"Invalid residual shape: {self.residual_size}"
                )

    def show_residual_history(self) -> None:

        # Calculate quantities to plot
        dt = (self.epochs - self.ref_epoch) / 3600.0
        dt_prop = (self.propagation_epochs - self.ref_epoch) / 3600.0
        dmars = self.dmars * 1e-7

        # Define setup for canvas
        source_id = self.get_directory_id(self.source_dir)
        canvas_setup = self.generate_canvas_setup(
            "closed-loop-residual-history.png"
        )
        canvas_setup.canvas_title = f"Estimation results :: {source_id}"

        # Generate mosaic based on size
        match self.residual_size:
            case 1:
                mosaic = "a;b"
                labels = ["Residual [mHz]"]
                scale = 1e3
            case 3:
                mosaic = "ab;cd"
                labels = ["dx [m]", "dy [m]", "dz [m]"]
                scale = 1
            case _:
                raise NotImplementedError(
                    f"Invalid residual size: {self.residual_size}"
                )

        # Generate common setup for subfigures
        with ng.Mosaic(mosaic, canvas_setup) as canvas:

            for idx in range(self.residual_size):

                # Generate setup for subfigure
                subfig_setup = self.generate_figure_setup(labels[idx])
                with canvas.subplot(
                    subfig_setup, generator=ng.DoubleAxis
                ) as subfig:

                    subfig.line(
                        dt_prop, dmars, color="black", alpha=0.2, axis="right"
                    )
                    for residual_set in self.residual_history:

                        subfig.line(
                            dt, residual_set[idx] * scale, fmt=".", markersize=2
                        )

            # Generate setup for final subfigure
            subfig_setup = self.generate_figure_setup("Final residual [mHz]")
            subfig_setup.axtitle = "Results for best iteration"
            with canvas.subplot(
                subfig_setup, generator=ng.DoubleAxis
            ) as subfig:

                subfig.line(
                    dt_prop, dmars, color="black", alpha=0.2, axis="right"
                )
                best = self.get_1d_residual_set(
                    self.residual_history[self.best_iteration]
                )

                subfig.line(dt, best * scale, fmt=".", markersize=2)


def show_doppler_residuals(user_input: CLInputFigure) -> None:

    # Define manager
    manager = ResidualHistoryManager(user_input, 1)

    # Calculate quanitities to plot
    dt = (manager.epochs - manager.ref_epoch) / 3600.0
    dt_prop = (manager.propagation_epochs - manager.ref_epoch) / 3600.0
    dmars = manager.dmars * 1e-7

    # Define setup for canvas
    source_id = manager.get_directory_id(manager.source_dir)
    canvas_setup = ng.PlotSetup(
        canvas_size=manager.default_canvas_size,
        canvas_title=f"Closed-loop Doppler residuals :: {source_id}",
        show=manager.user_input.show,
        save=manager.user_input.save,
        dir=manager.user_input.source_dir,
        name="closed-loop-residual-history.png",
        xlabel=f"Hours past {manager.ref_epoch_isot}",
        ylabel="Residual [mHz]",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
    )

    with ng.DoubleAxis(canvas_setup) as fig:

        fig.line(dt_prop, dmars, axis="right", color="black", alpha=0.2)
        for idx, residual_set in enumerate(manager.residual_history):
            fig.line(
                dt,
                residual_set[0] * 1e3,
                fmt=".",
                label=f"Iteration {idx + 1}",
            )

    # # Define common setup for subfigures
    # subfig_setup = ng.PlotSetup(
    #     xlabel=f"Hours past {manager.ref_epoch_isot}",
    #     rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
    #     scilimits=(-2, 3),
    # )

    # # Define figure
    # with ng.Mosaic("ab;cd", canvas_setup) as canvas:

    #     for idx, label in enumerate(("x", "y", "z")):

    #         # Define setup for current subfigure
    #         current_setup = subfig_setup.version(
    #             ylabel=rf"$\Delta {label}$ [m]"
    #         )

    #         # Create subfigure
    #         with canvas.subplot(
    #             current_setup, generator=ng.DoubleAxis
    #         ) as subfig:

    #             subfig.line(
    #                 dt_prop, dmars, color="black", alpha=0.2, axis="right"
    #             )

    #             for jdx, residual_set in enumerate(manager.residual_history):

    #                 subfig.line(dt, residual_set[idx])

    #     # Add subfigure with residual magnitude
    #     rsetup = subfig_setup.version(ylabel=r"$||\Delta \mathbf{r}||$ [m]")
    #     with canvas.subplot(rsetup, generator=ng.DoubleAxis) as subfig:

    #         subfig.line(dt_prop, dmars, color="black", alpha=0.2, axis="right")
    #         for jdx, residual_set in enumerate(manager.residual_history):

    #             subfig.line(
    #                 dt,
    #                 np.linalg.norm(residual_set, axis=0),
    #                 label=f"Iteration {jdx}",
    #             )

    return None
