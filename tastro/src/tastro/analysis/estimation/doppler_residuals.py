from ...io.command_line.core import (
    CLInputFigure,
    CLParserFigure,
    CommandLineInput,
)
from ..core import AnalysisFigureManagerBase
from ...io import EstimationResults, PropagationOutput
from nastro import types as nt, graphics as ng
import numpy as np
from ...logging import log
from pathlib import Path
from matplotlib import pyplot as plt
from ...config import CaseSetup


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

        # Define a mask for the propagation epochs based on observations
        self.epochs = estimation.epochs
        observation_period_mask = (
            propagation.epochs >= np.min(self.epochs)
        ) * (propagation.epochs <= np.max(self.epochs))

        # Get masked propagation epochs and distance to Mars
        self.propagation_epochs = propagation.epochs[observation_period_mask]
        self.dmars = nt.CartesianState(*propagation.rstate_j2000.T).r_mag[
            observation_period_mask
        ]

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
                    for jdx, residual_set in enumerate(self.residual_history):

                        subfig.line(
                            dt,
                            residual_set[idx] * scale,
                            fmt=".",
                            markersize=2,
                            label=f"Iteration {jdx + 1}",
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


def show_estimation_summary_single(source_dir: Path) -> None:

    # Define path to output file and show info
    output_file = source_dir / "estimation-summary.log"
    log.info(f"Writing estimation summary to {output_file}")

    # Load estimation output
    estimation = EstimationResults.from_file(source_dir / "estimation.pkl")

    with output_file.open("w") as buffer:

        # Loop over parameter and residual histories
        for step, residuals in enumerate(estimation.residual_history.T):

            # Get parameter update in current step
            previous_parameters = estimation.parameter_history[:, step]
            current_parameters = estimation.parameter_history[:, step + 1]
            parameter_update = current_parameters - previous_parameters
            update_string = " ".join(
                [f"{item:.5e}" for item in parameter_update]
            )

            # Get RMS of current residuals
            n_observations = len(residuals)
            residual_rms = np.linalg.norm(residuals) / np.sqrt(n_observations)

            # Display information
            buffer.write(
                f"Estimation step: {step + 1} - {n_observations} points\n"
            )
            buffer.write(f"Current residual: {residual_rms:.7f}\n")
            buffer.write(f"Parameter update: {update_string}\n")

        # Get final residuals
        final_residuals = estimation.residual_history[
            :, estimation.best_iteration_index
        ]
        final_rms = np.linalg.norm(final_residuals) / np.sqrt(
            len(final_residuals)
        )

        # Get final parameters
        final_parameters = estimation.parameter_history[
            :, estimation.best_iteration_index
        ]
        final_parameters_string = " ".join(
            [f"{item:.5e}" for item in final_parameters]
        )

        # Get difference between initial and final parameters
        initial_parameters = estimation.parameter_history[:, 0]
        final_parameter_update = final_parameters - initial_parameters
        final_parameter_update_string = " ".join(
            [f"{item:.5e}" for item in final_parameter_update]
        )

        # Show final residual and parameters
        buffer.write(f"Final residual: {final_rms:.7f}\n")
        buffer.write(
            f"Total parameter update: {final_parameter_update_string}\n"
        )
        buffer.write(f"Final parameters: {final_parameters_string}\n")


def show_correlation_matrix_single(source_dir: Path) -> None:

    # Load estimation results
    estimation = EstimationResults.from_file(source_dir / "estimation.pkl")

    # Configuration
    config = CaseSetup.from_config_file(source_dir / "configuration.yaml")
    param_config = config.estimation.parameters

    # Estimated parameters
    parameters = ["$x$", "$y$", "$z$", r"$\dot x$", r"$\dot y$", r"$\dot z$"]
    if param_config.drag_coefficient:
        parameters.append("$C_D$")
    if param_config.arcwise_drag_coefficient:
        parameters += ["$C_{D, 1}$", "$C_{D,2}$", "$C_{D,3}$", "$C_{D,4}$"]
    if param_config.radiation_pressure_coefficient:
        parameters.append("$K_1$")

    # assert estimation.covariance_history is not None
    # formal_errors = CartesianStateUncertainty.from_covariance_history(
    #     estimation.covariance_history
    # )
    # formal_error_epochs = np.array(list(estimation.covariance_history.keys()))

    # propagation = estimation.final_propagation_results
    # cstate = nt.CartesianState(*propagation.cstate_j2000.T)
    # rstate = nt.CartesianState(*propagation.rstate_j2000.T)
    # residual = cstate - rstate

    # with ng.SingleAxis() as fig:

    #     for term in ("x", "y", "z"):

    #         fig.line(propagation.epochs, getattr(residual, term), label=term)
    #         fig.line(
    #             formal_error_epochs, 0 * formal_error_epochs, color="w", alpha=0
    #         )
    #         fig.boundary(3.0 * getattr(formal_errors, term), alpha=0.8)

    # exit(0)

    correlation = estimation.covariance_matrix

    # print(covariance.shape)
    # exit(0)

    # # Get correlation matrix
    # correlation = estimation.correlation_matrix

    setup = ng.PlotSetup(
        aspect="equal",
        canvas_size=(6, 6),
        grid=False,
        minor_ticks=False,
        canvas_title="Correlation matrix",
    )

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("Correlation matrix")
    foo = ax.imshow(correlation, aspect="equal", cmap="GnBu")
    ax.set_xticks([idx for idx in range(len(parameters))], labels=parameters)
    ax.set_yticks([idx for idx in range(len(parameters))], labels=parameters)
    fig.colorbar(mappable=foo)
    for i in range(correlation.shape[0]):
        # for j in range(correlation.shape[1]):

        val = np.sqrt(correlation[i, i])
        # dcor = 1 - val

        ax.text(
            x=i,
            y=i,
            s=f"{val:.3e}",
            ha="center",
            va="center",
            color="k",
        )  # type: ignore

    fig.savefig(source_dir / "correlation_matrix.png")
    plt.show()

    # with ng.SingleAxis(setup) as fig:

    #     fig.imshow(np.abs(correlation))
    #     fig.axes["left"].set_xticks(
    #         [idx for idx in range(len(parameters))], labels=parameters
    #     )
    #     fig.axes["left"].set_yticks(
    #         [idx for idx in range(len(parameters))], labels=parameters
    #     )

    # for i in range(correlation.shape[0]):
    #     for j in range(correlation.shape[1]):

    #         fig.axes["left"].text(
    #             x=i,
    #             y=j,
    #             s=f"{correlation[i, j]:.2e}",
    #             ha="center",
    #             va="center",
    #             color="k",
    #         )  # type: ignore

    return None
