from ...io.command_line.core import CLInputFigure, CommandLineInput
from ..core import AnalysisFigureManagerBase
from ...io import EstimationResults, PrefitResults
from pathlib import Path
from ...logging import log
import traceback
from nastro import graphics as ng, types as nt
import numpy as np


class DopplerEstimationManager(AnalysisFigureManagerBase[CLInputFigure]):

    residual_scale: float = 1.0e3
    """Change of units for residuals"""

    def __init__(self, user_input: CLInputFigure) -> None:

        super().__init__(user_input)

        self.estimation_results: dict[Path, EstimationResults] = {}
        self.prefit_results: dict[Path, PrefitResults] = {}
        self.complete_sources: list[Path] = self.user_input.source_dirs

        for source in self.complete_sources:

            # Load estimation results if possible
            estimation_source = source / "estimation.pkl"
            if estimation_source.exists():
                self.estimation_results[source] = EstimationResults.from_file(
                    estimation_source
                )

            # Load pre-fit results if possible
            prefit_source = source / "prefit_results.pkl"
            if prefit_source.exists():
                self.prefit_results[source] = PrefitResults.from_file(
                    prefit_source
                )

        return None

    def reshape_residual_history(
        self, estimation: EstimationResults
    ) -> list[np.ndarray]:

        n_epochs = len(estimation.epochs)
        residual_size = int(len(estimation.final_residuals) // n_epochs)
        residual_history = [
            np.array(residual_set).reshape(n_epochs, residual_size).T
            for residual_set in estimation.residual_history.T
        ]

        return residual_history

    def show_residual_history_single_source(self, source: Path) -> None:

        # Ensure estimation results are available
        source_id = self.get_directory_id(source)
        if source not in self.estimation_results:
            log.fatal(f"Estimation results not available for {source_id}")
            exit(1)

        # Define setup for canvas
        canvas_setup = self.generate_canvas_setup(
            "doppler-residual-history.png", save_in_base=False
        )
        canvas_setup.canvas_title = f"Doppler residual history :: {source_id}"

        # Generate figure
        with ng.Mosaic("a;b", canvas_setup) as canvas:

            # Get estimation epochs and reshaped residuals
            estimation = self.estimation_results[source]
            dt = (estimation.epochs - self.ref_epoch) / 3600.0
            residuals = self.reshape_residual_history(estimation)

            # Define a mask for propagation epochs based on observations
            best_propagation = estimation.final_propagation_results
            observation_period_mask = (
                best_propagation.epochs >= np.min(estimation.epochs)
            ) * (best_propagation.epochs <= np.max(estimation.epochs))

            # Get masked propagation epochs and distance to Mars
            dt_prop = (
                best_propagation.epochs[observation_period_mask]
                - self.ref_epoch
            ) / 3600.0
            dmars = nt.CartesianState(*best_propagation.rstate_j2000.T).r_mag[
                observation_period_mask
            ]

            # Define settings for history subplot
            history_setup = self.generate_figure_setup(
                "Residual [mHz]", right_axis=True
            )

            # Generate residual-history subplot
            with canvas.subplot(
                history_setup, generator=ng.DoubleAxis
            ) as history_subfig:

                # Add distance to Mars in right axis
                history_subfig.line(
                    dt_prop,
                    dmars * 1e-7,
                    color="black",
                    alpha=0.2,
                    axis="right",
                )

                # Add residuals for current iteration
                for iteration, residual_set in enumerate(residuals):

                    history_subfig.line(
                        dt,
                        residual_set[0] * self.residual_scale,
                        fmt=".",
                        markersize=2,
                        label=f"Iteration {iteration + 1}",
                    )

            # Calculate statistics for final residual
            final_residual = (
                residuals[estimation.best_iteration_index][0]
                * self.residual_scale
            )
            mu_final = np.mean(final_residual)
            std_final = np.abs(np.std(final_residual))
            final_label = (
                rf"$\mu \pm \sigma$ = {mu_final:.2f} $\pm$ {std_final:.2f} mHz"
            )

            # Define settings for final-residual subplot
            final_setup = self.generate_figure_setup(
                "Final residuals [mHz]", right_axis=True
            )
            final_setup.axtitle = f"Best iteration :: {final_label}"

            # Generate subplot for final residual
            with canvas.subplot(
                final_setup, generator=ng.DoubleAxis
            ) as final_subfig:

                # Add distance to Mars in right axis
                final_subfig.line(
                    dt_prop,
                    dmars * 1e-7,
                    color="black",
                    alpha=0.2,
                    axis="right",
                )

                # Add final residual for current source
                final_subfig.line(
                    dt,
                    residuals[estimation.best_iteration_index][0]
                    * self.residual_scale,
                    fmt=".",
                    markersize=2,
                )

        return None

    def show_residual_history(self) -> None:

        for source in self.complete_sources:

            self.show_residual_history_single_source(source)

        return None

    def show_complete_estimation_results_grouped(self) -> None:

        # Define settings for canvas
        canvas_setup = self.generate_canvas_setup(
            "doppler-estimation-results.png",
            save_in_base=True,
        )
        canvas_setup.canvas_title = (
            f"Estimation results :: {self.user_input.base_dir.name}"
        )

        with ng.Mosaic("ab;ab;ab;cd;cd;cd;ee", canvas_setup) as canvas:

            # Initialize pre-fit subfigure
            prefit_setup = self.generate_figure_setup(
                "Pre-fit residuals [mHz]", right_axis=False
            )
            prefit_subfig = canvas.subplot(prefit_setup)

            # Initialize post-fit subfigure
            postfit_setup = self.generate_figure_setup(
                "Post-fit residuals [mHz]", right_axis=False
            )
            postfit_subfig = canvas.subplot(postfit_setup)

            # Initialize orbit subfigures
            dr_setup = self.generate_figure_setup(
                r"$|| \mathbf{r_{est}} - \mathbf{r_{eph}} ||$ [km]",
                right_axis=False,
            )
            dv_setup = self.generate_figure_setup(
                r"$|| \mathbf{v_{est}} - \mathbf{v_{eph}} ||$ [m/s]",
                right_axis=False,
            )
            dr_subfig = canvas.subplot(dr_setup)
            dv_subfig = canvas.subplot(dv_setup)

            # Populate figure with sources
            for source in self.complete_sources:

                # Ensure estimation and pre-fit results are available
                source_id = self.get_directory_id(source)
                if source not in self.estimation_results:
                    log.fatal(f"Estimation results not found for {source_id}")
                    exit(1)
                if source not in self.prefit_results:
                    log.fatal(f"Pre-fit results not found for {source_id}")
                    exit(1)

                # Get estimation, propagation, and pre-fit results
                estimation = self.estimation_results[source]
                propagation = estimation.final_propagation_results
                prefit = self.prefit_results[source]

                # Update pre-fit subfigure
                dt_prefit = (prefit.epochs - self.ref_epoch) / 3600.0
                prefit_subfig.line(
                    dt_prefit,
                    prefit.residual * self.residual_scale,
                    fmt=".",
                    markersize=2,
                    label=source_id,
                )

                # Update post-fit subfigure
                dt_postfit = (estimation.epochs - self.ref_epoch) / 3600.0
                postfit_subfig.line(
                    dt_postfit,
                    estimation.final_residuals * self.residual_scale,
                    fmt=".",
                    markersize=2,
                )

                # Update position and velocity residual subfigures
                dt_prop = (propagation.epochs - self.ref_epoch) / 3600.0
                cstate = nt.CartesianState(*propagation.cstate_j2000.T)
                rstate = nt.CartesianState(*propagation.rstate_j2000.T)
                dstate = cstate - rstate

                dr_subfig.line(dt_prop, dstate.r_mag * 1e-3)
                dv_subfig.line(dt_prop, dstate.v_mag)

            # Add legend subfigure
            legend_setup = ng.PlotSetup(legend_columns=4)
            with canvas.subplot(legend_setup, ng.Legend) as legend:

                legend.add_legend(prefit_subfig)

            # Post-process subfigures
            for subfig in (prefit_subfig, postfit_subfig, dr_subfig, dv_subfig):

                subfig.common_postprocessing()
                subfig.custom_postprocessing()

        return None

    def show_complete_estimation_results_single(self) -> None:

        for source in self.complete_sources:

            # Define settings for canvas
            canvas_setup = self.generate_canvas_setup(
                "doppler-estimation-results.png",
                save_in_base=False,
            )
            source_id = self.get_directory_id(source)
            canvas_setup.canvas_title = f"Estimation results :: {source_id}"

            # Ensure estimation and prefit results are available
            if source not in self.estimation_results:
                log.fatal(f"Estimation results not found for {source_id}")
                exit(1)
            if source not in self.prefit_results:
                log.fatal(f"Pre-fit results not found for {source_id}")
                exit(1)

            # Get estimation, propagation, and pre-fit results
            estimation = self.estimation_results[source]
            propagation = estimation.final_propagation_results
            prefit = self.prefit_results[source]

            # Get dt for prefit, postfit and propagation
            dt_prefit = (prefit.epochs - self.ref_epoch) / 3600.0
            dt_postfit = (estimation.epochs - self.ref_epoch) / 3600.0
            dt_prop = (propagation.epochs - self.ref_epoch) / 3600.0

            # Calculate residual between propagation and ephemerides
            cstate = nt.CartesianState(*propagation.cstate_j2000.T)
            rstate = nt.CartesianState(*propagation.rstate_j2000.T)
            dstate = cstate - rstate

            with ng.Mosaic("ab;cd", canvas_setup) as canvas:

                # Generate pre-fit subfigure
                prefit_setup = self.generate_figure_setup(
                    "Pre-fit residuals [mHz]", right_axis=False
                )
                with canvas.subplot(prefit_setup) as prefit_subfig:
                    prefit_subfig.line(
                        dt_prefit,
                        prefit.residual * self.residual_scale,
                        fmt=".",
                        markersize=2,
                    )

                # Generate post-fit subfigure
                postfit_setup = self.generate_figure_setup(
                    "Post-fit residuals [mHz]", right_axis=False
                )
                with canvas.subplot(postfit_setup) as postfit_subfig:
                    postfit_subfig.line(
                        dt_postfit,
                        estimation.final_residuals * self.residual_scale,
                        fmt=".",
                        markersize=2,
                    )

                # Generate position residual subfigure
                dr_setup = self.generate_figure_setup(
                    r"$|| \mathbf{r_{est}} - \mathbf{r_{eph}} ||$ [km]",
                    right_axis=False,
                )
                with canvas.subplot(dr_setup) as dr_subfig:
                    dr_subfig.line(dt_prop, dstate.r_mag * 1e-3)

                dv_setup = self.generate_figure_setup(
                    r"$|| \mathbf{v_{est}} - \mathbf{v_{eph}} ||$ [m/s]",
                    right_axis=False,
                )
                with canvas.subplot(dv_setup) as dv_subfig:
                    dv_subfig.line(dt_prop, dstate.v_mag)

        return None

    def show_complete_estimation_results(self) -> None:

        if self.user_input.group_in_base and (len(self.complete_sources) > 1):
            self.show_complete_estimation_results_grouped()
        else:
            self.show_complete_estimation_results_single()

        return None
