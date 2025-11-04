from ..core import AnalysisFigureManagerBase
from ...io.command_line.core import CLInputFigure
from ...logging import log
from ...io import PropagationOutput
from ...config import CaseSetup
from tudatpy.astro import time_representation as ttime
from nastro import graphics as ng, types as nt
import numpy as np
from scipy import interpolate as scint


class FixedIntegratorAnalysisManager(AnalysisFigureManagerBase[CLInputFigure]):

    def __init__(self, user_input: CLInputFigure) -> None:

        super().__init__(user_input)

        # Load each propagation result with its associated step size
        current_level = log.getEffectiveLevel()
        log.setLevel("ERROR")

        self.results: dict[str, dict[float, PropagationOutput]] = {}
        for source_dir in self.user_input.source_dirs:

            # Read step size and integrator from configuration
            config = CaseSetup.from_config_file(
                source_dir / "configuration.yaml"
            )
            rkf_fixed_setup = config.propagation.integrator.rkf_fixed
            step = rkf_fixed_setup.step_size.to_float()
            integrator = rkf_fixed_setup.coefficients.name

            if integrator not in self.results:
                self.results[integrator] = {}

            self.results[integrator][step] = PropagationOutput.from_file(
                source_dir / "results.pkl"
            )

        log.setLevel(current_level)

        return None

    def get_residuals_and_epochs(
        self, integrator: str, ref_step: float, current_step: float
    ) -> tuple[nt.CartesianState, np.ndarray]:

        # Load results
        r_results = self.results[integrator][ref_step]
        c_results = self.results[integrator][current_step]

        # Fail if steps are not half
        if not np.isclose(ref_step, current_step / 2):
            raise ValueError("Steps should follow half rule")

        # Get reference solution at current epochs
        epochs = c_results.epochs

        # Generate interpolator for reference results
        ref_interp = scint.Akima1DInterpolator(
            r_results.epochs,
            r_results.cstate_j2000,
            axis=0,
            extrapolate=False,
        )

        # Calculate residual
        reference = nt.CartesianState(*ref_interp(epochs).T)
        current = nt.CartesianState(*c_results.cstate_j2000.T)
        residual = current - reference

        return residual, epochs

    def generate_integrator_figure(self, integrator: str) -> None:

        integrator_results = self.results[integrator]
        integrator_steps = sorted(list(integrator_results.keys()))
        ref_step = integrator_steps[0]

        # Define setup for canvas
        canvas_setup = self.generate_canvas_setup(
            f"fixed-integrator-analysis-{integrator}.png"
        )
        canvas_setup = canvas_setup.version(
            canvas_title=f"Integrator analysis :: {integrator}",
            canvas_size=(8, 8),
        )

        # Subplot settings: Error as function of time
        err_setup = ng.PlotSetup(
            yscale="log",
            xlabel=f"Hours past {self.ref_epoch_isot}",
        )
        xerr_setup = err_setup.version(
            ylabel=r"$|| \mathbf{r}(t; \Delta t) - \mathbf{r}(t; \Delta t/2) ||$ [m]",
            ylim=(1e-10, 1e4),
        )
        verr_setup = err_setup.version(
            ylabel=(
                r"$||\mathbf{v}(t; \Delta t) - \mathbf{v}(t; \Delta t/2)||$[m/s]"
            ),
            ylim=(1e-11, 1e1),
        )

        # Subplot settings: Max error as function of step size
        max_err_setup = ng.PlotSetup(
            yscale="log",
            xscale="log",
            xlabel="Step size",
            legend_title=r"$\Delta t$ [s]",
            ylim=None,
        )
        max_xerr_setup = max_err_setup.version(
            ylabel=r"max$(|| \mathbf{r}(t; \Delta t) - \mathbf{r}(t; \Delta t/2) ||)$ [m]",
            ylim=(1e-4, 1e4),
        )
        max_verr_setup = max_err_setup.version(
            ylabel=r"max$(|| \mathbf{v}(t; \Delta t) - \mathbf{v}(t; \Delta t/2) ||)$ [m/s]",
            ylim=(1e-7, 1e1),
        )

        with ng.Mosaic("ab;cd", canvas_setup) as canvas:

            # Initialize subplots
            xerr = canvas.subplot(xerr_setup)
            xerr_max = canvas.subplot(max_xerr_setup)
            verr = canvas.subplot(verr_setup)
            verr_max = canvas.subplot(max_verr_setup)

            for current_step in integrator_steps[1:]:

                # Get residual and epochs
                residual, epochs = self.get_residuals_and_epochs(
                    integrator, ref_step, current_step
                )

                # Calculate dt
                dt = (epochs - self.ref_epoch) / 3600.0

                xerr.line(dt, residual.r_mag)
                verr.line(dt, residual.v_mag)
                xerr_max.line(
                    current_step,
                    np.max(residual.r_mag[np.isfinite(residual.r_mag)]),
                    fmt="o",
                    label=f"{current_step:.2f}",
                )
                verr_max.line(
                    current_step,
                    np.max(residual.v_mag[np.isfinite(residual.v_mag)]),
                    fmt="o",
                    label=f"{current_step:.2f}",
                )

                ref_step = current_step

            # Post-process subfigures
            for subfig in (xerr, verr, xerr_max, verr_max):

                subfig.custom_postprocessing()
                subfig.common_postprocessing()

        return None
