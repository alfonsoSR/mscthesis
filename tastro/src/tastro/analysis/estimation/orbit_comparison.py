from ..core import AnalysisFigureManagerBase
from ...io.command_line.comparison import CLInputFigure
from ...io import EstimationResults
from pathlib import Path
from nastro import graphics as ng, types as nt
from ..utils import get_propagation_start_epoch_from_config
from ...logging import log
import traceback
import numpy as np


class CartesianStateUncertainty[T: nt.Vector](nt.CartesianState):

    @staticmethod
    def from_covariance_history(
        covariance_history: dict[float, np.ndarray],
    ) -> "CartesianStateUncertainty":

        # Extract covariance matrix as function of time: (n_epochs, m, m)
        covariances = np.array(list(covariance_history.values()))

        # Extract standard deviation as function of time: (m, n_epochs)
        standard_deviations = np.sqrt(covariances.diagonal(0, 1, 2).T)

        return CartesianStateUncertainty(*standard_deviations)


class OrbitComparisonManager(AnalysisFigureManagerBase[CLInputFigure]):

    def __init__(self, user_input: "CLInputFigure", frame: str) -> None:

        super().__init__(user_input)

        # Initialize containers for calculated and reference states
        self.cstates: dict[Path, nt.CartesianState[nt.Vector]] = {}
        self.rstates: dict[Path, nt.CartesianState[nt.Vector]] = {}
        self.uncertainties: dict[
            Path, CartesianStateUncertainty[nt.Vector] | None
        ] = {}
        self.epochs: dict[Path, "np.ndarray"] = {}
        self.est_epochs: dict[Path, "np.ndarray"] = {}

        # Fill containers with results
        for source in self.user_input.source_dirs:

            # Load propagation results
            estimation = EstimationResults.from_file(source / "estimation.pkl")
            results = estimation.final_propagation_results

            # Get propagation output in correct frame
            if frame in ("j2000", "rsw"):
                self.cstates[source] = nt.CartesianState(
                    *getattr(results, f"cstate_{frame}").T
                )
                self.rstates[source] = nt.CartesianState(
                    *getattr(results, f"rstate_{frame}").T
                )
                self.epochs[source] = results.epochs
            else:
                log.fatal(
                    f"Requested propagation output in invalid frame: {frame}"
                )
                log.fatal(traceback.extract_stack()[-2])
                exit(1)

            # Get covariance history
            if estimation.covariance_history is None:
                log.warning(f"Covariance history not available for {source}")
                self.uncertainties[source] = None
                self.est_epochs[source] = estimation.epochs
            else:
                self.uncertainties[source] = (
                    CartesianStateUncertainty.from_covariance_history(
                        estimation.covariance_history
                    )
                )
                self.est_epochs[source] = estimation.epochs

        # Save frame choice
        self.frame = frame

        return None

    @property
    def vector_components(self) -> list[str]:

        match self.frame:

            case "j2000":
                return ["x", r"\dot x", "y", r"\dot y", "z", r"\dot z"]

            case "rsw":
                return ["r", r"\dot r", "s", r"\dot s", "w", r"\dot w"]

            case _:
                raise NotImplementedError("Unreachable")

    def propagation_residual_wrt_ephemerides_group(self) -> None:

        # Generate settings for canvas
        canvas_setup = self.generate_canvas_setup(
            f"estimated-orbit-vs-ephemerides-{self.frame}.png",
            save_in_base=True,
        )
        canvas_setup.canvas_title = (
            f"Estimated orbit vs ephemerides :: {self.frame.upper()}"
        )

        # Generate common setup for subfigures
        subfig_setup = self.generate_figure_setup("", right_axis=False)
        legend_setup = ng.PlotSetup(legend_columns=2)

        # Initialize canvas
        with ng.Mosaic("ab;ab;ab;cd;cd;cd;ef;ef;ef;gg", canvas_setup) as canvas:

            components = ["x", "dx", "y", "dy", "z", "dz"]
            scales: list[float] = []

            # Initialize subfigures
            subfigures: dict[str, ng.BaseFigure] = {}
            for idx, label in enumerate(self.vector_components):

                # Generate settings for current subplot
                unit = "mm/s" if ((idx % 2) != 0) else "m"
                scales.append(1e3 if ((idx % 2) != 0) else 1)
                current_setup = subfig_setup.version(
                    ylabel=rf"$\Delta {label}$ [{unit}]"
                )

                # Add subfigure to list
                subfigures[components[idx]] = canvas.subplot(current_setup)

            # Fill subfigures with data
            dmars_present: bool = False
            for source, epochs in self.epochs.items():

                # Get source id
                source_id = self.get_directory_id(source)

                # Get dt, residuals and distance to Mars
                dt = (epochs - self.ref_epoch) / 3600.0
                dt_est = (self.est_epochs[source] - self.ref_epoch) / 3600.0
                residual = self.cstates[source] - self.rstates[source]
                uncertainty = self.uncertainties[source]
                dmars = self.rstates[source].r_mag * 1e-7

                for idx, component in enumerate(subfigures):

                    # Get subfigure and scale
                    subfig = subfigures[component]
                    scale = scales[idx]

                    # Only add legend to x subfigure
                    subfig.line(
                        dt,
                        getattr(residual, component) * scale,
                        label=source_id if component == "x" else None,
                    )
                    # if uncertainty is not None:
                    #     subfig.line(dt_est, 0.0 * dt_est, color="w", alpha=0)
                    #     subfig.boundary(
                    #         getattr(uncertainty, component) * scale, alpha=0.7
                    #     )
                    # if not dmars_present:
                    #     subfig.line(
                    #         dt, dmars, axis="right", alpha=0.2, color="black"
                    #     )

                dmars_present = True

            # Add legend
            with canvas.subplot(legend_setup, ng.Legend) as legend:
                legend.add_legend(subfigures["x"])

            # Post-process subfigures
            for subfig in subfigures.values():
                subfig.common_postprocessing()
                subfig.custom_postprocessing()

        return None

    def propagation_residual_wrt_ephemerides_single(self) -> None:

        # Loop over sources
        for source, epochs in self.epochs.items():

            # Generate settings for canvas
            canvas_setup = self.generate_canvas_setup(
                f"estimated-orbit-vs-ephemerides-{self.frame}.png",
                save_in_base=False,
            )
            canvas_setup.canvas_title = (
                f"Estimated orbit vs ephemerides :: {self.frame.upper()} "
                f":: {self.get_directory_id(source)}"
            )

            # Generate common settings for subfigures
            subfig_setup = self.generate_figure_setup("")

            # Get dt, residual and distance to Mars
            dt = (epochs - self.ref_epoch) / 3600.0
            residual = self.cstates[source] - self.rstates[source]
            dmars = self.rstates[source].r_mag * 1e-7

            # Generate figure
            with ng.Mosaic("ab;cd;ef", canvas_setup) as canvas:

                components = ["x", "dx", "y", "dy", "z", "dz"]

                for idx, label in enumerate(self.vector_components):

                    # Generate settings for current subplot
                    unit = "mm/s" if ((idx % 2) != 0) else "m"
                    scale = 1e3 if ((idx % 2) != 0) else 1
                    current_setup = subfig_setup.version(
                        ylabel=rf"$\Delta {label}$ [{unit}]"
                    )

                    # Generate subfigure
                    with canvas.subplot(current_setup, ng.DoubleAxis) as subfig:

                        subfig.line(
                            dt, getattr(residual, components[idx]) * scale
                        )
                        subfig.line(
                            dt, dmars, axis="right", alpha=0.2, color="black"
                        )

        return None

    def propagation_residual_wrt_ephemerides(self) -> None:

        if self.user_input.group_in_base:
            self.propagation_residual_wrt_ephemerides_group()
        else:
            self.propagation_residual_wrt_ephemerides_single()

        return None


# def show_orbits(user_input: "CLInputFigure") -> None:

#     # Initialize manager
#     manager = OrbitComparisonManager(user_input)

#     with ng.CompareRswStates() as figure:

#         # Loop over results
#         for source, results in manager.estimation_results.items():

#             # Get orbit for best iteration
#             simulation = results.final_propagation_results

#             # Calculate dt and difference wrt ephemerides
#             dt = (simulation.epochs - manager.ref_epoch) / 3600.0
#             cstate = nt.CartesianState(*simulation.cstate_rsw.T)
#             rstate = nt.CartesianState(*simulation.rstate_rsw.T)
#             residual = cstate - rstate

#             # Add residual to figure
#             figure.compare_states(
#                 dt,
#                 cstate,
#                 rstate,
#                 is_dt=True,
#                 label=manager.get_directory_id(source),
#             )

#     return None
