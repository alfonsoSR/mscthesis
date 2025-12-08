from ...io import PrefitResults
from ...io.command_line.core import CLInputFigure
from ..core import AnalysisFigureManagerBase
from nastro import graphics as ng
from pathlib import Path


class Prefit1DFigureManager(AnalysisFigureManagerBase[CLInputFigure]):

    def __init__(self, user_input: CLInputFigure) -> None:

        super().__init__(user_input)

        # Load results from pre-fit analysis
        self.results: dict[Path, PrefitResults] = {
            source_dir: PrefitResults.from_file(
                source_dir / "prefit_results.pkl"
            )
            for source_dir in self.user_input.source_dirs
        }

        return None

    def update_figure(
        self, figure: ng.BaseFigure, x, y, label: str | None = None
    ) -> None:

        figure.line(x, y, fmt=".", markersize=2, label=label)

        return None


def show_doppler_prefit_residuals(user_input: "CLInputFigure") -> None:

    # Initialize manager
    manager = Prefit1DFigureManager(user_input)

    # Define figure setup
    output_dir = (
        user_input.source_dir
        if len(user_input.source_dirs) == 1
        else user_input.base_dir
    )
    figure_setup = ng.PlotSetup(
        canvas_title="Pre-fit residuals",
        canvas_size=(6, 3),
        xlabel=f"Hours since {manager.ref_epoch_isot}",
        ylabel="Frequency residual [mHz]",
        legend_columns=3,
        scilimits_y=(-3, 3),
        show=user_input.show,
        save=user_input.save,
        dir=output_dir,
        name=f"doppler-prefits.png",
    )

    # Generate figure
    with ng.SingleAxis(figure_setup) as fig:

        for source_dir in user_input.source_dirs:

            # Load prefit results
            prefits = PrefitResults.from_file(source_dir / "prefit_results.pkl")

            # Define epochs and residuals
            dt = (prefits.epochs - manager.ref_epoch) / 3600.0
            residuals = prefits.residual * 1e3

            # Define label from source directory
            source_label = manager.get_directory_id(source_dir)

            # Add data to figure
            manager.update_figure(fig, dt, residuals, source_label)


# def show_doppler_prefit_residuals(user_input: "CLInputFigure") -> None:

#     # Initialize manager
#     manager = Prefit1DFigureManager(user_input)

#     # Define setup for canvas and subfigures
#     canvas_setup = manager.generate_canvas_setup("prefit-residuals-doppler.png")
#     canvas_setup = canvas_setup.version(
#         canvas_title="Closed-loop pre-fit analysis"
#     )
#     obs_setup = manager.generate_figure_setup(
#         ylabel="Observable [Hz]",
#         right_axis=False,
#     )
#     residual_setup = manager.generate_figure_setup(
#         ylabel="Residual [mHz]",
#         right_axis=False,
#     )

#     # Generate figure
#     with ng.Mosaic("a;b", canvas_setup) as canvas:

#         observations_subfig = canvas.subplot(obs_setup)
#         residuals_subfig = canvas.subplot(residual_setup)

#         # Add results from sources
#         for source_dir, source_data in manager.results.items():

#             # Get id for current source
#             id = manager.get_directory_id(source_dir)

#             # Calculate dt for current source
#             dt = (source_data.epochs - manager.ref_epoch) / 3600.0

#             # If first source, add measurements to observations' figure
#             if source_dir == manager.source_dir:

#                 manager.update_figure(
#                     observations_subfig,
#                     dt,
#                     source_data.measured,
#                     label="Observations",
#                 )

#             # Add observations
#             manager.update_figure(
#                 observations_subfig,
#                 dt,
#                 source_data.simulated,
#                 label=f"Simulated {id}",
#             )
#             manager.update_figure(
#                 residuals_subfig, dt, source_data.residual * 1e3, label=id
#             )

#         # Post-process figures
#         for subfig in (observations_subfig, residuals_subfig):

#             subfig.custom_postprocessing()
#             subfig.common_postprocessing()

#     return None


def show_doppler_prefit_residuals_with_observations(
    user_input: "CLInputFigure",
) -> None:

    # Initialize manager
    manager = Prefit1DFigureManager(user_input)

    # Define setup for canvas and subfigures
    canvas_setup = manager.generate_canvas_setup("prefit-residuals-doppler.png")
    canvas_setup = canvas_setup.version(
        canvas_title="Closed-loop pre-fit analysis"
    )
    obs_setup = manager.generate_figure_setup(
        ylabel="Observable [Hz]",
        right_axis=False,
    )
    residual_setup = manager.generate_figure_setup(
        ylabel="Residual [mHz]",
        right_axis=False,
    )

    # Generate figure
    with ng.Mosaic("a;b", canvas_setup) as canvas:

        observations_subfig = canvas.subplot(obs_setup)
        residuals_subfig = canvas.subplot(residual_setup)

        # Add results from sources
        for source_dir, source_data in manager.results.items():

            # Get id for current source
            id = manager.get_directory_id(source_dir)

            # Calculate dt for current source
            dt = (source_data.epochs - manager.ref_epoch) / 3600.0

            # If first source, add measurements to observations' figure
            if source_dir == manager.source_dir:

                manager.update_figure(
                    observations_subfig,
                    dt,
                    source_data.measured,
                    label="Observations",
                )

            # Add observations
            manager.update_figure(
                observations_subfig,
                dt,
                source_data.simulated,
                label=f"Simulated {id}",
            )
            manager.update_figure(
                residuals_subfig, dt, source_data.residual * 1e3, label=id
            )

        # Post-process figures
        for subfig in (observations_subfig, residuals_subfig):

            subfig.custom_postprocessing()
            subfig.common_postprocessing()

    return None
