"""Residuals

Functionality to plot pre-fit and post-fit residuals
"""

from pathlib import Path
from nastro import graphics as ng
from ...io import PrefitResults, EstimationResults
from .core import VisualizationManager, ComparisonManager
import numpy as np
from dataclasses import dataclass


@dataclass
class FigureData:

    x: np.ndarray
    y: np.ndarray
    fmt: str = "."
    markersize: int = 3
    color: str | None = None
    label: str | None = None
    alpha: float = 0.7

    def __sub__(self, other: "FigureData") -> "FigureData":

        # Fail if epochs are not the same
        if len(self.x) != len(other.x):
            raise ValueError(
                "Attempted to subtract data with different lengths"
            )
        if not all(np.isclose(self.x, other.x, rtol=0, atol=1e-12)):
            raise ValueError("Attempted to subtract data non-coincident epochs")

        values = self.__dict__.copy()
        values["y"] = self.y - other.y
        return FigureData(**values)


def pre_fit_residuals_from_source_dir(
    source_dir: Path, ref_epoch: float
) -> FigureData:

    # Load pre-fit results
    prefit = PrefitResults.from_file(source_dir / "prefit_results.pkl")

    # Extract epochs and residuals
    dt = (prefit.epochs - ref_epoch) / 3600.0
    residual = prefit.residual * 1e3

    return FigureData(dt, residual)


def post_fit_residuals_from_source_dir(
    source_dir: Path, ref_epoch: float
) -> FigureData:

    # Load pre-fit results
    estimation = EstimationResults.from_file(source_dir / "estimation.pkl")

    # Extract epochs and residuals
    dt = (estimation.epochs - ref_epoch) / 3600.0
    residual = (estimation.final_residuals * 1e3).reshape(len(dt), 1).T[0]

    return FigureData(dt, residual)


def post_fit_residual_history_from_source_dir(
    source_dir: Path, ref_epoch: float
) -> list[FigureData]:

    # Load pre-fit results
    estimation = EstimationResults.from_file(source_dir / "estimation.pkl")

    # Extract epochs and residuals
    dt = (estimation.epochs - ref_epoch) / 3600.0

    return [
        FigureData(
            x=dt, y=(np.array(current_residuals) * 1e3).reshape(len(dt), 1).T[0]
        )
        for current_residuals in estimation.residual_history.T
    ]


def update_figure[T: ng.BaseFigure](fig: T, data: FigureData) -> T:

    fig.line(
        data.x,
        data.y,
        fmt=data.fmt,
        color=data.color,
        markersize=data.markersize,
        label=data.label,
        alpha=data.alpha,
    )

    return fig


class DopplerVisualizationManager(VisualizationManager):

    # def execute_function(self, function_name: str) -> None:

    #     if (not self._single_source) and self.group:
    #         getattr(self, function_name)(self.source_dirs)
    #     else:
    #         for source_dir in self.source_dirs:
    #             getattr(self, function_name)([source_dir])

    #     return None

    def pre_fit_residuals(self, sources: list[Path]) -> None:

        # Define settings for figure
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="pre-fit-residuals.png",
            title="Pre-fit residuals",
            canvas_size=(6, 4),
        )
        figure_setup = ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel="Pre-fit residual [mHz]",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        with ng.Mosaic("a;a;a;a;b", canvas_setup) as canvas:

            with canvas.subplot(figure_setup) as fig:

                for source in sources:

                    data = pre_fit_residuals_from_source_dir(
                        source, self.ref_epoch
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(fig, data)

            with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

                for source in sources:

                    leg.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def post_fit_residuals(self, sources: list[Path]) -> None:

        # Define settings for figure
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="post-fit-residuals.png",
            title="Post-fit residuals",
            canvas_size=(6, 4),
        )
        figure_setup = ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel="Post-fit residual [mHz]",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data per source
        data_per_source = {}
        for source in sources:
            data = post_fit_residuals_from_source_dir(source, self.ref_epoch)
            data.label = rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
            data_per_source[source] = data

        with ng.Mosaic("a;a;a;a;c", canvas_setup) as canvas:

            with canvas.subplot(figure_setup) as fig:
                for source in sources:
                    update_figure(fig, data_per_source[source])

            with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

                for source in sources:

                    leg.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def residuals(self, sources: list[Path]) -> None:

        # Define settings for canvas
        figure_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="doppler-residuals.png",
            title=f"Doppler residuals",
            canvas_size=(6, 6.5),
        )

        # Define settings for subfigures
        prefit_setup = ng.PlotSetup(
            ylabel="Pre-fit residual [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        postfit_setup = ng.PlotSetup(
            ylabel="Post-fit residual [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Generate figure
        with ng.Mosaic("a;a;a;a;b;b;b;b;c", figure_setup) as canvas:

            with canvas.subplot(prefit_setup) as prefit_fig:
                for source in sources:
                    data = pre_fit_residuals_from_source_dir(
                        source, self.ref_epoch
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(prefit_fig, data)

            with canvas.subplot(postfit_setup) as postfit_fig:
                for source in sources:
                    data = post_fit_residuals_from_source_dir(
                        source, self.ref_epoch
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(postfit_fig, data)

            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def residuals_with_distribution(self, sources: list[Path]) -> None:

        # Define settings for canvas
        figure_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="doppler-residuals.png",
            title=f"Doppler residuals",
            canvas_size=(6, 6.5),
        )

        # Define settings for subfigures
        prefit_setup = ng.PlotSetup(
            ylabel="Pre-fit residual [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        postfit_setup = ng.PlotSetup(
            ylabel="Post-fit residual [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        distribution_setup = ng.PlotSetup(
            show_tick_labels_y=False,
            xlabel="PDF",
            xlim=(0, 0.1),
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Pre-fit and post-fit data per source
        pre_fit_per_source: dict[Path, FigureData] = {}
        post_fit_per_source: dict[Path, FigureData] = {}
        for source in sources:

            # Pre-fit residuals
            pre_fit_data = pre_fit_residuals_from_source_dir(
                source, self.ref_epoch
            )
            pre_fit_data.label = (
                rf"{np.average(pre_fit_data.y):.2f}"
                rf" $\pm$ {np.std(pre_fit_data.y):.2f}"
            )
            pre_fit_per_source[source] = pre_fit_data

            # Post-fit residuals
            post_fit_data = post_fit_residuals_from_source_dir(
                source, self.ref_epoch
            )
            post_fit_data.label = (
                rf"{np.average(post_fit_data.y):.2f}"
                rf" $\pm$ {np.std(post_fit_data.y):.2f}"
            )
            post_fit_per_source[source] = post_fit_data

        # Generate figure
        mosaic_str = "aaab;aaab;aaab;aaab;cccd;cccd;cccd;cccd;eeee"
        with ng.Mosaic(mosaic_str, figure_setup) as canvas:

            with canvas.subplot(prefit_setup) as prefit_fig:
                for source in sources:
                    update_figure(prefit_fig, pre_fit_per_source[source])

            with canvas.subplot(distribution_setup) as prefit_dist:
                for source in sources:
                    prefit_dist.hist(
                        pre_fit_per_source[source].y,
                        bins=20,
                        normalize=True,
                        orientation="horizontal",
                    )

            with canvas.subplot(postfit_setup) as postfit_fig:
                for source in sources:
                    update_figure(postfit_fig, post_fit_per_source[source])

            with canvas.subplot(distribution_setup) as postfit_dist:
                for source in sources:
                    postfit_dist.hist(
                        post_fit_per_source[source].y,
                        bins=20,
                        normalize=True,
                        orientation="horizontal",
                    )

            with canvas.subplot(legend_setup, ng.Legend) as legend:
                for source in sources:
                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def post_fit_residuals_with_distribution(self, sources: list[Path]) -> None:

        # Define settings for figure
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="post-fit-residuals.png",
            title="Post-fit residuals",
            canvas_size=(6, 4),
        )
        figure_setup = ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel="Post-fit residual [mHz]",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        histo_setup = ng.PlotSetup(
            show_tick_labels_y=False, xlabel="PDF", xlim=(0, 0.1)
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data per source
        data_per_source = {}
        for source in sources:
            data = post_fit_residuals_from_source_dir(source, self.ref_epoch)
            data.label = rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
            data_per_source[source] = data

        with ng.Mosaic("aaab;aaab;aaab;aaab;cccc", canvas_setup) as canvas:

            with canvas.subplot(figure_setup) as fig:
                for source in sources:
                    update_figure(fig, data_per_source[source])

            with canvas.subplot(histo_setup) as histo:
                for source_data in data_per_source.values():

                    histo.hist(
                        source_data.y,
                        bins=20,
                        orientation="horizontal",
                        normalize=True,
                        cumulative=False,
                    )

            with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

                for source in sources:

                    leg.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def pre_fit_residuals_with_distribution(self, sources: list[Path]) -> None:

        # Define settings for figure
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="pre-fit-residuals.png",
            title="Pre-fit residuals",
            canvas_size=(6, 4),
        )
        figure_setup = ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel="Pre-fit residual [mHz]",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        histo_setup = ng.PlotSetup(
            show_tick_labels_y=False,
            xlabel="PDF",
            xlim=(0, 0.1),
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data per source
        data_per_source = {}
        for source in sources:
            data = pre_fit_residuals_from_source_dir(source, self.ref_epoch)
            data.label = rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
            data_per_source[source] = data

        with ng.Mosaic("aaab;aaab;aaab;aaab;cccc", canvas_setup) as canvas:

            with canvas.subplot(figure_setup) as fig:
                for source in sources:
                    update_figure(fig, data_per_source[source])

            with canvas.subplot(histo_setup) as histo:
                for source_data in data_per_source.values():
                    histo.hist(
                        source_data.y,
                        bins=20,
                        orientation="horizontal",
                        normalize=True,
                    )

            with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

                for source in sources:

                    leg.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def post_fit_residual_history(self, sources: list[Path]) -> None:

        if len(sources) > 1:
            raise ValueError("Residual history only with one source at a time")

        # Define settings for figure
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="post-fit-residual-history.png",
            title=(
                f"Post-fit residual history ::"
                f" {self.get_id_for_directory(sources[0])}"
            ),
            canvas_size=(6, 4),
        )
        figure_setup = ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel="Post-fit residual [mHz]",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data for source
        data_series = post_fit_residual_history_from_source_dir(
            sources[0], self.ref_epoch
        )
        for idx, data in enumerate(data_series):
            data_series[idx].label = (
                rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
            )

        with ng.Mosaic("a;a;a;a;c", canvas_setup) as canvas:

            with canvas.subplot(figure_setup) as fig:
                for data in data_series:
                    update_figure(fig, data)

            with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

                for idx, _ in enumerate(data_series):
                    leg.line(-1, -1, fmt="o", label=f"Step {idx + 1}")

        return None


class DopplerComparisonManager(ComparisonManager):

    # def execute_function(self, function_name: str) -> None:

    #     if (not self._single_source) and self.group:
    #         getattr(self, function_name)(self.source_dirs)
    #     else:
    #         for source_dir in self.source_dirs:
    #             getattr(self, function_name)([source_dir])

    #     return None

    # def pre_fit_residuals(self, sources: list[Path]) -> None:

    #     # Get ID for reference
    #     ref_id = self.get_id_for_directory(self.ref_dir)

    #     # Load data for reference
    #     reference = pre_fit_residuals_from_source_dir(
    #         self.ref_dir, self.ref_epoch
    #     )

    #     # Define settings for figure
    #     canvas_setup = self.generate_default_canvas_setup(
    #         source=self.outdir,
    #         figure_name="delta-pre-fit-residuals.png",
    #         title=f"Difference wrt. {ref_id}",
    #         canvas_size=(6, 4),
    #     )
    #     figure_setup = ng.PlotSetup(
    #         xlabel=f"Hours past {self.ref_epoch_isot}",
    #         ylabel=r"Pre-fit :: $\rho - \rho_{ref}$ [mHz]",
    #         legend_title=r"$\mu \pm \sigma$ [mHz]",
    #         legend_location=self.legend_location,
    #     )
    #     legend_setup = ng.PlotSetup(
    #         legend_columns=self.legend_cols,
    #     )

    #     with ng.Mosaic("a;a;a;a;b", canvas_setup) as canvas:

    #         with canvas.subplot(figure_setup) as fig:

    #             for source in sources:

    #                 data = (
    #                     pre_fit_residuals_from_source_dir(
    #                         source, self.ref_epoch
    #                     )
    #                     - reference
    #                 )
    #                 data.label = (
    #                     rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
    #                 )
    #                 update_figure(fig, data)

    #         with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

    #             for source in sources:

    #                 leg.line(
    #                     -1, -1, fmt="o", label=self.get_id_for_directory(source)
    #                 )

    #     return None

    # def post_fit_residuals(self, sources: list[Path]) -> None:

    #     # Get ID for reference
    #     ref_id = self.get_id_for_directory(self.ref_dir)

    #     # Load data for reference
    #     reference = post_fit_residuals_from_source_dir(
    #         self.ref_dir, self.ref_epoch
    #     )

    #     # Define settings for figure
    #     canvas_setup = self.generate_default_canvas_setup(
    #         source=self.outdir,
    #         figure_name="delta-post-fit-residuals.png",
    #         title=f"Difference wrt. {ref_id}",
    #         canvas_size=(6, 4),
    #     )
    #     figure_setup = ng.PlotSetup(
    #         xlabel=f"Hours past {self.ref_epoch_isot}",
    #         ylabel=r"Post-fit :: $\rho - \rho_{ref}$ [mHz]",
    #         legend_title=r"$\mu \pm \sigma$ [mHz]",
    #         legend_location=self.legend_location,
    #     )
    #     legend_setup = ng.PlotSetup(
    #         legend_columns=self.legend_cols,
    #     )

    #     with ng.Mosaic("a;a;a;a;b", canvas_setup) as canvas:

    #         with canvas.subplot(figure_setup) as fig:

    #             for source in sources:

    #                 data = (
    #                     post_fit_residuals_from_source_dir(
    #                         source, self.ref_epoch
    #                     )
    #                     - reference
    #                 )
    #                 data.label = (
    #                     rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
    #                 )
    #                 update_figure(fig, data)

    #         with canvas.subplot(setup=legend_setup, generator=ng.Legend) as leg:

    #             for source in sources:

    #                 leg.line(
    #                     -1, -1, fmt="o", label=self.get_id_for_directory(source)
    #                 )

    #     return None

    def residuals(self, sources: list[Path]) -> None:

        # Get ID for reference
        ref_id = self.get_id_for_directory(self.ref_dir)

        # Load data for reference
        prefit_ref = pre_fit_residuals_from_source_dir(
            self.ref_dir, self.ref_epoch
        )
        postfit_ref = post_fit_residuals_from_source_dir(
            self.ref_dir, self.ref_epoch
        )

        # Define settings for canvas
        figure_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="delta-doppler-residuals.png",
            title=f"Difference wrt. {ref_id}",
            canvas_size=(6, 6.5),
        )

        # Define settings for subfigures
        prefit_setup = ng.PlotSetup(
            ylabel=r"Pre-fit :: $\rho - \rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
        )
        postfit_setup = ng.PlotSetup(
            ylabel=r"Post-fit :: $\rho - \rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Generate figure
        with ng.Mosaic("a;a;a;a;b;b;b;b;c", figure_setup) as canvas:

            with canvas.subplot(prefit_setup) as prefit_fig:
                for source in sources:
                    data = (
                        pre_fit_residuals_from_source_dir(
                            source, self.ref_epoch
                        )
                        - prefit_ref
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(prefit_fig, data)

            with canvas.subplot(postfit_setup) as postfit_fig:
                for source in sources:
                    data = (
                        post_fit_residuals_from_source_dir(
                            source, self.ref_epoch
                        )
                        - postfit_ref
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(postfit_fig, data)

            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def pre_fit_residuals(self, sources: list[Path]) -> None:

        # Get ID for reference
        ref_id = self.get_id_for_directory(self.ref_dir)

        # Load data for reference
        prefit_ref = pre_fit_residuals_from_source_dir(
            self.ref_dir, self.ref_epoch
        )

        # Define settings for canvas
        figure_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="delta-pre-fit-residuals.png",
            title=f"Difference wrt. {ref_id}",
            canvas_size=(6, 6.5),
        )

        # Define settings for subfigures
        prefit_setup = ng.PlotSetup(
            ylabel=r"Pre-fit :: $\rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
        )
        delta_setup = ng.PlotSetup(
            ylabel=r"Pre-fit :: $\rho - \rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Generate figure
        with ng.Mosaic("a;a;a;a;b;b;b;b;c", figure_setup) as canvas:

            with canvas.subplot(prefit_setup) as prefit_fig:

                data = prefit_ref
                data.label = (
                    rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                )
                update_figure(prefit_fig, data)

            with canvas.subplot(delta_setup) as postfit_fig:
                for source in sources:
                    data = (
                        pre_fit_residuals_from_source_dir(
                            source, self.ref_epoch
                        )
                        - prefit_ref
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(postfit_fig, data)

            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None

    def post_fit_residuals(self, sources: list[Path]) -> None:

        # Get ID for reference
        ref_id = self.get_id_for_directory(self.ref_dir)

        # Load data for reference
        ref = post_fit_residuals_from_source_dir(self.ref_dir, self.ref_epoch)

        # Define settings for canvas
        figure_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name="delta-post-fit-residuals.png",
            title=f"Difference wrt. {ref_id}",
            canvas_size=(6, 6.5),
        )

        # Define settings for subfigures
        reference_setup = ng.PlotSetup(
            ylabel=r"Post-fit :: $\rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
        )
        delta_setup = ng.PlotSetup(
            ylabel=r"Post-fit :: $\rho - \rho_{ref}$ [mHz]",
            xlabel=f"Hours since {self.ref_epoch_isot}",
            legend_title=r"$\mu \pm \sigma$ [mHz]",
            legend_location=self.legend_location,
            legend_columns=self.legend_cols,
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Generate figure
        with ng.Mosaic("a;a;a;a;b;b;b;b;c", figure_setup) as canvas:

            with canvas.subplot(reference_setup) as reference_fig:

                data = ref
                data.label = (
                    rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                )
                update_figure(reference_fig, data)

            with canvas.subplot(delta_setup) as delta_fig:
                for source in sources:
                    data = (
                        post_fit_residuals_from_source_dir(
                            source, self.ref_epoch
                        )
                        - ref
                    )
                    data.label = (
                        rf"{np.average(data.y):.2f} $\pm$ {np.std(data.y):.2f}"
                    )
                    update_figure(delta_fig, data)

            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None
