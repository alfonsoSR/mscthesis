from tastro.io.command_line.plotters import DEFAULT_CANVAS_SIZE
from .core import AnalysisManager, ComparisonManager, VisualizationManager
from ...io.command_line.plotters import PlotterInput
from pathlib import Path
from ...io import EstimationResults
import yaml
from nastro import graphics as ng
import numpy as np
from ...io.command_line import (
    ResultComparisonInput,
    ResultVisualizationInput,
    ParameterDistributionInput,
)
from matplotlib import pyplot as plt


def get_labels_of_estimated_parameters(
    parameter_config: dict[str, bool],
) -> list[str]:

    params = ["$x$", "$y$", "$z$", r"$\dot x$", r"$\dot y$", r"$\dot z$"]
    if parameter_config["drag_coefficient"]:
        params.append("$C_D$")
    if parameter_config["radiation_pressure_coefficient"]:
        params.append("$C_R$")
    if parameter_config["k2_radiation_coefficient"]:
        params.append("$K_2$")
    if parameter_config["gm_phobos"]:
        params.append(r"$\mu_{ph}$")
    if parameter_config["gm_deimos"]:
        params.append(r"$\mu_{de}$")
    if parameter_config["arcwise_drag_coefficient"]:
        params += ["$C_{D1}$", "$C_{D2}$", "$C_{D3}$"]

    return params


def get_parameter_labels_from_source(source: Path) -> list[str]:

    # Load parameter configuration
    with (source / "configuration.yaml").open("rb") as buffer:
        __raw_config = yaml.safe_load(buffer)
    __param_config = __raw_config["Estimation"]["parameters"]

    # Add potentially missing parameters
    __potentially_missing = ["gm_deimos", "open_loop_biases"]
    for parameter in __potentially_missing:
        if parameter not in __param_config:
            __param_config[parameter] = False

    # Get list of parameters without biases
    params = get_labels_of_estimated_parameters(__param_config)

    # Handle biases
    if not __param_config["open_loop_biases"]:
        return params

    # Load estimation results
    estimation = EstimationResults.from_file(source / "estimation.pkl")
    biases = len(estimation.final_parameters) - len(params)
    for idx in range(biases):
        params.append(rf"$b_{idx}$")

    return params


class MatrixVisualizationManager(AnalysisManager[ResultVisualizationInput]):

    def generate_default_canvas_setup(
        self,
        source: Path,
        figure_name: str,
        title: str,
        canvas_size: tuple[float, float] = (6, 5),
    ) -> ng.PlotSetup:

        canvas_setup = super().generate_default_canvas_setup(
            source, figure_name, title, canvas_size
        )

        canvas_setup.custom_ticks = True
        canvas_setup.aspect = "equal"
        canvas_setup.grid = False

        return canvas_setup

    def correlation_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")

        # Get parameters to estimate from configuration
        params = get_parameter_labels_from_source(source)

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="correlation-matrix.png",
            title=f"Correlation matrix :: {self.get_id_for_directory(source)}",
            canvas_size=(6, 5),
        )

        quantity = estimation.correlation_matrix

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, vmin=-1, vmax=1, cmap="coolwarm")
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None

    def absolute_correlation_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")

        # Get parameters to estimate from configuration
        params = get_parameter_labels_from_source(source)

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="absolute-correlation-matrix.png",
            title=f"Absolute correlation matrix :: {self.get_id_for_directory(source)}",
            canvas_size=(6, 5),
        )

        quantity = np.abs(estimation.correlation_matrix)

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, vmin=0, vmax=1)
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None

    def anticorrelation_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")

        # Estimated parameters
        params = get_parameter_labels_from_source(source)

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="anticorrelation-matrix.png",
            title=f"Anticorrelation matrix :: {self.get_id_for_directory(source)}",
            canvas_size=(6, 5),
        )

        quantity = 1.0 - np.abs(estimation.correlation_matrix)

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, vmin=0, vmax=1)
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None

    def covariance_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")

        # Get parameters to estimate from configuration
        params = get_parameter_labels_from_source(source)

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="covariance-matrix.png",
            title=(
                f"log10(sqrt(abs(Covariance))) :: "
                f"{self.get_id_for_directory(source)}"
            ),
            canvas_size=(6, 5),
        )

        quantity = np.log10(np.sqrt(np.abs(estimation.covariance_matrix)))

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, cmap="coolwarm", vmin=-6, vmax=6)
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None

    def average_absolute_correlation(self, sources: list[Path]) -> None:

        # Define settings for figure
        self.group = True
        outdir = self.get_output_directory(sources[0])
        canvas_setup = self.generate_default_canvas_setup(
            source=outdir,
            figure_name="average-absolute-correlation.png",
            title=f"Average absolute correlation :: {outdir.name}",
            canvas_size=(6, 5),
        )

        # Iterate over sources
        for idx, source in enumerate(sources):

            # Load absolute correlation from estimation results
            _estimation = EstimationResults.from_file(source / "estimation.pkl")
            abscorr = np.abs(_estimation.correlation_matrix)

            # Get labels of parameters to estimate
            if idx == 0:
                params = get_parameter_labels_from_source(source)
            else:
                __params = get_parameter_labels_from_source(source)
                pass

        return None


class MatrixComparisonManager(ComparisonManager):

    def correlation_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")
        ref_estimation = EstimationResults.from_file(
            self.ref_dir / "estimation.pkl"
        )

        # Get parameters to estimate from configuration
        with (source / "configuration.yaml").open("rb") as buffer:
            raw_config = yaml.safe_load(buffer)
        current_param_config: dict[str, bool] = raw_config["Estimation"][
            "parameters"
        ]

        # Get parameters to estimate in reference
        with (self.ref_dir / "configuration.yaml").open("rb") as buffer:
            ref_raw_config = yaml.safe_load(buffer)
        ref_param_config: dict[str, bool] = ref_raw_config["Estimation"][
            "parameters"
        ]

        # Find largest collection
        param_config = ref_param_config
        if ref_param_config != current_param_config:
            raise ValueError("Cannot compare inconsistent configurations")
        # if len(ref_param_config) >= len(current_param_config):
        #     param_config = ref_param_config
        # else:
        #     param_config = current_param_config

        # Add potentially missing
        potentially_missing = ["gm_deimos"]
        for parameter in potentially_missing:
            if parameter not in param_config:
                param_config[parameter] = False

        # Estimated parameters
        params = get_labels_of_estimated_parameters(param_config)

        # Define title for figure
        canvas_title = (
            f"Correlation difference :: {self.get_id_for_directory(source)} "
            f"vs {self.get_id_for_directory(self.ref_dir)}"
        )
        if len(canvas_title) > 80:
            canvas_title = canvas_title.replace("vs ", r"vs\n")

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="delta-correlation-matrix.png",
            title=canvas_title,
            canvas_size=(6, 5),
        )
        canvas_setup.custom_ticks = True
        canvas_setup.aspect = "equal"
        canvas_setup.grid = False

        # Calculate difference
        quantity = (
            estimation.correlation_matrix - ref_estimation.correlation_matrix
        )

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, vmin=-1, vmax=1, cmap="coolwarm")
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None

    def absolute_correlation_matrix(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")
        ref_estimation = EstimationResults.from_file(
            self.ref_dir / "estimation.pkl"
        )

        # Get parameters to estimate from configuration
        with (source / "configuration.yaml").open("rb") as buffer:
            raw_config = yaml.safe_load(buffer)
        current_param_config: dict[str, bool] = raw_config["Estimation"][
            "parameters"
        ]

        # Get parameters to estimate in reference
        with (self.ref_dir / "configuration.yaml").open("rb") as buffer:
            ref_raw_config = yaml.safe_load(buffer)
        ref_param_config: dict[str, bool] = ref_raw_config["Estimation"][
            "parameters"
        ]

        # Find largest collection
        param_config = ref_param_config
        if ref_param_config != current_param_config:
            raise ValueError("Cannot compare inconsistent configurations")
        # if len(ref_param_config) >= len(current_param_config):
        #     param_config = ref_param_config
        # else:
        #     param_config = current_param_config

        # Add potentially missing
        potentially_missing = ["gm_deimos"]
        for parameter in potentially_missing:
            if parameter not in param_config:
                param_config[parameter] = False

        # Estimated parameters
        params = get_labels_of_estimated_parameters(param_config)

        # Define title for figure
        canvas_title = (
            f"Correlation difference :: {self.get_id_for_directory(source)} "
            f"vs {self.get_id_for_directory(self.ref_dir)}"
        )
        if len(canvas_title) > 70:
            canvas_title = canvas_title.replace("vs ", "vs\n")
        else:
            print(f"{len(canvas_title)=}")

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="delta-absolute-correlation-matrix.png",
            title=canvas_title,
            canvas_size=(8.5, 3.8),
        )
        # canvas_setup.canvas_title = None

        subfig_setup = ng.PlotSetup(
            custom_ticks=True,
            aspect="equal",
            grid=False,
        )
        abscorr_setup = subfig_setup.version(axtitle=(r"$|\rho_{ref}|$"))
        dabscorr_setup = subfig_setup.version(
            axtitle=(r"$|\rho| - |\rho_{ref}|$")
        )

        # Calculate difference
        abs_rcorr = np.abs(ref_estimation.correlation_matrix)
        abs_ccorr = np.abs(estimation.correlation_matrix)
        quantity = abs_ccorr - abs_rcorr

        with ng.Mosaic("ab", canvas_setup) as canvas:

            with canvas.subplot(abscorr_setup) as fig:

                fig.imshow(abs_rcorr, vmin=0, vmax=1)
                fig.axes["left"].set_xticks(
                    [idx for idx in range(len(params))], labels=params
                )
                fig.axes["left"].set_yticks(
                    [idx for idx in range(len(params))], labels=params
                )

                for i, _ in enumerate(abs_rcorr):
                    for j, value in enumerate(abs_rcorr[i]):
                        fig.axes["left"].text(
                            x=j,
                            y=i,
                            s=f"{value:.1f}",
                            ha="center",
                            va="center",
                            color="k",
                        )  # type: ignore

            with canvas.subplot(dabscorr_setup) as fig:

                fig.imshow(quantity, vmin=-1, vmax=1, cmap="coolwarm")
                fig.axes["left"].set_xticks(
                    [idx for idx in range(len(params))], labels=params
                )
                fig.axes["left"].set_yticks(
                    [idx for idx in range(len(params))], labels=params
                )

                for i, _ in enumerate(quantity):
                    for j, value in enumerate(quantity[i]):
                        fig.axes["left"].text(
                            x=i,
                            y=j,
                            s=f"{value:.1f}",
                            ha="center",
                            va="center",
                            color="k",
                        )  # type: ignore

        return None

    def relative_difference_absolute_correlation(self, source: Path) -> None:

        # Load estimation results
        estimation = EstimationResults.from_file(source / "estimation.pkl")
        ref_estimation = EstimationResults.from_file(
            self.ref_dir / "estimation.pkl"
        )

        # Get parameters to estimate from configuration
        with (source / "configuration.yaml").open("rb") as buffer:
            raw_config = yaml.safe_load(buffer)
        current_param_config: dict[str, bool] = raw_config["Estimation"][
            "parameters"
        ]

        # Get parameters to estimate in reference
        with (self.ref_dir / "configuration.yaml").open("rb") as buffer:
            ref_raw_config = yaml.safe_load(buffer)
        ref_param_config: dict[str, bool] = ref_raw_config["Estimation"][
            "parameters"
        ]

        # Find largest collection
        param_config = ref_param_config
        if ref_param_config != current_param_config:
            raise ValueError("Cannot compare inconsistent configurations")
        # if len(ref_param_config) >= len(current_param_config):
        #     param_config = ref_param_config
        # else:
        #     param_config = current_param_config

        # Add potentially missing
        potentially_missing = ["gm_deimos"]
        for parameter in potentially_missing:
            if parameter not in param_config:
                param_config[parameter] = False

        # Estimated parameters
        params = get_labels_of_estimated_parameters(param_config)

        # Define title for figure
        canvas_title = (
            f"Relative difference in absolute correlation\n"
            f"{self.get_id_for_directory(source)} "
            f"vs {self.get_id_for_directory(self.ref_dir)}"
        )
        # if len(canvas_title) > 70:
        #     canvas_title = canvas_title.replace("vs ", "vs\n")
        # else:
        #     print(f"{len(canvas_title)=}")

        # Canvas setup
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(source),
            figure_name="relative-delta-absolute-correlation-matrix.png",
            title=canvas_title,
            canvas_size=(6, 5),
        )
        canvas_setup.custom_ticks = True
        canvas_setup.aspect = "equal"
        canvas_setup.grid = False

        # Calculate difference
        quantity = 1.0 - (
            np.abs(ref_estimation.correlation_matrix)
            / (np.abs(estimation.correlation_matrix) + 1e-15)
        )

        with ng.SingleAxis(canvas_setup) as fig:

            fig.imshow(quantity, cmap="coolwarm")
            fig.axes["left"].set_xticks(
                [idx for idx in range(len(params))], labels=params
            )
            fig.axes["left"].set_yticks(
                [idx for idx in range(len(params))], labels=params
            )

            for i, _ in enumerate(quantity):
                for j, value in enumerate(quantity[i]):
                    fig.axes["left"].text(
                        x=i,
                        y=j,
                        s=f"{value:.2f}",
                        ha="center",
                        va="center",
                        color="k",
                    )  # type: ignore

        return None


class DistributionVisualizationManager(
    AnalysisManager[ParameterDistributionInput]
):

    def __init__(self, user_input: ParameterDistributionInput) -> None:

        super().__init__(user_input)

    # def postfit_residuals(self, source: Path) -> None:

    #     # Load estimation results
    #     estimation = EstimationResults.from_file(source / "estimation.pkl")
    #     rfig_setup = ng.PlotSetup(
    #         show_axes=False,
    #     )

    #     with ng.Mosaic("aaab") as canvas:

    #         with canvas.subplot() as lfig:
    #             dt = (estimation.epochs - self.ref_epoch) / 3600.0
    #             lfig.line(
    #                 dt, estimation.final_residuals * 1e3, fmt=".", markersize=2
    #             )

    #         with canvas.subplot(rfig_setup) as rfig:

    #             rfig.hist(
    #                 estimation.final_residuals,
    #                 bins=30,
    #                 orientation="horizontal",
    #             )

    #     return None

    def parameter_distribution(
        self, sources: list[Path], requested: int, bins: int
    ) -> None:

        values = []
        for source in sources:

            estimation = EstimationResults.from_file(source / "estimation.pkl")
            parameters = estimation.final_parameters

            # Get names of estimated parameters
            parameter_names = get_parameter_labels_from_source(source)

            # Get index of requested parameter
            values.append(estimation.final_parameters[requested])

        average = np.average(values)
        std = np.std(values)

        statistics = rf"{average:.2e} $\pm$ {std:.2e}"

        with ng.SingleAxis() as fig:

            fig.hist(values, label=statistics, bins=bins)

        return None
