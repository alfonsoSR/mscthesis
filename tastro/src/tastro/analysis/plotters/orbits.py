from dataclasses import dataclass
from nastro import graphics as ng, types as nt
import numpy as np
from ...io.command_line import ResultVisualizationInput, ResultComparisonInput
from .core import AnalysisManager
from pathlib import Path
from ...io import PropagationOutput, EstimationResults
import enum
from ..utils import get_propagaton_start_epoch_from_config_file
from tudatpy.astro import (
    time_representation as ttime,
    frame_conversion as tframe,
)
from ...logging import log
import traceback
from tudatpy.interface import spice
from tudatpy.math import interpolators as tinterp


class AnalysisType(enum.IntEnum):

    PROPAGATION = 0
    ESTIMATION = 1


class OrbitType(enum.StrEnum):

    SIMULATION = "cstate"
    EPHEMERIDES = "rstate"


def state_config_from_string(config_str: str) -> tuple[AnalysisType, OrbitType]:

    if len(config_str) != 2:
        raise ValueError(f"String must be of size 2: {config_str}")

    match config_str:
        case "pp":
            return (AnalysisType.PROPAGATION, OrbitType.SIMULATION)
        case "pe":
            return (AnalysisType.PROPAGATION, OrbitType.EPHEMERIDES)
        case "ep":
            return (AnalysisType.ESTIMATION, OrbitType.SIMULATION)
        case "ee":
            return (AnalysisType.ESTIMATION, OrbitType.EPHEMERIDES)
        case _:
            raise ValueError(f"Invalid state config string: {config_str}")


def name_component_from_state_config(
    state_config: tuple[AnalysisType, OrbitType],
) -> str:

    match state_config:
        case (AnalysisType.PROPAGATION, OrbitType.SIMULATION):
            return "pp"
        case (AnalysisType.PROPAGATION, OrbitType.EPHEMERIDES):
            return "pe"
        case (AnalysisType.ESTIMATION, OrbitType.SIMULATION):
            return "ep"
        case (AnalysisType.ESTIMATION, OrbitType.EPHEMERIDES):
            return "ee"
        case _:
            raise ValueError(f"Invalid state config: {state_config}")


class ReferenceFrame(enum.StrEnum):

    J2000 = "j2000"
    RSW = "rsw"
    TNW = "tnw"

    @classmethod
    def from_str(cls, frame: str) -> "ReferenceFrame":

        match frame:
            case "j2000":
                return cls.J2000
            case "rsw":
                return cls.RSW
            case "tnw":
                return cls.TNW
            case _:
                raise ValueError("Invalid reference frame")


@dataclass
class FigureData:

    dt: np.ndarray
    state: nt.CartesianState[nt.Vector]
    fmt: str = "-"
    markersize: int = 3
    color: str | None = None
    label: str | None = None

    def __sub__(self, other: "FigureData") -> "FigureData":

        # Fail if epochs are not the same
        if len(self.dt) != len(other.dt):
            print(f"{len(self.dt)=} :: {len(other.dt)=}")
            raise ValueError(
                "Attempted to subtract data with different lengths"
            )
        if not all(np.isclose(self.dt, other.dt, rtol=0, atol=1e-12)):
            raise ValueError("Attempted to subtract data non-coincident epochs")

        values = self.__dict__.copy()
        values["state"] = self.state - other.state
        return FigureData(**values)


# def orbit_from_source_dir(
#     source: Path,
#     ref_epoch: float,
#     analysis_type: AnalysisType,
#     orbit_type: OrbitType,
#     reference_frame: ReferenceFrame,
#     scale: float,
# ) -> FigureData:

#     # Load trajectory
#     match analysis_type:
#         case AnalysisType.PROPAGATION:
#             results = PropagationOutput.from_file(source / "results.pkl")
#         case AnalysisType.ESTIMATION:
#             results = EstimationResults.from_file(
#                 source / "estimation.pkl"
#             ).final_propagation_results
#         case _:
#             raise TypeError(f"Invalid analysis type: {analysis_type}")

#     # Extract cartesian states and epochs
#     state = getattr(results, f"{orbit_type}_{reference_frame}") * scale
#     dt = (results.epochs - ref_epoch) / 3600.0

#     return FigureData(dt, nt.CartesianState(*state.T))


# def state_residual_from_source_dir(
#     source: Path,
#     ref_epoch: float,
#     state_1: tuple[AnalysisType, OrbitType],
#     state_2: tuple[AnalysisType, OrbitType],
#     reference_frame: ReferenceFrame,
# ) -> FigureData:

#     # Load orbits
#     orbit_1 = orbit_from_source_dir(
#         source, ref_epoch, state_1[0], state_1[1], reference_frame
#     )
#     orbit_2 = orbit_from_source_dir(
#         source, ref_epoch, state_2[0], state_2[1], reference_frame
#     )

#     # Ensure that epochs are the same
#     if not np.all(np.isclose(orbit_1.dt, orbit_2.dt, rtol=0, atol=1e-12)):
#         raise ValueError("Epochs should always be the same?")

#     # Calculate state residual and return
#     return orbit_1 - orbit_2


def update_figure[T: ng.BaseFigure](
    fig: T, data: FigureData, component: str, scale
) -> T:

    fig.line(
        data.dt,
        getattr(data.state, component) * scale,
        fmt=data.fmt,
        color=data.color,
        markersize=data.markersize,
        label=data.label,
    )

    return fig


class StateResidualVisualizationInput(ResultVisualizationInput):

    reference_frame: str
    target_state: str
    reference_state: str
    formal_errors: bool = False
    show_variability: bool = True


class StateResidualComparisonInput(ResultComparisonInput):

    reference_frame: str
    target_state: str
    reference_state: str


class StateResidualVisualizationManager(
    AnalysisManager[StateResidualVisualizationInput]
):

    def __init__(self, user_input: StateResidualVisualizationInput) -> None:

        super().__init__(user_input)

        # Extract reference epoch from first configuration file
        ref_epoch = get_propagaton_start_epoch_from_config_file(
            config_file=(self.source_dirs[0] / "configuration.yaml")
        )
        self.ref_epoch = ref_epoch.to_float()
        self.ref_epoch_isot = ttime.DateTime.from_epoch_time_object(
            ref_epoch
        ).to_iso_string(add_T=True, number_of_digits_seconds=0)

        # Define reference frame and configurations
        self.reference_frame = ReferenceFrame.from_str(
            user_input.reference_frame
        )
        self.state_1 = state_config_from_string(user_input.target_state)
        self.state_2 = state_config_from_string(user_input.reference_state)

        # Define wether to plot formal errors or not
        self.formal_errors = user_input.formal_errors
        self.show_variability = user_input.show_variability

        return None

    def __orbit_from_source_dir(
        self,
        source: Path,
        analysis_type: AnalysisType,
        orbit_type: OrbitType,
        reference_frame: ReferenceFrame,
    ) -> FigureData:

        # Load trajectory
        match analysis_type:
            case AnalysisType.PROPAGATION:
                results = PropagationOutput.from_file(source / "results.pkl")
            case AnalysisType.ESTIMATION:
                results = EstimationResults.from_file(
                    source / "estimation.pkl"
                ).final_propagation_results
            case _:
                raise TypeError(f"Invalid analysis type: {analysis_type}")

        # Extract cartesian states and epochs
        state = getattr(results, f"{orbit_type}_{reference_frame}")
        dt = (results.epochs - self.ref_epoch) / 3600.0

        return FigureData(dt, nt.CartesianState(*state.T))

    def __formal_errors_from_source_dir(
        self,
        source: Path,
    ) -> FigureData:

        # Ensure this is an estimation
        if self.state_1[0] != AnalysisType.ESTIMATION:
            log.fatal("Requested fromal errors for propagation analysis")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Load estimation output and ensure covariance is available
        estimation = EstimationResults.from_file(source / "estimation.pkl")
        propagation = estimation.final_propagation_results
        if estimation.covariance_history is None:
            log.fatal(
                "Requested formal errors for estimation without "
                "covariance history"
            )
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Get covariance history in correct reference frame
        match self.reference_frame:

            case ReferenceFrame.J2000:

                covariance_history = np.array(
                    list(estimation.covariance_history.values())
                )

            case ReferenceFrame.RSW:

                covariance_history_j2000 = np.array(
                    list(estimation.covariance_history.values())
                )

                # Compose rotation matrix
                partial_matrix = np.array(
                    [
                        tframe.inertial_to_rsw_rotation_matrix(rstate_i)
                        for rstate_i in propagation.rstate_j2000
                    ]
                )
                rotation_matrix = (
                    np.ones((partial_matrix.shape[0]))[:, None, None]
                    * np.identity(6)[None, :, :]
                )
                rotation_matrix[:, :3, :3] = partial_matrix
                rotation_matrix[:, 3:6, 3:6] = partial_matrix

                # Rotate covariance history
                covariance_history = np.array(
                    [
                        matrix_i @ covariance_i @ matrix_i.T
                        for matrix_i, covariance_i in zip(
                            rotation_matrix, covariance_history_j2000
                        )
                    ]
                )

        # Extract formal errors from covariance matrix
        formal_errors = nt.CartesianState(
            *np.sqrt(np.diagonal(covariance_history, axis1=1, axis2=2).T)
        )

        # covariance_history = (
        #     np.array(list(estimation.covariance_history.values()))
        #     .swapaxes(0, 1)
        #     .swapaxes(1, 2)
        # )
        # _formal_errors = np.array(
        #     [
        #         np.sqrt(covariance_history[idx][idx])
        #         for idx in range(covariance_history.shape[0])
        #     ]
        # )

        # Get epochs at which formal errors are known
        _formal_error_epochs = np.array(
            list(estimation.covariance_history.keys())
        )

        # # Rotate if requested
        # match self.reference_frame:
        #     case ReferenceFrame.J2000:

        #         formal_errors = nt.CartesianState(*_formal_errors)

        #     case ReferenceFrame.RSW:

        #         raise NotImplementedError("RSW formal errors not supported")

        #         # Get state in J2000
        #         rstate_j2000 = estimation.final_propagation_results.rstate_j2000
        #         rstate_history: dict[ttime.Time, np.ndarray] = {
        #             ttime.Time(epoch): state
        #             for epoch, state in zip(
        #                 estimation.final_propagation_results.epochs,
        #                 rstate_j2000,
        #             )
        #         }

        #         # Interpolate state at observation epochs
        #         __interp_settings = tinterp.lagrange_interpolation(8)
        #         rstate_interp = tinterp.create_one_dimensional_vector_interpolator_time_object(
        #             data_to_interpolate=rstate_history,
        #             interpolator_settings=__interp_settings,
        #         )

        #         # Get states in RSW frame (ephemerides)
        #         formal_errors_rsw = np.zeros((_formal_error_epochs.shape[0], 6))
        #         for idx, (epoch, formal_error) in enumerate(
        #             zip(rstate_history, _formal_errors.T)
        #         ):

        #             # # Get state of MEX wrt Mars at epoch
        #             # _rstate = spice.get_body_cartesian_state_at_epoch(
        #             #     "MEX", "Mars", "J2000", "NONE", epoch
        #             # )

        #             # Calculate rotation matrix
        #             rotation_matrix = tframe.inertial_to_rsw_rotation_matrix(
        #                 rstate_interp.interpolate(epoch)
        #             )

        #             # Reference state in RSW
        #             _refpos_rsw = rotation_matrix @ formal_error[:3]
        #             _refvel_rsw = rotation_matrix @ formal_error[3:]
        #             formal_errors_rsw[idx] = np.array(
        #                 [_refpos_rsw, _refvel_rsw]
        #             ).flatten()

        #         # Pack rotated formal errors
        #         formal_errors = nt.CartesianState(*formal_errors_rsw.T)

        #     case _:
        #         log.error("Frame option invalid or not implemented")
        #         log.error(traceback.extract_stack()[-2])
        #         exit(1)

        # Extract cartesian states and epochs
        dt = (_formal_error_epochs - self.ref_epoch) / 3600.0

        return FigureData(dt, formal_errors)

    def __state_residual_from_source_dir(self, source: Path) -> "FigureData":

        # Load orbits
        orbit_1 = self.__orbit_from_source_dir(
            source, self.state_1[0], self.state_1[1], self.reference_frame
        )
        orbit_2 = self.__orbit_from_source_dir(
            source, self.state_2[0], self.state_2[1], self.reference_frame
        )

        # Ensure that epochs are the same
        if not np.all(np.isclose(orbit_1.dt, orbit_2.dt, rtol=0, atol=1e-12)):
            raise ValueError("Epochs should always be the same?")

        # Calculate state residual and return
        return orbit_1 - orbit_2

    def state_residual(self, sources: list[Path]) -> None:

        # Define name of figure from state configurations
        name_1 = name_component_from_state_config(self.state_1)
        name_2 = name_component_from_state_config(self.state_2)
        figure_name = (
            f"state-residual-{self.reference_frame}-{name_1}-{name_2}.png"
        )

        # Define settings for canvas and legend
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name=figure_name,
            title=f"State residual :: {name_1.upper()} vs. {name_2.upper()}",
            canvas_size=(6, 6.5),
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data to plot per source
        data_per_source: dict[Path, FigureData] = {
            source: self.__state_residual_from_source_dir(source)
            for source in sources
        }

        # Get formal errors per source if requested
        if self.formal_errors:
            try:
                spice.load_kernel(str(sources[0] / "metak.tm"))
                formal_errors_per_source: dict[Path, FigureData] = {
                    source: self.__formal_errors_from_source_dir(source)
                    for source in sources
                }
            finally:
                spice.clear_kernels()

        # Generate figure
        mosaic = ("ab;" * 3) + ("cd;" * 3) + ("ef;" * 3) + "gg"
        with ng.Mosaic(mosaic, canvas_setup) as canvas:

            components = ["x", "dx", "y", "dy", "z", "dz"]
            match self.reference_frame:
                case ReferenceFrame.J2000:
                    label_components = components
                case ReferenceFrame.RSW:
                    label_components = ["r", "dr", "s", "ds", "w", "dw"]
                case _:
                    raise ValueError("Frame not implemented")
            labels = [
                rf"\dot {component[1]}" if component[0] == "d" else component[0]
                for component in label_components
            ]
            units = [
                "mm/s" if component[0] == "d" else "m"
                for component in components
            ]
            scales = [
                1e3 if component[0] == "d" else 1 for component in components
            ]

            for idx, component in enumerate(components):

                label = labels[idx]
                unit = units[idx]
                scale = scales[idx]

                # Define setup for subfigure
                subfig_setup = ng.PlotSetup(
                    xlabel=f"Hours past {self.ref_epoch_isot}",
                    ylabel=rf"$\Delta {label}$ [{unit}]",
                    scilimits=(-3, 3),
                )

                # Generate subfigure
                with canvas.subplot(subfig_setup) as subfig:

                    for source in sources:

                        update_figure(
                            subfig, data_per_source[source], component, scale
                        )

                        if self.formal_errors:

                            current_errors = formal_errors_per_source[source]
                            subfig.line(
                                current_errors.dt,
                                current_errors.dt * 0,
                                alpha=0,
                                color="black",
                            )
                            subfig.boundary(
                                getattr(current_errors.state, component)
                                * scale,
                                alpha=0.5,
                            )

            # Add legend
            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )
                    if self.formal_errors:
                        legend.line(-1, -1, fmt="o")

        return None

    def state_formal_errors(self, sources: list[Path]) -> None:

        # Define name of figure from state configurations
        target_name = name_component_from_state_config(self.state_1)
        figure_name = (
            f"state-formal-errors-{self.reference_frame}-{target_name}.png"
        )

        # Define settings for canvas and legend
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name=figure_name,
            title=f"State formal errors :: {target_name}",
            canvas_size=(6, 6.5),
        )
        legend_setup = ng.PlotSetup(legend_columns=self.legend_cols)

        # Get formal errors per source
        formal_errors_per_source: dict[Path, FigureData] = {
            source: self.__formal_errors_from_source_dir(source)
            for source in sources
        }

        # Generate figure
        mosaic = ("ab;" * 3) + ("cd;" * 3) + ("ef;" * 3) + "gg"
        with ng.Mosaic(mosaic, canvas_setup) as canvas:

            components = ["x", "dx", "y", "dy", "z", "dz"]
            match self.reference_frame:
                case ReferenceFrame.J2000:
                    label_components = components
                case ReferenceFrame.RSW:
                    label_components = ["r", "dr", "s", "ds", "w", "dw"]
                case _:
                    raise ValueError("Frame not implemented")
            labels = [
                rf"\dot {component[1]}" if component[0] == "d" else component[0]
                for component in label_components
            ]
            units = [
                "mm/s" if component[0] == "d" else "m"
                for component in components
            ]
            scales = [
                1e3 if component[0] == "d" else 1 for component in components
            ]

            for idx, component in enumerate(components):

                label = labels[idx]
                unit = units[idx]
                scale = scales[idx]

                # Define setup for subfigure
                subfig_setup = ng.PlotSetup(
                    xlabel=f"Hours past {self.ref_epoch_isot}",
                    ylabel=rf"$\Delta {label}$ [{unit}]",
                    scilimits=(-3, 3),
                )

                # Generate subfigure
                with canvas.subplot(subfig_setup) as subfig:

                    for source in sources:
                        update_figure(
                            subfig,
                            formal_errors_per_source[source],
                            component,
                            scale,
                        )

            # Add legend
            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1,
                        -1,
                        fmt="o",
                        label=self.get_id_for_directory(source),
                    )

        return None

    def average_state_residual(self, sources: list[Path]) -> None:

        # Force -g flag on user input
        self.group = True

        # Define name of figure from state configurations
        name_1 = name_component_from_state_config(self.state_1)
        name_2 = name_component_from_state_config(self.state_2)
        figure_name = (
            "average-state-residual"
            f"-{self.reference_frame}-{name_1}-{name_2}.png"
        )

        # Define settings for canvas and legend
        outdir = self.get_output_directory(sources[0])
        canvas_setup = self.generate_default_canvas_setup(
            source=outdir,
            figure_name=figure_name,
            title=(
                f"Average state residual :: {outdir.name} :: "
                f"{name_1.upper()} vs. {name_2.upper()}"
            ),
            canvas_size=(6, 6.5),
        )

        legend_setup = ng.PlotSetup(
            legend_columns=2,
        )

        # Get data to plot per source
        data_per_source: dict[Path, FigureData] = {
            source: self.__state_residual_from_source_dir(source)
            for source in sources
        }

        # Generate figure
        mosaic = ("ab;" * 3) + ("cd;" * 3) + ("ef;" * 3) + "gg"
        with ng.Mosaic(mosaic, canvas_setup) as canvas:

            components = ["x", "dx", "y", "dy", "z", "dz"]
            match self.reference_frame:
                case ReferenceFrame.J2000:
                    label_components = components
                case ReferenceFrame.RSW:
                    label_components = ["r", "dr", "s", "ds", "w", "dw"]
                case _:
                    raise ValueError("Frame not implemented")
            labels = [
                rf"\dot {component[1]}" if component[0] == "d" else component[0]
                for component in label_components
            ]
            units = [
                "mm/s" if component[0] == "d" else "m"
                for component in components
            ]
            scales = [
                1e3 if component[0] == "d" else 1 for component in components
            ]

            for idx, component in enumerate(components):

                label = labels[idx]
                unit = units[idx]
                scale = scales[idx]

                # Define setup for subfigure
                subfig_setup = ng.PlotSetup(
                    xlabel=f"Hours past {self.ref_epoch_isot}",
                    ylabel=rf"$\Delta {label}$ [{unit}]",
                    scilimits=(-3, 3),
                )

                # Generate subfigure
                with canvas.subplot(subfig_setup) as subfig:

                    # Get average of component
                    avg_component = np.average(
                        [
                            getattr(_data.state, component)
                            for _data in data_per_source.values()
                        ],
                        axis=0,
                    )

                    # Get light and dark blue
                    dark = subfig.next_color()
                    light = subfig.next_color()

                    for source in sources:

                        # Set color to light blue
                        data_per_source[source].color = light

                        update_figure(
                            subfig, data_per_source[source], component, scale
                        )

                    subfig.line(
                        data_per_source[source].dt,
                        avg_component * scale,
                        color=dark,
                    )

            # Add legend
            with canvas.subplot(legend_setup, ng.Legend) as legend:

                legend.line(-1, -1, fmt="o", label="Average")
                legend.line(-1, -1, fmt="o", label="Individual solutions")

        return None

    def average_state_residual_magnitude(self, bases: list[Path]) -> None:

        # Force -g flag on user input
        self.group = True

        # Define name of figure from state configurations
        figure_name = "average-state-residual-magnitude.png"

        # Define settings for canvas and legend
        outdir = self.get_output_directory(bases[0])
        canvas_setup = self.generate_default_canvas_setup(
            source=outdir,
            figure_name=figure_name,
            title=(f"Average magnitude of state residuals :: {outdir.name}"),
            canvas_size=(6, 3),
        )

        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get data to plot per source
        data_per_base: dict[Path, dict[Path, FigureData]] = {}
        for base in bases:

            __base_data = {}
            for source in base.iterdir():

                if (not source.is_dir()) or (source.name.startswith(".")):
                    continue

                __base_data[source] = self.__state_residual_from_source_dir(
                    source
                )

            data_per_base[base] = __base_data

        # Generate figure
        components = ["r_mag", "v_mag"]
        scales = [1, 1e3]
        units = ["m", "mm/s"]
        labels = [r"|\Delta \mathbf{r}|", r"|\Delta \mathbf{v}|"]
        with ng.Mosaic("ab;ab;ab;cc", canvas_setup) as canvas:

            for label, unit, scale, component in zip(
                labels, units, scales, components
            ):

                # Define setup for subfigure
                subfig_setup = ng.PlotSetup(
                    xlabel=f"Hours past {self.ref_epoch_isot}",
                    ylabel=rf"${label}$ [{unit}]",
                    scilimits=(-3, 3),
                )

                with canvas.subplot(subfig_setup) as subfig:

                    averages = {}
                    for base, data_per_source in data_per_base.items():

                        # Get light and dark blue
                        dark = subfig.next_color()
                        light = subfig.next_color()

                        # Get average of component
                        avg_component = np.average(
                            [
                                getattr(_data.state, component)
                                for _data in data_per_source.values()
                            ],
                            axis=0,
                        )

                        averages[base] = (dark, avg_component)

                        if not self.show_variability:
                            continue

                        for source in data_per_source:

                            # Set color to light blue
                            data_per_source[source].color = light

                            subfig.line(
                                data_per_source[source].dt,
                                getattr(
                                    data_per_source[source].state, component
                                )
                                * scale,
                                color=light,
                                alpha=0.3,
                            )

                            # update_figure(
                            #     subfig,
                            #     data_per_source[source],
                            #     component,
                            #     scale,
                            # )

                        # subfig.line(
                        #     data_per_source[source].dt,
                        #     avg_component * scale,
                        #     color=dark,
                        # )

                    for base, data_per_source in data_per_base.items():

                        __source = list(data_per_source.keys())[0]
                        if self.show_variability:
                            subfig.line(
                                data_per_source[__source].dt,
                                averages[base][1] * scale,
                                color=averages[base][0],
                            )
                        else:
                            subfig.line(
                                data_per_source[__source].dt,
                                averages[base][1] * scale,
                            )

            with canvas.subplot(legend_setup, generator=ng.Legend) as lfig:

                if self.show_variability:

                    for base in bases:
                        dark = lfig.next_color()
                        light = lfig.next_color()
                        lfig.line(-1, -1, fmt="o", color=dark, label=base.name)

                else:
                    for base in bases:
                        lfig.line(-1, -1, fmt="o", label=base.name)

        # mosaic = ("ab;" * 3) + "cc"
        # with ng.Mosaic(mosaic, canvas_setup) as canvas:

        #     components = ["x", "dx", "y", "dy", "z", "dz"]
        #     match self.reference_frame:
        #         case ReferenceFrame.J2000:
        #             label_components = components
        #         case ReferenceFrame.RSW:
        #             label_components = ["r", "dr", "s", "ds", "w", "dw"]
        #         case _:
        #             raise ValueError("Frame not implemented")
        #     labels = [
        #         rf"\dot {component[1]}" if component[0] == "d" else component[0]
        #         for component in label_components
        #     ]
        #     units = [
        #         "mm/s" if component[0] == "d" else "m"
        #         for component in components
        #     ]
        #     scales = [
        #         1e3 if component[0] == "d" else 1 for component in components
        #     ]

        #     for idx, component in enumerate(components):

        #         label = labels[idx]
        #         unit = units[idx]
        #         scale = scales[idx]

        #         # Define setup for subfigure
        #         subfig_setup = ng.PlotSetup(
        #             xlabel=f"Hours past {self.ref_epoch_isot}",
        #             ylabel=rf"$\Delta {label}$ [{unit}]",
        #             scilimits=(-3, 3),
        #         )

        #         # Generate subfigure
        #         with canvas.subplot(subfig_setup) as subfig:

        #             # Get light and dark blue
        #             dark = subfig.next_color()
        #             light = subfig.next_color()

        #             # Get average of component
        #             avg_component = np.average(
        #                 [
        #                     getattr(_data.state, component)
        #                     for _data in data_per_source.values()
        #                 ],
        #                 axis=0,
        #             )

        #             for source in sources:

        #                 # Set color to light blue
        #                 data_per_source[source].color = light

        #                 update_figure(
        #                     subfig, data_per_source[source], component, scale
        #                 )

        #             subfig.line(
        #                 data_per_source[source].dt,
        #                 avg_component * scale,
        #                 color=dark,
        #             )

        #     # Add legend
        #     with canvas.subplot(legend_setup, ng.Legend) as legend:

        #         legend.line(-1, -1, fmt="o", label="Average")
        #         legend.line(-1, -1, fmt="o", label="Individual solutions")

        return None


class StateResidualComparisonManager(
    AnalysisManager[StateResidualComparisonInput]
):

    def __init__(self, user_input: StateResidualComparisonInput) -> None:

        super().__init__(user_input)

        # Directory with reference solution
        self.ref_dir = user_input.reference

        # Extract reference epoch from first configuration file
        ref_epoch = get_propagaton_start_epoch_from_config_file(
            config_file=(self.ref_dir / "configuration.yaml")
        )
        self.ref_epoch = ref_epoch.to_float()
        self.ref_epoch_isot = ttime.DateTime.from_epoch_time_object(
            ref_epoch
        ).to_iso_string(add_T=True, number_of_digits_seconds=0)

        # Define reference frame and configurations
        self.reference_frame = ReferenceFrame.from_str(
            user_input.reference_frame
        )
        self.target_config = state_config_from_string(user_input.target_state)
        self.ref_config = state_config_from_string(user_input.reference_state)

        return None

    def __orbit_from_source_dir(
        self,
        source: Path,
        analysis_type: AnalysisType,
        orbit_type: OrbitType,
        reference_frame: ReferenceFrame,
    ) -> FigureData:

        # Load trajectory
        match analysis_type:
            case AnalysisType.PROPAGATION:
                results = PropagationOutput.from_file(source / "results.pkl")
            case AnalysisType.ESTIMATION:
                results = EstimationResults.from_file(
                    source / "estimation.pkl"
                ).final_propagation_results
            case _:
                raise TypeError(f"Invalid analysis type: {analysis_type}")

        # Extract cartesian states and epochs
        state = getattr(results, f"{orbit_type}_{reference_frame}")
        dt = (results.epochs - self.ref_epoch) / 3600.0

        return FigureData(dt, nt.CartesianState(*state.T))

    def __state_residual_from_source_dir(self, source: Path) -> "FigureData":

        # Load orbits
        orbit_1 = self.__orbit_from_source_dir(
            source,
            self.target_config[0],
            self.target_config[1],
            self.reference_frame,
        )
        orbit_2 = self.__orbit_from_source_dir(
            source,
            self.ref_config[0],
            self.ref_config[1],
            self.reference_frame,
        )

        # Ensure that epochs are the same
        if not np.all(np.isclose(orbit_1.dt, orbit_2.dt, rtol=0, atol=1e-12)):
            raise ValueError("Epochs should always be the same?")

        # Calculate state residual and return
        return orbit_1 - orbit_2

    def state_residual(self, sources: list[Path]) -> None:

        # raise NotImplementedError(
        #     "Missing handling of different propagation epochs"
        # )

        # Define name of figure from state configurations
        target_id = name_component_from_state_config(self.target_config)
        reference_id = name_component_from_state_config(self.ref_config)
        figure_name = (
            f"state-residual-{self.reference_frame}-"
            f"{target_id}-{reference_id}.png"
        )

        # Define settings for canvas and legend
        canvas_setup = self.generate_default_canvas_setup(
            source=self.get_output_directory(sources[0]),
            figure_name=figure_name,
            title=(
                f"State residual :: {target_id.upper()} "
                f"vs. {reference_id.upper()}"
            ),
            canvas_size=(6, 6.5),
        )
        legend_setup = ng.PlotSetup(
            legend_columns=self.legend_cols,
        )

        # Get state residual for reference case
        reference_residual = self.__state_residual_from_source_dir(self.ref_dir)

        # Get data to plot per source
        residual_difference_per_source: dict[Path, FigureData] = {
            source: (
                self.__state_residual_from_source_dir(source)
                - reference_residual
            )
            for source in sources
        }

        # Generate figure
        mosaic = ("ab;" * 3) + ("cd;" * 3) + ("ef;" * 3) + "gg"
        with ng.Mosaic(mosaic, canvas_setup) as canvas:

            components = ["x", "dx", "y", "dy", "z", "dz"]
            match self.reference_frame:
                case ReferenceFrame.J2000:
                    label_components = components
                case ReferenceFrame.RSW:
                    label_components = ["r", "dr", "s", "ds", "w", "dw"]
                case _:
                    raise ValueError("Frame not implemented")
            labels = [
                rf"\dot {component[1]}" if component[0] == "d" else component[0]
                for component in label_components
            ]
            units = [
                "mm/s" if component[0] == "d" else "m"
                for component in components
            ]
            scales = [
                1e3 if component[0] == "d" else 1 for component in components
            ]

            for idx, component in enumerate(components):

                label = labels[idx]
                unit = units[idx]
                scale = scales[idx]

                # Define setup for subfigure
                subfig_setup = ng.PlotSetup(
                    xlabel=f"Hours past {self.ref_epoch_isot}",
                    ylabel=rf"$\Delta {label}$ [{unit}]",
                    scilimits=(-3, 3),
                )

                # Generate subfigure
                with canvas.subplot(subfig_setup) as subfig:

                    for source in sources:

                        update_figure(
                            subfig,
                            residual_difference_per_source[source],
                            component,
                            scale,
                        )

            # Add legend
            with canvas.subplot(legend_setup, ng.Legend) as legend:

                for source in sources:

                    legend.line(
                        -1, -1, fmt="o", label=self.get_id_for_directory(source)
                    )

        return None


# class OrbitVisualizationManager(VisualizationManager):
