from typing import Any
from ...io.command_line import (
    ResultVisualizationParser,
    ResultComparisonParser,
    ParameterDistributionInput,
)
from ...io.command_line.plotters import CommandLinePlotterParser
from .residuals import DopplerVisualizationManager, DopplerComparisonManager
from .orbits import (
    StateResidualVisualizationManager,
    StateResidualVisualizationInput,
    StateResidualComparisonInput,
    StateResidualComparisonManager,
)
from .statistics import (
    MatrixVisualizationManager,
    MatrixComparisonManager,
    DistributionVisualizationManager,
)
from pathlib import Path


class StateResidualVisualizationParser(
    CommandLinePlotterParser[StateResidualVisualizationInput]
):

    namespace = StateResidualVisualizationInput

    def __init__(self) -> None:

        super().__init__()

        target_orbits = self.add_argument_group("Target orbits")
        target_orbits.add_argument(
            "-f",
            "--frame",
            dest="reference_frame",
            default="rsw",
            choices=["rsw", "tnw", "j2000"],
            help="Reference frame in which orbits are described",
        )
        target_orbits.add_argument(
            "-r",
            "--reference-state",
            dest="reference_state",
            choices=["pp", "pe", "ep", "ee"],
            help="Orbit to use as reference",
            required=True,
        )
        target_orbits.add_argument(
            "-t",
            "--target-state",
            dest="target_state",
            choices=["pp", "pe", "ep", "ee"],
            help="Orbit to use as target",
            required=True,
        )
        target_orbits.add_argument(
            "--formal-errors",
            dest="formal_errors",
            help="Whether to show formal errors",
            action="store_true",
        )

        self.add_argument(
            "--show-variability",
            dest="show_variability",
            action="store_true",
        )

        return None

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)

        # Add reference frame and state configurations
        for argument in [
            "reference_frame",
            "reference_state",
            "target_state",
            "formal_errors",
            "show_variability",
        ]:
            arguments[argument] = getattr(defaults, argument)

        return arguments


class StateResidualComparisonParser(
    CommandLinePlotterParser[StateResidualComparisonInput]
):

    namespace = StateResidualComparisonInput

    def __init__(self) -> None:

        super().__init__()

        # Source of reference results
        reference = self.add_argument_group("Reference results")
        reference.add_argument(
            "-r",
            dest="reference",
            required=True,
            help="Source directory with reference results",
        )

        target_orbits = self.add_argument_group("Target orbits")
        target_orbits.add_argument(
            "-f",
            "--frame",
            dest="reference_frame",
            default="rsw",
            choices=["rsw", "tnw", "j2000"],
            help="Reference frame in which orbits are described",
        )
        target_orbits.add_argument(
            "--reference-state",
            dest="reference_state",
            choices=["pp", "pe", "ep", "ee"],
            help="Orbit to use as reference",
            required=True,
        )
        target_orbits.add_argument(
            "--target-state",
            dest="target_state",
            choices=["pp", "pe", "ep", "ee"],
            help="Orbit to use as target",
            required=True,
        )

        return None

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)

        # Parse directory with reference results
        reference = Path(defaults.reference).absolute()
        self._ensure_path_exists(reference, "directory with reference results")
        arguments["reference"] = reference

        # Add reference frame and state configurations
        for argument in [
            "reference_frame",
            "reference_state",
            "target_state",
        ]:
            arguments[argument] = getattr(defaults, argument)

        return arguments


def doppler_prefits() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "pre_fit_residuals"
    )

    return None


def doppler_prefits_distribution() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "pre_fit_residuals_with_distribution"
    )

    return None


def doppler_postfits() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "post_fit_residuals"
    )

    return None


def doppler_postfits_history() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "post_fit_residual_history"
    )

    return None


def doppler_postfits_distribution() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "post_fit_residuals_with_distribution"
    )

    return None


def doppler_residuals() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function("residuals")

    return None


def doppler_residuals_distribution() -> None:

    user_input = ResultVisualizationParser().parse_args()
    DopplerVisualizationManager(user_input).execute_function(
        "residuals_with_distribution"
    )

    return None


def delta_doppler_prefits() -> None:

    user_input = ResultComparisonParser().parse_args()
    DopplerComparisonManager(user_input).execute_function("pre_fit_residuals")

    return None


def delta_doppler_postfits() -> None:

    user_input = ResultComparisonParser().parse_args()
    DopplerComparisonManager(user_input).execute_function("post_fit_residuals")

    return None


def delta_doppler_residuals() -> None:

    user_input = ResultComparisonParser().parse_args()
    DopplerComparisonManager(user_input).execute_function("residuals")

    return None


def state_residuals() -> None:

    user_input = StateResidualVisualizationParser().parse_args()
    StateResidualVisualizationManager(user_input).execute_function(
        "state_residual"
    )

    return None


def state_formal_errors() -> None:

    user_input = StateResidualVisualizationParser().parse_args()
    StateResidualVisualizationManager(user_input).execute_function(
        "state_formal_errors"
    )

    return None


def average_state_residual() -> None:

    user_input = StateResidualVisualizationParser().parse_args()
    StateResidualVisualizationManager(user_input).average_state_residual(
        user_input.source_dirs
    )

    return None


def average_state_residual_magnitude() -> None:

    user_input = StateResidualVisualizationParser().parse_args()
    StateResidualVisualizationManager(
        user_input
    ).average_state_residual_magnitude(user_input.source_dirs)

    return None


def delta_state_residuals() -> None:

    user_input = StateResidualComparisonParser().parse_args()
    StateResidualComparisonManager(user_input).execute_function(
        "state_residual"
    )

    return None


def correlation_matrix() -> None:

    user_input = ResultVisualizationParser().parse_args()
    MatrixVisualizationManager(user_input).concurrent_execution(
        "correlation_matrix"
    )


def covariance_matrix() -> None:

    user_input = ResultVisualizationParser().parse_args()
    MatrixVisualizationManager(user_input).concurrent_execution(
        "covariance_matrix"
    )


def absolute_correlation_matrix() -> None:

    user_input = ResultVisualizationParser().parse_args()
    MatrixVisualizationManager(user_input).concurrent_execution(
        "absolute_correlation_matrix"
    )


def anticorrelation_matrix() -> None:

    user_input = ResultVisualizationParser().parse_args()
    MatrixVisualizationManager(user_input).concurrent_execution(
        "anticorrelation_matrix"
    )


def delta_correlation_matrix() -> None:

    user_input = ResultComparisonParser().parse_args()
    MatrixComparisonManager(user_input).concurrent_execution(
        "correlation_matrix"
    )


def delta_absolute_correlation_matrix() -> None:

    user_input = ResultComparisonParser().parse_args()
    MatrixComparisonManager(user_input).concurrent_execution(
        "absolute_correlation_matrix"
    )


def relative_difference_absolute_correlation() -> None:

    user_input = ResultComparisonParser().parse_args()
    MatrixComparisonManager(user_input).concurrent_execution(
        "relative_difference_absolute_correlation"
    )

    return None


class ParameterDistributionParser(
    CommandLinePlotterParser[ParameterDistributionInput]
):

    namespace = ParameterDistributionInput

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "--requested",
            dest="requested",
            required=True,
            help="Index of the requested parameter",
        )
        self.add_argument(
            "-b",
            dest="bins",
            help="Number of bins",
            default=5,
        )

        return None

    def local_parser(
        self, defaults, arguments: dict[str, Any]
    ) -> dict[str, Any]:

        arguments = super().local_parser(defaults, arguments)
        arguments["requested"] = int(defaults.requested)
        arguments["bins"] = int(defaults.bins)

        return arguments


def parameter_distribution() -> None:

    user_input = ParameterDistributionParser().parse_args()
    DistributionVisualizationManager(user_input).parameter_distribution(
        user_input.source_dirs, user_input.requested, user_input.bins
    )

    return None
