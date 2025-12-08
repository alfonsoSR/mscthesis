from .cartesian_residuals import show_cartesian_residuals
from .doppler_residuals import show_doppler_residuals
from ...io.command_line.core import CLParserFigure, CLParser
from .prefit import show_doppler_prefit_residuals
from .doppler_residuals import (
    ResidualHistoryManager,
    show_estimation_summary_single,
    show_correlation_matrix_single,
)
from .orbit_comparison import OrbitComparisonManager
from .doppler import DopplerEstimationManager
from ..core import ExecutionManager


def show_cartesian_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    ResidualHistoryManager(user_input, 3).show_residual_history()
    # show_cartesian_residuals(user_input)

    return None


def show_doppler_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    DopplerEstimationManager(user_input).show_residual_history()
    # ResidualHistoryManager(user_input, 1).show_residual_history()
    # show_doppler_residuals(user_input)

    return None


def show_doppler_estimation_results_cli() -> None:

    user_input = CLParserFigure().parse_args()
    DopplerEstimationManager(user_input).show_complete_estimation_results()

    return None


def show_doppler_prefit_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    show_doppler_prefit_residuals(user_input)

    return None


def show_estimation_summary_cli() -> None:

    user_input = CLParser().parse_args()
    ExecutionManager(user_input).concurrent_execution(
        show_estimation_summary_single
    )

    return None


def show_correlation_matrix_cli() -> None:

    user_input = CLParser().parse_args()
    ExecutionManager(user_input).concurrent_execution(
        show_correlation_matrix_single
    )

    return None


def estimated_orbit_vs_ephemerides_rsw_cli() -> None:

    user_input = CLParserFigure().parse_args()
    OrbitComparisonManager(
        user_input, "rsw"
    ).propagation_residual_wrt_ephemerides()

    return None


def estimated_orbit_vs_ephemerides_j2000_cli() -> None:

    user_input = CLParserFigure().parse_args()
    OrbitComparisonManager(
        user_input, "j2000"
    ).propagation_residual_wrt_ephemerides()

    return None
