from .cartesian_residuals import show_cartesian_residuals
from .doppler_residuals import show_doppler_residuals
from ...io.command_line.core import CLParserFigure
from .prefit import show_doppler_prefit_residuals
from .doppler_residuals import ResidualHistoryManager


def show_cartesian_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    ResidualHistoryManager(user_input, 3).show_residual_history()
    # show_cartesian_residuals(user_input)

    return None


def show_doppler_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    ResidualHistoryManager(user_input, 1).show_residual_history()
    # show_doppler_residuals(user_input)

    return None


def show_doppler_prefit_residuals_cli() -> None:

    user_input = CLParserFigure().parse_args()
    show_doppler_prefit_residuals(user_input)

    return None
