from .residual_ephemerides import propagation_residual_wrt_ephemerides
from ...io.command_line.core import CLParserFigure
from ...io.command_line.comparison import CLParserComparison
from .propagation_comparison import compare_propagation_results


def propagation_residual_wrt_ephemerides_rsw_cli() -> None:

    user_input = CLParserFigure().parse_args()
    propagation_residual_wrt_ephemerides(user_input, "rsw")

    return None


def propagation_residual_wrt_ephemerides_j2000_cli() -> None:

    user_input = CLParserFigure().parse_args()
    propagation_residual_wrt_ephemerides(user_input, "j2000")

    return None


def compare_propagation_results_cli() -> None:

    user_input = CLParserComparison().parse_args()
    compare_propagation_results(user_input)

    return None
