from ...io.command_line.runner import CommandLineParserRunner
from ...io.command_line.core import CLParser
from .main import SimulationManager
from .prefit import calculate_prefit_residuals


def runner_cli() -> None:

    user_input = CommandLineParserRunner().parse_args()
    SimulationManager(user_input).run_simulations()

    return None


def calculate_prefit_residuals_cli() -> None:

    user_input = CLParser().parse_args()
    calculate_prefit_residuals(user_input)

    return None
