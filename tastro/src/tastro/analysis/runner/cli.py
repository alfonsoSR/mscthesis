from ...io.command_line.runner import CommandLineParserRunner
from ...io.command_line.core import CLParser
from .main import SimulationManager
from .prefit import PrefitManager


def runner_cli() -> None:

    user_input = CommandLineParserRunner().parse_args()
    SimulationManager(user_input).run_simulations()

    return None


def calculate_prefit_residuals_cli() -> None:

    user_input = CLParser().parse_args()
    PrefitManager(user_input).calculate_prefit_residuals()

    return None
