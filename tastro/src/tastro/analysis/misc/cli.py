from .integrator_analysis import FixedIntegratorAnalysisManager
from ...io.command_line.core import CLParserFigure


def compare_fixed_step_integrators_cli() -> None:

    user_input = CLParserFigure().parse_args()
    manager = FixedIntegratorAnalysisManager(user_input)

    for integrator in manager.results:

        manager.generate_integrator_figure(integrator)

    return None
