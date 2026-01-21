from .core import (
    propagation_from_source_directory,
    estimation_from_source_directory,
    prefit_analysis_from_source_directory,
)
from ..core import ExecutionManager
from ...io.command_line import CLParser
from ...logging import log, FileHandler, StdoutHandler


def tpropagator_cli() -> None:

    user_input = CLParser().parse_args()
    ExecutionManager(user_input).concurrent_execution(
        propagation_from_source_directory
    )

    return None


def testimator_cli() -> None:

    user_input = CLParser().parse_args()
    ExecutionManager(user_input).concurrent_execution(
        estimation_from_source_directory
    )

    return None


def tprefit_cli() -> None:

    user_input = CLParser().parse_args()
    ExecutionManager(user_input).concurrent_execution(
        prefit_analysis_from_source_directory
    )

    return None
