from ...io.command_line.core import CLInputFigure
from ..core import AnalysisFigureManagerBase
import typing
from ...io import PropagationOutput
from ..utils import get_propagation_start_epoch_from_config
from ...config import CaseSetup
from tudatpy.astro import time_representation as ttime
from nastro import types as nt, graphics as ng
from ...logging import log
import traceback


class EphemeridesResidualManager(AnalysisFigureManagerBase[CLInputFigure]):

    def __init__(self, user_input: CLInputFigure, frame: str) -> None:

        super().__init__(user_input)

        # Load propagation results
        results = PropagationOutput.from_file(self.source_dir / "results.pkl")

        # Get propagation output in correct frame
        if frame in ("j2000", "rsw"):

            self.cstate = nt.CartesianState[nt.Vector](
                *getattr(results, f"cstate_{frame}").T
            )
            self.rstate = nt.CartesianState[nt.Vector](
                *getattr(results, f"rstate_{frame}").T
            )

        else:

            log.fatal(f"Requested propagation output in invalid frame: {frame}")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Save frame choice
        self.frame = frame

        # Get initial epoch from configuration
        self.epochs = results.epochs

        return None

    @property
    def vector_components(self) -> list[str]:

        match self.frame:

            case "j2000":
                return ["x", r"\dot x", "y", r"\dot y", "z", r"\dot z"]

            case "rsw":
                return ["r", r"\dot r", "s", r"\dot s", "w", r"\dot w"]

            case _:
                raise NotImplementedError("Unreachable")


def propagation_residual_wrt_ephemerides(
    user_input: "CLInputFigure",
    frame: typing.Literal["rsw", "j2000"],
) -> None:

    # Initialize manager
    manager = EphemeridesResidualManager(user_input, frame)

    # Define quantities to plot
    dt = (manager.epochs - manager.ref_epoch) / 3600.0
    residual = manager.cstate - manager.rstate
    dmars = manager.rstate.r_mag * 1e-7

    # Define settings for main figure
    canvas_setup = ng.PlotSetup(
        canvas_size=(8, 6),
        canvas_title=(
            f"Propagation vs ephemerides :: {frame.upper()} "
            f":: {manager.get_directory_id(manager.source_dir)}"
        ),
        dir=manager.source_dir,
        name=f"propagation-vs-ephemerides-{frame}.png",
        save=manager.user_input.save,
        show=manager.user_input.show,
    )

    # Define common settings for subfigure
    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {manager.ref_epoch_isot}",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
        scilimits=(-2, 3),
    )

    # Generate figure
    with ng.Mosaic("ab;cd;ef", canvas_setup) as canvas:

        components = ["x", "dx", "y", "dy", "z", "dz"]

        for idx, label in enumerate(manager.vector_components):

            # Generate settings for current subplot
            unit = "m/s" if ((idx % 2) != 0) else "m"
            current_setup = subfig_setup.version(
                ylabel=rf"$\Delta {label}$ [{unit}]"
            )

            # Generate subfigure
            with canvas.subplot(
                setup=current_setup, generator=ng.DoubleAxis
            ) as subfig:

                subfig.line(dt, getattr(residual, components[idx]))
                subfig.line(dt, dmars, axis="right", alpha=0.5)

    return None
