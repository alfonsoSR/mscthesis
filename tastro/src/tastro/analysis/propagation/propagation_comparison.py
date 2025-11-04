from ...io.command_line.comparison import CLInputComparison, CLParserComparison
from ..core import AnalysisFigureManagerBase
from ...io import PropagationOutput
from nastro import types as nt, graphics as ng
import numpy as np
from ...logging import log


class PropagationComparisonManager(
    AnalysisFigureManagerBase[CLInputComparison]
):

    def __init__(self, user_input: CLInputComparison) -> None:

        super().__init__(user_input)

        # Load results for current and reference configurations
        current = PropagationOutput.from_file(
            self.user_input.source_dir / "results.pkl"
        )
        reference = PropagationOutput.from_file(
            self.user_input.ref_dir / "results.pkl"
        )

        # Extract epochs and states
        self.current = nt.CartesianState[nt.Vector](*current.cstate_j2000.T)
        self.reference = nt.CartesianState[nt.Vector](*reference.cstate_j2000.T)
        self.epochs = current.epochs

        # Ensure dimensions and epochs are the same
        if not np.all(
            np.isclose(self.epochs, reference.epochs, atol=0.0, rtol=1e-15)
        ):
            log.fatal("Cannot compare propagations: epochs do not match")
            exit(1)

        # Define IDs for reference and current directories
        self.current_id = self.get_directory_id(self.user_input.source_dir)
        self.ref_id = self.get_directory_id(self.user_input.ref_dir)

        return None

    @property
    def vector_components(self) -> list[str]:

        return ["x", r"\dot x", "y", r"\dot y", "z", r"\dot z"]


def compare_propagation_results(user_input: CLInputComparison) -> None:

    # Define manager
    manager = PropagationComparisonManager(user_input)

    # Compute quantities to plot
    dt = (manager.epochs - manager.ref_epoch) / 3600.0
    residual = manager.current - manager.reference
    dmars = manager.current.r_mag * 1e-7

    # Define setup for figure
    canvas_setup = ng.PlotSetup(
        canvas_size=manager.default_canvas_size,
        canvas_title=f"Residual between {manager.current_id} and {manager.ref_id} [R]",
        show=manager.user_input.show,
        save=manager.user_input.save,
        dir=manager.user_input.source_dir,
        name=f"propagation-residual-{manager.ref_id.replace("/", "_")}.png",
    )

    # Define common settings for subfigures
    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {manager.ref_epoch_isot}",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
        scilimits=(-2, 3),
    )

    # Generate figure
    with ng.Mosaic("ab;cd;ef", canvas_setup) as canvas:

        components = ["x", "dx", "y", "dy", "z", "dz"]

        for idx, label in enumerate(manager.vector_components):

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
