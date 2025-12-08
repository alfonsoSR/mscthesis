from ..logging import log
from ..io.command_line.core import CommandLineInput, CLInputFigure
from .utils import get_propagation_start_epoch_from_config
from ..config import CaseSetup
from tudatpy.astro import time_representation as ttime
from pathlib import Path
from nastro import graphics as ng
from typing import Callable
import multiprocessing as mp


class ExecutionManager[T: "CommandLineInput"]:
    """Manager for the concurrent execution of analyses"""

    def __init__(self, user_input: T) -> None:

        # Save source directories
        self.source_dirs: list[Path] = user_input.source_dirs

        return None

    def concurrent_execution(self, function: Callable) -> None:

        if len(self.source_dirs) > 1:
            with mp.Pool() as pool:
                pool.map(function, self.source_dirs)
        else:
            function(self.source_dirs[0])

        return None


class AnalysisManagerBase[T: "CommandLineInput"]:

    def __init__(self, user_input: T) -> None:

        # Save user input and source directory as attributes
        self.user_input = user_input
        self.source_dir = user_input.source_dir

        # Adjust verbosity based on user input
        if not self.user_input.verbose:
            log.setLevel("INFO")

        return None


class AnalysisFigureManagerBase[T: "CLInputFigure"](AnalysisManagerBase[T]):

    def __init__(self, user_input: T) -> None:

        super().__init__(user_input)

        # Load configuration in source directory
        log.setLevel("ERROR")
        self.config = CaseSetup.from_config_file(
            self.source_dir / "configuration.yaml"
        )
        log.setLevel("DEBUG")

        # Define reference propagation epoch as float and ISO string
        __epoch = get_propagation_start_epoch_from_config(self.config)
        self.ref_epoch = __epoch.to_float()
        self.ref_epoch_isot = ttime.DateTime.from_epoch_time_object(
            __epoch
        ).to_iso_string(add_T=True, number_of_digits_seconds=0)

        # Define defaults for configuration of figures
        self.default_canvas_size = (8, 6)

        # Define iterator for source directories
        self.sources_iterator = iter(self.user_input.source_dirs)

        return None

    def get_directory_id(self, directory: Path) -> str:

        return str(directory.relative_to(self.user_input.base_dir))

    def generate_canvas_setup(
        self, name: str, save_in_base: bool = True
    ) -> ng.PlotSetup:

        # Adapt output directory based on number of sources
        if (len(self.user_input.source_dirs) > 1) and save_in_base:
            outdir = self.user_input.base_dir
        else:
            outdir = next(self.sources_iterator)

        # Update name with modifier
        if self.user_input.name_modifier != "":
            __name_path = Path(name)
            name = (
                f"{__name_path.stem}-{self.user_input.name_modifier}"
                f"{__name_path.suffix}"
            )

        return ng.PlotSetup(
            # canvas_size=self.default_canvas_size,
            canvas_size=(6, 6),
            show=self.user_input.show,
            save=self.user_input.save,
            dir=outdir,
            name=name,
        )

    def generate_figure_setup(
        self, ylabel: str, right_axis: bool = True
    ) -> ng.PlotSetup:

        return ng.PlotSetup(
            xlabel=f"Hours past {self.ref_epoch_isot}",
            ylabel=ylabel,
            rlabel=r"d$_{mars}$" if right_axis else None,
            scilimits=(-2, 3),
            show_tick_labels_r=False if right_axis else True,
        )
