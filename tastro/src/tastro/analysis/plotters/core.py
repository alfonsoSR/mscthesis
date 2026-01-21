from ...io.command_line import ResultComparisonInput, ResultVisualizationInput
from ...io.command_line.plotters import DEFAULT_CANVAS_SIZE, PlotterInput
from pathlib import Path
from ..utils import get_propagaton_start_epoch_from_config_file
from tudatpy.astro import time_representation as ttime
from nastro import graphics as ng
import multiprocessing as mp


class AnalysisManager[T: PlotterInput]:

    def __init__(self, user_input: T) -> None:

        # Extract source directories
        self.source_dirs: list[Path] = user_input.source_dirs
        self._single_source = len(self.source_dirs) == 1

        # User-defined canvas settings
        self.canvas_size = user_input.canvas_size
        self.legend_cols = user_input.legend_columns
        self.legend_location = user_input.legend_location

        # Figure processing options
        self.save = user_input.save
        self.show = user_input.show
        self.group = user_input.group
        self.name_modifier = user_input.name_modifier
        self.__id_refdir = user_input.refdir
        self.outdir = user_input.outdir

        return None

    def concurrent_execution(self, function_name: str) -> None:

        if len(self.source_dirs) > 1:
            with mp.Pool(6) as pool:
                pool.map(getattr(self, function_name), self.source_dirs)
        else:
            getattr(self, function_name)(self.source_dirs[0])

        return None

    def execute_function(self, function_name: str) -> None:

        if (not self._single_source) and self.group:
            getattr(self, function_name)(self.source_dirs)
        else:
            for source_dir in self.source_dirs:
                getattr(self, function_name)([source_dir])

        return None

    def get_canvas_size(
        self, custom: tuple[float, float]
    ) -> tuple[float, float]:

        if self.canvas_size != DEFAULT_CANVAS_SIZE:
            return self.canvas_size
        else:
            return custom

    def get_id_for_directory(self, directory: Path) -> str:

        return str(directory.relative_to(self.__id_refdir))

    def get_output_directory(self, source: Path) -> Path:

        # If group requested and more than one source
        if (not self._single_source) and self.group:
            return self.outdir

        # Otherwise, save in source directory
        return source

    def generate_default_canvas_setup(
        self,
        source: Path,
        figure_name: str,
        title: str,
        canvas_size: tuple[float, float] = DEFAULT_CANVAS_SIZE,
    ) -> ng.PlotSetup:

        # Get output directory
        outdir = self.get_output_directory(source)

        # Update figure name with modifier
        if self.name_modifier != "":
            __name = Path(figure_name)
            figure_name = f"{__name.stem}-{self.name_modifier}{__name.suffix}"

        return ng.PlotSetup(
            canvas_title=title,
            canvas_size=self.get_canvas_size(canvas_size),
            show=self.show,
            save=self.save,
            dir=outdir,
            name=figure_name,
            legend_columns=self.legend_cols,
            legend_location=self.legend_location,
        )


class VisualizationManager(AnalysisManager[ResultVisualizationInput]):

    def __init__(self, user_input: ResultVisualizationInput) -> None:

        super().__init__(user_input)

        # Extract reference epoch from first configuration file
        ref_epoch = get_propagaton_start_epoch_from_config_file(
            config_file=(self.source_dirs[0] / "configuration.yaml")
        )
        self.ref_epoch = ref_epoch.to_float()
        self.ref_epoch_isot = ttime.DateTime.from_epoch_time_object(
            ref_epoch
        ).to_iso_string(add_T=True, number_of_digits_seconds=0)

        return None


class ComparisonManager(AnalysisManager[ResultComparisonInput]):

    def __init__(self, user_input: ResultComparisonInput) -> None:

        super().__init__(user_input)

        # Directory with reference solution
        self.ref_dir = user_input.reference

        # Extract reference epoch from first configuration file
        ref_epoch = get_propagaton_start_epoch_from_config_file(
            config_file=(self.ref_dir / "configuration.yaml")
        )
        self.ref_epoch = ref_epoch.to_float()
        self.ref_epoch_isot = ttime.DateTime.from_epoch_time_object(
            ref_epoch
        ).to_iso_string(add_T=True, number_of_digits_seconds=0)

        return None
