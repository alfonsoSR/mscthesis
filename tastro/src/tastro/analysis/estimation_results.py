from ..io import PropagationOutput, EstimationResults
from ..config import CaseSetup
from nastro import graphics as ng, types as nt
from dataclasses import dataclass
from pathlib import Path
import argparse
import numpy as np
import pickle
from .utils import get_propagation_start_epoch_from_config
from tudatpy.astro import time_representation as ttime
from .prefit import PrefitResults


@dataclass
class ResidualCLI:

    source_dir: Path
    show: bool
    save: bool


class ResidualInputParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "source_dir", help="Source directory with configuration and results"
        )

        self.add_argument("-s", dest="save", action="store_true", help="Save figure")

        self.add_argument(
            "-x", dest="show", action="store_false", help="Do not show figure "
        )

        return None

    def parse_args(self) -> ResidualCLI:

        defaults = super().parse_args()

        source_dir = Path(defaults.source_dir).absolute()
        if not source_dir.exists():
            raise FileExistsError("Invalid source directory")

        return ResidualCLI(
            source_dir=source_dir,
            show=defaults.show,
            save=defaults.save,
        )


def show_residual_history() -> None:

    user_input = ResidualInputParser().parse_args()
    show_residual_history_no_cli(user_input)

    return None


def show_residual_history_no_cli(user_input: ResidualCLI) -> None:

    # user_input = ResidualInputParser().parse_args()
    config = CaseSetup.from_config_file(user_input.source_dir / "configuration.yaml")

    # if (
    #     config.estimation.parameters.drag_coefficient
    #     or config.estimation.parameters.radiation_pressure_coefficient
    # ):
    #     raise NotImplementedError(
    #         "Incompatible with estimation of anything not initial states"
    #     )

    estimation = EstimationResults.from_file(user_input.source_dir / "estimation.pkl")
    propagation = PropagationOutput.from_config_file(
        user_input.source_dir / "configuration.yaml"
    )

    # Ephemerides
    rstate = nt.CartesianState(*propagation.rstate_j2000.T)

    dir_name = user_input.source_dir.relative_to(user_input.source_dir.parents[2])

    match config.propagation.integrator.general.starting_point:

        case "start":
            ref_epoch = propagation.epochs[0]
        case "end":
            ref_epoch = propagation.epochs[-1]
        case "middle":
            ref_epoch = propagation.epochs[len(propagation.epochs) // 2]
        case "custom":
            _ref_epoch = config.propagation.integrator.general.custom_start_epoch
            assert _ref_epoch is not None
            ref_epoch = _ref_epoch.to_float()
        case _:
            raise ValueError("Invalid starting point")

    dt = (propagation.epochs - ref_epoch) / 3600.0

    # Get reference from configuration
    ref = f"Ref: {config.estimation.observations.cartesian.sources[0].path.parent.name}"
    if config.estimation.observations.cartesian.sources[0].use_ephemerides:
        ref += " [Ephemerides]"

    canvas_setup = ng.PlotSetup(
        canvas_size=(12, 6),
        canvas_title=f"Cartesian residuals :: {ref} :: {dir_name}",
        show=user_input.show,
        save=user_input.save,
        dir=user_input.source_dir,
        name="cartesian-residuals.png",
    )

    subfig_setup = ng.PlotSetup(
        xlabel="Hours past start of propagation", legend_location="upper left"
    )

    with ng.Mosaic("ab;cd", canvas_setup) as canvas:

        cfigs = []
        for _label in ("x", "y", "z"):

            _setup = subfig_setup.version(ylabel=(r"$\Delta " + _label + "$ [m]"))
            cfigs.append(canvas.subplot(_setup, generator=ng.DoubleAxis))

        rsetup = subfig_setup.version(ylabel=r"$||\Delta \mathbf{r}||$ [m]")
        rfig = canvas.subplot(rsetup, generator=ng.DoubleAxis)

        # Get length of propagation period
        nepochs = propagation.epochs.shape[0]
        has_cd = config.estimation.parameters.drag_coefficient
        has_cr = config.estimation.parameters.radiation_pressure_coefficient

        print(estimation.final_parameters)

        for idx, residual_set in enumerate(estimation.residual_history.T):

            residuals = np.array(residual_set).reshape(nepochs, 3).T

            residuals_norm = np.linalg.norm(residuals, axis=0)

            for jdx, subfig in enumerate(cfigs):

                subfig.line(dt, residuals[jdx])

            rfig.line(dt, residuals_norm, label=f"Iteration {idx + 1}")

        count = 0
        for subfig in (*cfigs, rfig):

            if count == 3:
                label = "$d_{mars}$ (Right axes)"
            else:
                label = None

            subfig.line(
                dt,
                rstate.r_mag,
                alpha=0.3,
                axis="right",
                color="black",
                fmt="--",
                label=label,
            )

            subfig.__exit__(0, 0, 0)

            count += 1


def show_doppler_residual_history_no_cli(user_input: ResidualCLI) -> None:

    # user_input = ResidualInputParser().parse_args()
    config = CaseSetup.from_config_file(user_input.source_dir / "configuration.yaml")

    estimation = EstimationResults.from_file(user_input.source_dir / "estimation.pkl")
    propagation = PropagationOutput.from_config_file(
        user_input.source_dir / "configuration.yaml"
    )

    # Ephemerides
    rstate = nt.CartesianState(*propagation.rstate_j2000.T)
    with (user_input.source_dir / "epochs.pkl").open("rb") as buffer:

        epochs = np.array(pickle.load(buffer))

    assert isinstance(epochs, np.ndarray)
    assert isinstance(epochs[0], float)

    dir_name = user_input.source_dir.relative_to(user_input.source_dir.parents[2])

    t0 = get_propagation_start_epoch_from_config(config)
    t0_iso = ttime.DateTime.from_epoch_time_object(t0).to_iso_string(
        add_T=True, number_of_digits_seconds=0
    )
    dt = (epochs - t0.to_float()) / 3600.0

    # Load pre-fit residuals
    prefit = PrefitResults.from_file(user_input.source_dir / "prefit_results.pkl")

    # # Get reference from configuration
    # ref = f"Ref: {config.estimation.observations.cartesian.sources[0].path.parent.name}"
    # if config.estimation.observations.cartesian.sources[0].use_ephemerides:
    #     ref += " [Ephemerides]"

    canvas_setup = ng.PlotSetup(
        canvas_size=(8, 6),
        canvas_title=f"Doppler residuals :: {dir_name}",
        show=user_input.show,
        save=user_input.save,
        dir=user_input.source_dir,
        name="doppler-residuals.png",
        xlabel=f"Hours past {t0_iso}",
        ylabel="Frequency residual [mHz]",
    )

    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {t0_iso}",
        ylabel="Residual [mHz]",
    )

    with ng.Mosaic("a;b", canvas_setup) as canvas:

        with canvas.subplot(subfig_setup) as history_fig:

            for idx, residual_set in enumerate(estimation.residual_history.T):

                history_fig.line(
                    dt, residual_set * 1e3, fmt=".", label=f"Iteration {idx}"
                )

        with canvas.subplot(subfig_setup) as comparison_fig:

            comparison_fig.line(
                dt,
                estimation.final_residuals * 1e3,
                fmt=".",
                label="Estimated [Best]",
            )
            comparison_fig.line(
                dt,
                prefit.residual * 1e3,
                fmt=".",
                label="Pre-fit ephemerides",
            )


def show_doppler_residual_history() -> None:

    user_input = ResidualInputParser().parse_args()
    show_doppler_residual_history_no_cli(user_input)

    return None
