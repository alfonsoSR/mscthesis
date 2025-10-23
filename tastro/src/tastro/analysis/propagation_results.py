from ..io import PropagationOutput, EstimationResults, UserInputParser
from nastro import graphics as ng, types as nt
from dataclasses import dataclass
from pathlib import Path
import argparse
import numpy as np
from tudatpy.astro import time_representation as ttime


@dataclass
class ComparisonCommandLineArguments:

    reference: Path
    current: Path
    use_ephemerides: bool


class ComparisonInputParser(argparse.ArgumentParser):

    def __init__(self) -> None:

        super().__init__()

        self.add_argument(
            "-r",
            dest="reference",
            required=True,
            help="Name of folder with reference results",
        )
        self.add_argument(
            "-c",
            dest="current",
            required=True,
            help="Name of folder with current results",
        )
        self.add_argument(
            "-e",
            dest="use_ephemerides",
            action="store_true",
            help="Compare to ephemerides of reference directory",
        )

        self.add_argument("-b", "--base", dest="base_dir", default=".")

        return None

    def parse_args(self) -> ComparisonCommandLineArguments:

        # Get default output parse_args
        defaults = super().parse_args()

        # Define path to base directory
        base_dir = Path(defaults.base_dir).absolute()
        if not base_dir.exists():
            raise FileNotFoundError(f"Invalid base directory: {base_dir}")
        if not base_dir.is_dir():
            raise ValueError(f"Not a directory: {base_dir}")

        # Ensure reference and current directories exist
        current_dir: Path = base_dir / defaults.current
        if not current_dir.exists():
            raise FileNotFoundError(f"Invalid current directory: {current_dir}")
        reference_dir: Path = base_dir / defaults.reference
        if not reference_dir.exists():
            raise FileNotFoundError(
                f"Invalid reference directory: {reference_dir}"
            )

        # Package output in dataclass
        return ComparisonCommandLineArguments(
            reference=reference_dir,
            current=current_dir,
            use_ephemerides=defaults.use_ephemerides,
        )


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

        self.add_argument(
            "-s", dest="save", action="store_true", help="Save figure"
        )

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


def compare_propagations() -> None:

    comparison_input = ComparisonInputParser().parse_args()

    current = PropagationOutput.from_config_file(
        comparison_input.current / "configuration.yaml"
    )
    reference = PropagationOutput.from_config_file(
        comparison_input.reference / "configuration.yaml"
    )

    epochs = current.epochs
    c_cstate = nt.CartesianState(*current.cstate_j2000.T)
    c_rstate = nt.CartesianState(*current.rstate_j2000.T)
    r_cstate = nt.CartesianState(*reference.cstate_j2000.T)
    r_rstate = nt.CartesianState(*reference.rstate_j2000.T)

    with ng.CompareCartesianStates() as fig:

        if comparison_input.use_ephemerides:
            fig.compare_states(epochs, c_cstate, r_rstate)
        else:
            fig.compare_states(epochs, c_cstate, r_cstate)

    return None


def compare_to_ephemerides_rsw() -> None:

    user_input = ResidualInputParser().parse_args()
    compare_to_ephemerides_rsw_no_cli(user_input)

    return None


def compare_to_ephemerides_rsw_no_cli(user_input: ResidualCLI) -> None:

    # user_input = ResidualInputParser().parse_args()
    results = PropagationOutput.from_file(user_input.source_dir / "results.pkl")

    cstate = nt.CartesianState[nt.Vector](*results.cstate_rsw.T)
    rstate = nt.CartesianState[nt.Vector](*results.rstate_rsw.T)
    dt = (results.epochs - results.epochs[0]) / 3600.0

    # Get ISO-T representation of initial epoch
    initial_epoch_isot = ttime.DateTime.from_epoch(
        results.epochs[0]
    ).to_iso_string(add_T=True, number_of_digits_seconds=0)

    # Get relative path to configuration directory
    dir_relpath = user_input.source_dir.relative_to(
        user_input.source_dir.parents[2]
    )

    fig_setup = ng.PlotSetup(
        canvas_size=(12, 6),
        canvas_title=f"Propagation vs ephemerides: Mars-RSW :: {dir_relpath}",
        # xlabel=f"Hours past {initial_epoch_isot}",
        dir=user_input.source_dir,
        name="propagation-vs-ephemerides-rsw.png",
        save=user_input.save,
        show=user_input.show,
    )

    diff = cstate - rstate
    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {initial_epoch_isot}",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
        scilimits=(-2, 3),
    )

    with ng.Mosaic("ab;cd;ef", fig_setup) as canvas:

        # Delta r
        dr_setup = subfig_setup.version(ylabel=r"$\Delta r$ [m]")
        with canvas.subplot(setup=dr_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.x)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta dr
        ddr_setup = subfig_setup.version(ylabel=r"$\Delta \dot r$ [m/s]")
        with canvas.subplot(setup=ddr_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dx)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta s
        ds_setup = subfig_setup.version(ylabel=r"$\Delta s$ [m]")
        with canvas.subplot(setup=ds_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.y)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta ds
        dds_setup = subfig_setup.version(ylabel=r"$\Delta \dot s$ [m/s]")
        with canvas.subplot(setup=dds_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dy)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta w
        dw_setup = subfig_setup.version(ylabel=r"$\Delta w$ [m]")
        with canvas.subplot(setup=dw_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.z)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta dw
        ddw_setup = subfig_setup.version(ylabel=r"$\Delta \dot w$ [m/s]")
        with canvas.subplot(setup=ddw_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dz)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

    # with ng.CompareRswStates(fig_setup) as fig:

    #     fig.compare_states(dt, cstate, rstate, is_dt=True)

    #     # Add distances to periapsis
    #     __subplots = [getattr(fig, f"q{i}_subplot") for i in range(1, 7)]
    #     for __subplot in __subplots:
    #         __subplot.line(dt, rstate.r_mag, )

    return None


def compare_to_ephemerides_j2000() -> None:

    user_input = ResidualInputParser().parse_args()
    compare_to_ephemerides_j2000_no_cli(user_input)

    return None


def compare_to_ephemerides_j2000_no_cli(user_input: ResidualCLI) -> None:

    # user_input = ResidualInputParser().parse_args()
    results = PropagationOutput.from_file(user_input.source_dir / "results.pkl")

    cstate = nt.CartesianState[nt.Vector](*results.cstate_j2000.T)
    rstate = nt.CartesianState[nt.Vector](*results.rstate_j2000.T)
    dt = (results.epochs - results.epochs[0]) / 3600.0

    # Get ISO-T representation of initial epoch
    initial_epoch_isot = ttime.DateTime.from_epoch(
        results.epochs[0]
    ).to_iso_string(add_T=True, number_of_digits_seconds=0)

    # Get relative path to configuration directory
    dir_relpath = user_input.source_dir.relative_to(
        user_input.source_dir.parents[2]
    )

    fig_setup = ng.PlotSetup(
        canvas_size=(12, 6),
        canvas_title=f"Propagation vs ephemerides: Mars-J2000 :: {dir_relpath}",
        # xlabel=f"Hours past {initial_epoch_isot}",
        dir=user_input.source_dir,
        name="propagation-vs-ephemerides-j2000.png",
        save=user_input.save,
        show=user_input.show,
    )

    diff = cstate - rstate
    subfig_setup = ng.PlotSetup(
        xlabel=f"Hours past {initial_epoch_isot}",
        rlabel=r"$d_{mars}$ [$x10^{-7}$ m]",
        scilimits=(-2, 3),
    )

    with ng.Mosaic("ab;cd;ef", fig_setup) as canvas:

        # Delta r
        dr_setup = subfig_setup.version(ylabel=r"$\Delta x$ [m]")
        with canvas.subplot(setup=dr_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.x)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta dr
        ddr_setup = subfig_setup.version(ylabel=r"$\Delta \dot x$ [m/s]")
        with canvas.subplot(setup=ddr_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dx)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta s
        ds_setup = subfig_setup.version(ylabel=r"$\Delta y$ [m]")
        with canvas.subplot(setup=ds_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.y)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta ds
        dds_setup = subfig_setup.version(ylabel=r"$\Delta \dot y$ [m/s]")
        with canvas.subplot(setup=dds_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dy)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta w
        dw_setup = subfig_setup.version(ylabel=r"$\Delta z$ [m]")
        with canvas.subplot(setup=dw_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.z)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

        # Delta dw
        ddw_setup = subfig_setup.version(ylabel=r"$\Delta \dot z$ [m/s]")
        with canvas.subplot(setup=ddw_setup, generator=ng.DoubleAxis) as subfig:
            subfig.line(dt, diff.dz)
            subfig.line(dt, rstate.r_mag * 1e-7, axis="right", alpha=0.5)

    return None
