from nastro import graphics as ng, types as nt
import numpy as np
from proptools import config as pconfig, io as pio
from pathlib import Path
import argparse
from astropy import time

parser = argparse.ArgumentParser()
parser.add_argument("config_file", help="Path to configuration file")


def get_plot_version(filename_base: str, outdir: Path) -> str:

    # Get version of plot
    previous_versions = [xi for xi in outdir.glob(f"*{filename_base}*.png")]
    if len(previous_versions) != 0:
        cvrsn = f"v{len(previous_versions)}"
    else:
        cvrsn = "v0"

    return cvrsn


if __name__ == "__main__":

    # Get path to configuration file
    config_path = Path(parser.parse_args().config_file).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    config = pconfig.PropSettings(config_path)

    # Load results from file
    results = pio.PropagationOutput.from_file(
        config_path.parent / "results.pkl"
    )

    # Cartesian states in Mars-centered J2000
    cstate_j2000 = nt.CartesianState[nt.Double](*results.cstate_j2000.T)
    rstate_j2000 = nt.CartesianState[nt.Double](*results.rstate_j2000.T)

    # Position difference in RSW
    cstate_rsw = nt.CartesianState[nt.Double](*results.cstate_rsw.T)
    rstate_rsw = nt.CartesianState[nt.Double](*results.rstate_rsw.T)
    diff = cstate_rsw - rstate_rsw

    # ISO string representation of the initial epoch
    j2000_tdb = time.Time("2000-01-01T12:00:00", scale="tt").tdb
    initial_epoch = j2000_tdb + time.TimeDelta(
        results.epochs[0], scale="tdb", format="sec"
    )
    isot_t0: str = initial_epoch.isot  # type: ignore

    # Figure setup: RSW error
    if config.plots.rsw_error:

        filename_base = "rsw-comparison"
        version = get_plot_version(filename_base, config_path.parent)
        canvas_setup = ng.PlotSetup(
            canvas_size=(12, 7),
            canvas_title=f"Difference in RSW state: Propagated vs ephemerides",
            show=config.plots.show,
            save=config.plots.save,
            dir=config_path.parent,
            name=f"rsw-comparison-{version}.png",
            xlabel=f"Days past {isot_t0} [TDB]",
        )
        with ng.CompareRswStates(canvas_setup) as fig:

            fig.compare_states(results.epochs, cstate_rsw, rstate_rsw)

    # 3D orbit
    if config.plots.orbit:

        with ng.PlotOrbit() as fig:

            fig.add_orbit(cstate_j2000, label="Current")
            fig.add_orbit(rstate_j2000, label="Reference")
