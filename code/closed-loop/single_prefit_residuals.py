from nastro import graphics as ng
from pathlib import Path
import csv
import numpy as np
from astropy import time
from prefit import paths as ppaths

source_dir = ppaths.outdir / "prefit-closed-loop"
figdir = ppaths.figdir / "prefit-closed-loop"
figdir.mkdir(exist_ok=True)
LETTERS = iter([x for x in "abcdefghijklmnopqrstuvwxyz"])


def get_comparison_version(savefig_flag: bool) -> int:

    cache_file = Path(__file__).parent / "comparison_version.txt"
    with cache_file.open() as buffer:
        comparison_version = int(buffer.read(1))

    if not savefig_flag:
        return comparison_version

    # Update version
    with cache_file.open("w") as buffer:
        buffer.write(str(comparison_version + 1))

    return comparison_version


if __name__ == "__main__":

    # Get comparison version from file
    savefig: bool = True
    filename: str = "earth_spherical"
    cvrsn = get_comparison_version(savefig)
    skip = ["MALARGUE", "CEBREROS", "DSS65"]

    # Data files
    data_files = [
        item for item in source_dir.glob("*.npy") if item.stem not in skip
    ]

    # Reference epoch
    j2000_tdb = time.Time("2000-01-01T12:00:00", scale="tt").tdb
    ref_epoch = 441526849.68383
    et_0 = time.TimeDelta(ref_epoch, scale="tdb", format="sec")
    initial_epoch_utc = (j2000_tdb + et_0).utc.isot  # type: ignore

    # Figure setup
    figure_setup = ng.PlotSetup(
        canvas_title="Pre-fit residuals: Closed-loop Doppler",
        xlabel=f"Hours past {initial_epoch_utc} [UTC]",
        ylabel="Residual [Hz]",
        ylim=(-0.06, 0.06),
        canvas_size=(12, 4),
        legend_location="upper right",
    )

    with ng.SingleAxis(figure_setup) as fig:

        for file in data_files:

            # Get station name
            station: str = file.name.split(".")[0]

            # Get contents
            cepochs, cobs, robs, repochs, rres, cres = np.load(file)
            cdt = (cepochs - ref_epoch) / 3600

            # Plot residuals
            fig.line(cdt, cres, fmt=".", label=station)
