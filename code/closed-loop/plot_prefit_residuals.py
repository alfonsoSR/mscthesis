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
    savefig: bool = False
    filename: str = "email_martin"
    cvrsn = get_comparison_version(savefig)
    skip = ["MALARGUE", "CEBREROS", "DSS65"]

    # Data files
    data_files = [
        item for item in source_dir.glob("*.npy") if item.stem not in skip
    ]

    # Generate mosaic based on number of stations
    # if len(data_files) == 0:
    #     raise ValueError("Nothing to plot")
    # mosaic = ""
    # for item in data_files:
    #     mosaic += f"{next(LETTERS)}{next(LETTERS)};"
    # mosaic = mosaic[:-1]

    mosaic = ";".join([f"{next(LETTERS)}{next(LETTERS)}" for _ in data_files])

    # Canvas setup
    figdir.mkdir(exist_ok=True)
    canvas_setup = ng.PlotSetup(
        canvas_title=f"Variation with respect to IFMS residuals: v{cvrsn}",
        canvas_size=(12, 8),
        save=savefig,
        show=True,
        dir=figdir,
        name=f"{filename}_v{cvrsn}.png",
    )

    with ng.Mosaic(mosaic, setup=canvas_setup) as fig:

        for file in data_files:

            # Get station name
            station: str = file.name.split(".")[0]

            # Get Luigi's residuals
            luigi: np.ndarray | None = None
            try:
                luigi = np.load(source_dir / f"luigi/{station}.npy")
            except FileNotFoundError:
                print(f"Luigi's pre-fit not available for {station}")

            # Get contents
            cepochs, cobs, robs, repochs, rres, cres = np.load(file)
            cdt = (cepochs - cepochs[0]) / 3600
            rdt = (repochs - repochs[0]) / 3600

            # Check that epochs are coincident
            print(f"Station {station}: Check on epoch coincidence")
            print(f"Number of mismatches: {np.sum(np.abs(cdt - rdt) > 2e-15)}")

            # Get iso string representation of initial epoch
            j2000_tdb = time.Time("2000-01-01T12:00:00", scale="tt").tdb
            et_0 = time.TimeDelta(cepochs[0], scale="tdb", format="sec")
            initial_epoch_utc = (j2000_tdb + et_0).utc.isot  # type: ignore

            # Use specific ylimits for NWNORCIA
            comp_ylim: tuple[float, float] | None = (-0.15, 0.15)
            if station == "NWNORCIA":
                comp_ylim = (-0.15, 0.15)

            # Set up subplot
            error_setup = ng.PlotSetup(
                axtitle=f"Station: {station}",
                xlabel=f"Hours past {initial_epoch_utc} [TDB]",
                ylabel="IFMS res - computed res [Hz]",
                ylim=comp_ylim,
            )

            comp_setup = ng.PlotSetup(
                axtitle=f"Station: {station}",
                xlabel=f"Hours past {initial_epoch_utc} [TDB]",
                ylabel="Residual [Hz]",
                ylim=comp_ylim,
            )

            with fig.subplot(setup=comp_setup) as subfig:
                subfig.line(cdt, cres, fmt=".", label="Computed")
                subfig.line(rdt, rres, fmt=".", label="From IFMS")
                # if luigi is not None:
                #     subfig.line(rdt, luigi, fmt=".", label="Luigi")
                # subfig.line(rdt, cres, fmt=".", label="Predicted")
                # subfig.line(cdt, luigi, fmt=".", label="Luigi")
                # subfig.line(cdt, cobs, fmt=".", label="Predicted")

            with fig.subplot(setup=error_setup) as subfig:
                subfig.line(cdt, rres - cres, fmt=".")
                # subfig.line(cdt, luigi - rres, fmt=".", label="Luigi")
