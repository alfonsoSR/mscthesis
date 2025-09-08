from nastro import graphics as ng
from pathlib import Path
import csv
import numpy as np
from astropy import time
from prefit import paths as ppaths

source_dir = ppaths.outdir / "prefit-closed-loop"
figdir = ppaths.figdir / "prefit-closed-loop"
figdir.mkdir(exist_ok=True)


def get_content_of_residuals_file(source_file: Path) -> np.ndarray:

    # Read epochs and residuals into array
    with source_file.open() as buffer:

        # Skip header (First two lines)
        for _ in range(2):
            buffer.readline()

        epochs: list[float] = []
        residuals: list[float] = []
        for line in buffer:
            epoch, _, residual = line.strip().split(",")
            epochs.append(float(epoch))
            residuals.append(float(residual))

    out = np.array([epochs, residuals])

    # Return sorted by epoch
    return out[:, np.argsort(out[0])]


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
    filename: str = "no_com_correction"
    cvrsn = get_comparison_version(savefig)

    # Data files
    data_files = [item for item in source_dir.glob("*.npy")]

    # Generate mosaic based on number of stations
    mosaic = ";".join([f"{i}{i+1}" for i in range(1, len(data_files) + 3, 2)])

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
            luigi = np.load(source_dir / f"luigi/{station}.npy")

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

            # Set up subplot
            error_setup = ng.PlotSetup(
                axtitle=f"Station: {station}",
                xlabel=f"Hours past {initial_epoch_utc} UTC",
                ylabel="Difference in residuals [Hz]",
            )

            # Use specific ylimits for NWNORCIA
            comp_ylim: tuple[float, float] | None = None
            if station == "NWNORCIA":
                comp_ylim = (-0.15, 0.15)

            comp_setup = ng.PlotSetup(
                axtitle=f"Station: {station}",
                xlabel=f"Hours past {initial_epoch_utc} [TDB]",
                ylabel="Residual [Hz]",
                ylim=comp_ylim,
            )

            with fig.subplot(setup=error_setup) as subfig:
                subfig.line(cdt, cres - rres, fmt=".")

            with fig.subplot(setup=comp_setup) as subfig:
                subfig.line(rdt, rres, fmt=".", label="Reference")
                subfig.line(rdt, luigi, fmt=".", label="Luigi")
                subfig.line(cdt, cres, fmt=".", label="Current")
