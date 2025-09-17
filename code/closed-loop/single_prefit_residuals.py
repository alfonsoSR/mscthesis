from nastro import graphics as ng
from pathlib import Path
import csv
import numpy as np
from astropy import time
from prefit import paths as ppaths, io as pio
import argparse
import datetime

current_epoch = datetime.datetime.now().isoformat(sep="T")


parser = argparse.ArgumentParser()
parser.add_argument("config_file", help="Configuration file")
LETTERS = iter([x for x in "abcdefghijklmnopqrstuvwxyz"])


def get_comparison_version(savefig_flag: bool) -> str:

    cache_file = Path(__file__).parent / "comparison_version.txt"
    with cache_file.open() as buffer:
        comparison_version = int(buffer.readline().strip())

    if not savefig_flag:
        return f"v{comparison_version}"

    # Update version
    with cache_file.open("w") as buffer:
        buffer.write(str(comparison_version + 1))

    return f"v{comparison_version}"


if __name__ == "__main__":

    # Get configuration
    config_path = Path(parser.parse_args().config_file).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Not found: {config_path}")

    config = pio.load_configuration(config_path)
    plot_config = config.plotting
    config_relpath = config_path.relative_to(ppaths.outdir)

    # Define paths
    source_dir = config_path.parent
    figdir = source_dir
    figdir.mkdir(exist_ok=True)

    # Get version of plot
    filename_base: str = "residuals-single"
    previous_versions = [xi for xi in figdir.glob("*residuals-single*.png")]
    if len(previous_versions) != 0:
        cvrsn = f"v{len(previous_versions)}"
    else:
        cvrsn = "v0"

    # savefig = bool(args.save_fig)
    # filename = str(args.label)
    # cvrsn = get_comparison_version(savefig)

    # source_dir = ppaths.outdir / "prefit-closed-loop"
    # figdir = ppaths.figdir / "prefit-closed-loop"
    # figdir.mkdir(exist_ok=True)

    # Data files
    skip = plot_config["ignore_stations"]
    data_files = [
        item for item in source_dir.glob("*.npy") if item.stem not in skip
    ]

    # Load data per station
    data_per_station = {
        file.name.split(".")[0]: np.load(file) for file in data_files
    }

    # Reference epoch
    all_epochs = np.concatenate([item[0] for item in data_per_station.values()])
    all_epochs.sort()
    ref_epoch = all_epochs[0]
    j2000_tdb = time.Time("2000-01-01T12:00:00", scale="tt").tdb
    et_0 = time.TimeDelta(ref_epoch, scale="tdb", format="sec")
    initial_epoch_utc = (j2000_tdb + et_0).utc.isot  # type: ignore

    # Figure setup
    figure_setup = ng.PlotSetup(
        canvas_title=f"Closed-loop pre-fit residuals :: {config_relpath}",
        xlabel=f"Hours past {initial_epoch_utc} [UTC]",
        ylabel="Residual [Hz]",
        ylim=(-0.025, 0.025),
        canvas_size=(12, 4),
        legend_location="upper right",
        save=config.plotting["save"],
        show=config.plotting["show"],
        dir=figdir,
        name=f"{filename_base}-{cvrsn}.png",
    )

    with ng.SingleAxis(figure_setup) as fig:

        for station, station_data in data_per_station.items():

            # Get contents
            cepochs, cobs, robs, repochs, rres, cres = station_data
            cdt = (cepochs - ref_epoch) / 3600

            # Plot residuals
            fig.line(cdt, cres, fmt=".", label=station)
