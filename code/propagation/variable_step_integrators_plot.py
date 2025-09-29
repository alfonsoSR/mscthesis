from proptools import config as pcon, io as pio
from prefit import paths as ppaths
from nastro import graphics as ng, types as nt
from pathlib import Path
import numpy as np
from scipy import interpolate as scint
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("integrator", help="Integrator to analyze")
parser.add_argument("-s", dest="save", action="store_true", help="Save figure")


if __name__ == "__main__":

    # Get integrator from command-line input
    args = parser.parse_args()
    integrator: str = str(args.integrator)
    save_fig: bool = bool(args.save)

    # Select benchmark with fixed step size
    benchmark_dir: Path = (
        ppaths.outdir / "propagation-cl/integrator/benchmark_fixed/rkf89/0_25"
    )
    benchmark = pio.PropagationOutput.from_config_file(
        benchmark_dir / "prop-config.yaml"
    )
    benchmark_interp = scint.Akima1DInterpolator(
        x=benchmark.epochs,
        y=benchmark.cstate_j2000,
        axis=0,
        extrapolate=False,
    )

    # Define list of configuration files from which to extract data
    base_dir: Path = (
        ppaths.outdir
        / "propagation-cl/integrator/benchmark_variable"
        / integrator
    )
    config_files: dict[float, Path] = {
        int(file.parent.name): file
        for file in base_dir.rglob(f"prop-config.yaml")
    }

    # Sort configuration files by key
    config_files = dict(sorted(config_files.items()))

    # Figure settings
    canvas_setup = ng.PlotSetup(
        canvas_title=f"Estimated integration error as function of rtol: {integrator.upper()}",
        canvas_size=(12, 7),
        show=True,
        save=save_fig,
        dir=base_dir,
        name="error-vs-tol.png",
    )

    # Subplot settings: Error as function of time
    err_setup = ng.PlotSetup(
        yscale="log",
        xlabel="Days past initial epoch",
    )
    xerr_setup = err_setup.version(
        ylabel=r"$|| \mathbf{r}(t; tol) - \mathbf{r}_{bench}(t) ||$ [m]",
    )
    verr_setup = err_setup.version(
        ylabel=(r"$||\mathbf{v}(t; tol) - \mathbf{v}(t)_{bench}||$[m/s]"),
    )

    # Subplot settings: Max error as function of step size
    max_err_setup = ng.PlotSetup(
        yscale="log",
        xscale="log",
        xlabel="rtol",
        legend_title="rtol [-log(10)]",
    )
    max_xerr_setup = max_err_setup.version(
        ylabel=r"max$(|| \mathbf{r}(t; \Delta t) - \mathbf{r}(t; \Delta t/2) ||)$ [m]",
    )
    max_verr_setup = max_err_setup.version(
        ylabel=r"max$(|| \mathbf{v}(t; \Delta t) - \mathbf{v}(t; \Delta t/2) ||)$ [m/s]",
    )

    # Initialize figure
    with ng.Mosaic("ab;cd", canvas_setup) as canvas:

        # Initialize subplots
        xerr = canvas.subplot(xerr_setup)
        xerr_max = canvas.subplot(max_xerr_setup)
        verr = canvas.subplot(verr_setup)
        verr_max = canvas.subplot(max_verr_setup)

        # Initialize lists for max_xerr, and max_verr
        max_xerr_list: list[float] = []
        max_verr_list: list[float] = []
        tol_list: list[int] = []
        feval_list: list[float] = []

        # Loop over configuration files and add data to figure
        for rtol, config_file in config_files.items():

            # Load configuration and results for current tolerance
            cconfig = pcon.PropSettings(config_file)
            cresults = pio.PropagationOutput.from_config_file(config_file)

            benchmark_filtered = benchmark_interp(cresults.epochs)

            # Calculate error between results
            cstate_step = nt.CartesianState[nt.Vector](*cresults.cstate_j2000.T)
            cstate_pstep = nt.CartesianState[nt.Vector](*benchmark_filtered.T)
            cerror = cstate_step - cstate_pstep

            # Update subplots: Error as function of time
            clabel = f"1e-{config_file.parent.name}"
            dt = (cresults.epochs - cresults.epochs[0]) / 86400.0
            xerr.line(dt, cerror.r_mag)
            verr.line(dt, cerror.v_mag)

            # Update lists with max_err data
            feval_list.append(cresults.number_of_function_evaluations)
            tol_list.append(int(config_file.parent.name))
            max_xerr_list.append(max(cerror.r_mag))
            max_verr_list.append(max(cerror.v_mag))

        # Initialize max_err subplots with grey line connecting dots
        xerr_max.line(feval_list, max_xerr_list, color="grey", alpha=0.4)
        verr_max.line(feval_list, max_verr_list, color="grey", alpha=0.4)

        # Add data to max_err subplots
        for toli, stepi, xmaxi, vmaxi in zip(
            tol_list, feval_list, max_xerr_list, max_verr_list
        ):
            xerr_max.line(stepi, xmaxi, fmt="o", label=f"{toli}")
            verr_max.line(stepi, vmaxi, fmt="o", label=f"{toli}")

        # Post-process subplots
        xerr.__exit__(0, 0, 0)
        xerr_max.__exit__(0, 0, 0)
        verr.__exit__(0, 0, 0)
        verr_max.__exit__(0, 0, 0)
