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
parser.add_argument(
    "-x", dest="show", action="store_false", help="Do now show the figure"
)


if __name__ == "__main__":

    # Get integrator from command-line input
    args = parser.parse_args()
    integrator: str = str(args.integrator)
    save_fig: bool = bool(args.save)
    show_fig: bool = bool(args.show)

    # Define list of configuration files from which to extract data
    base_dir: Path = (
        ppaths.outdir / "propagation-cl/integrator/benchmark_fixed" / integrator
    )
    config_files: dict[float, Path] = {
        float(file.parent.name.replace("_", ".")): file
        for file in base_dir.rglob(f"prop-config.yaml")
    }

    # Sort configuration files by key
    config_files = dict(sorted(config_files.items()))

    # Get smallest step size in list
    step_list: np.ndarray = np.array(list(config_files.keys()))
    min_step: float = min(step_list)

    # Figure settings
    canvas_setup = ng.PlotSetup(
        canvas_title=f"Estimated integration error as function of step size: {integrator.upper()}",
        canvas_size=(12, 7),
        show=show_fig,
        save=save_fig,
        dir=base_dir,
        name="error-vs-step-size.png",
    )

    # Subplot settings: Error as function of time
    err_setup = ng.PlotSetup(
        yscale="log",
        xlabel="Days past initial epoch",
    )
    xerr_setup = err_setup.version(
        ylabel=r"$|| \mathbf{r}(t; \Delta t) - \mathbf{r}(t; \Delta t/2) ||$ [m]",
        ylim=(1e-10, 1e4),
    )
    verr_setup = err_setup.version(
        ylabel=(
            r"$||\mathbf{v}(t; \Delta t) - \mathbf{v}(t; \Delta t/2)||$[m/s]"
        ),
        ylim=(1e-11, 1e1),
    )

    # Subplot settings: Max error as function of step size
    max_err_setup = ng.PlotSetup(
        yscale="log",
        xscale="log",
        xlabel="Step size",
        legend_title=r"$\Delta t$ [s]",
    )
    max_xerr_setup = max_err_setup.version(
        ylabel=r"max$(|| \mathbf{r}(t; \Delta t) - \mathbf{r}(t; \Delta t/2) ||)$ [m]",
        ylim=(1e-4, 1e4),
    )
    max_verr_setup = max_err_setup.version(
        ylabel=r"max$(|| \mathbf{v}(t; \Delta t) - \mathbf{v}(t; \Delta t/2) ||)$ [m/s]",
        ylim=(1e-7, 1e1),
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

        # Loop over configuration files and add data to figure
        for step, config_file in config_files.items():

            # Skip the smallest step
            if np.isclose(step, min_step, rtol=0, atol=0.01):
                continue

            # Load configuration and results for current step
            cconfig = pcon.PropSettings(config_file)
            cresults = pio.PropagationOutput.from_config_file(config_file)

            # Load configuration and results for step / 2
            pconfig_file = config_files[
                step_list[np.argmin(np.abs(step_list - (step / 2)))]
            ]
            presults = pio.PropagationOutput.from_config_file(pconfig_file)

            # Interpolate output with smallest step size
            pinterp = scint.Akima1DInterpolator(
                x=presults.epochs,
                y=presults.cstate_j2000,
                axis=0,
                extrapolate=False,
            )
            presults_filtered = pinterp(cresults.epochs)

            # Calculate error between results
            cstate_step = nt.CartesianState[nt.Vector](*cresults.cstate_j2000.T)
            cstate_pstep = nt.CartesianState[nt.Vector](*presults_filtered.T)
            cerror = cstate_step - cstate_pstep

            # Update subplots: Error as function of time
            clabel = config_file.parent.name
            dt = (cresults.epochs - cresults.epochs[0]) / 86400.0
            xerr.line(dt, cerror.r_mag)
            verr.line(dt, cerror.v_mag)

            # Update lists with max_err data
            max_xerr_list.append(max(cerror.r_mag))
            max_verr_list.append(max(cerror.v_mag))

        # Initialize max_err subplots with grey line connecting dots
        xerr_max.line(step_list[1:], max_xerr_list, color="grey", alpha=0.4)
        verr_max.line(step_list[1:], max_verr_list, color="grey", alpha=0.4)

        # Add data to max_err subplots
        for stepi, xmaxi, vmaxi in zip(
            step_list[1:], max_xerr_list, max_verr_list
        ):
            xerr_max.line(stepi, xmaxi, fmt="o", label=f"{stepi:.2f}")
            verr_max.line(stepi, vmaxi, fmt="o", label=f"{stepi:.2f}")

        # Post-process subplots
        xerr.__exit__(0, 0, 0)
        xerr_max.__exit__(0, 0, 0)
        verr.__exit__(0, 0, 0)
        verr_max.__exit__(0, 0, 0)
