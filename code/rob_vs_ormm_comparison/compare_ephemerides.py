from nastro import graphics as ng, types as nt, catalog as nc
import numpy as np
from pathlib import Path
from prefit import paths as ppaths

sourcedir = ppaths.outdir / "rob_ormm_comparison"
outdir = ppaths.figdir / "rob_ormm_comparison"
outdir.mkdir(exist_ok=True)

xpos_mex_mars_rsw_ormm = np.load(sourcedir / "ormm_rsw.npy").T
xpos_mex_mars_rsw_rob = np.load(sourcedir / "rob_rsw.npy").T

cstate_mex_mars_j2000_ormm = nt.CartesianState(
    *np.load(sourcedir / "ormm_icrf.npy").T
)
cstate_mex_mars_j2000_rob = nt.CartesianState(
    *np.load(sourcedir / "rob_icrf.npy").T
)

epochs = np.load(sourcedir / "epochs.npy")
dt = (epochs - epochs[0]) / 86400.0


# Get error between mars-ICRF states
error_icrf = cstate_mex_mars_j2000_ormm - cstate_mex_mars_j2000_rob
icrf_setup = ng.PlotSetup(
    canvas_title="Difference between MEX coordinates in ICRF: ORMM vs ROB",
    xlabel="Days past 2013-12-28T00:00:00",
    ylabel="Magnitude of position error [m]",
    save=True,
    show=False,
    dir=outdir,
    name="icrf-error-rob-ormm.pdf",
)
with ng.SingleAxis(icrf_setup) as fig:
    fig.line(dt, error_icrf.r_mag)


# Error in RSW frame
rsw_labels = ("R", "S", "W")
rsw_setup = ng.PlotSetup(
    canvas_size=(6, 7),
    canvas_title="Difference between MEX coordinates in RSW: ORMM vs ROB",
    save=True,
    show=False,
    dir=outdir,
    name="rsw-error-rob-ormm.pdf",
)
with ng.Mosaic("a;b;c", setup=rsw_setup) as fig:

    rsw_labels = ("R", "S", "w")

    for idx, label in enumerate(rsw_labels):

        setup = ng.PlotSetup(
            xlabel="Days past 2013-12-28T00:00:00",
            ylabel=r"$\epsilon" + f"_{label}$",
        )

        with fig.subplot(setup=setup) as subfig:

            subfig.line(
                dt, xpos_mex_mars_rsw_ormm[idx] - xpos_mex_mars_rsw_rob[idx]
            )
