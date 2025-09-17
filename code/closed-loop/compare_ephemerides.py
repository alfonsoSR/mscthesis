from pathlib import Path
from prefit import paths as ppaths, io as pio
from nastro import graphics as ng
import spiceypy as spice
import numpy as np

basedir = ppaths.outdir / "prefit-cl/relativistic-correction"
enddir = "single-bodies/sun"

refdir = basedir / "ormm" / enddir
comparedir = basedir / "rob" / enddir


results_per_station: dict[str, tuple[pio.PrefitResults, pio.PrefitResults]] = {}
for rfile in refdir.glob("*.npy"):

    cfile = comparedir / rfile.name

    results_per_station[rfile.stem] = (
        pio.PrefitResults(cfile),
        pio.PrefitResults(rfile),
    )


poserr_per_station = {}
for station, (cdata, rdata) in results_per_station.items():

    # ROB state
    with spice.KernelPool(
        str(ppaths.kerneldir / "MEX_ROB_130101_131231_001.BSP")
    ):

        robpos = np.array(
            spice.spkezr("MEX", cdata.cepochs, "J2000", "NONE", "Mars")[0]
        ).T

    # ORMM state
    with spice.KernelPool(
        str(ppaths.kerneldir / "ORMM_T19_131201000000_01033.BSP")
    ):

        ormmpos = np.array(
            spice.spkezr("MEX", cdata.cepochs, "J2000", "NONE", "Mars")[0]
        ).T

    poserr_per_station[station] = np.linalg.norm((ormmpos - robpos)[:3], axis=0)


lim = 0.08
figure_setup = ng.PlotSetup(
    canvas_title=f"Relativistic effects: {comparedir.name}",
    canvas_size=(12, 6),
)

residual_setup = ng.PlotSetup(
    ylabel="delta residual [Hz]",
    xlabel="epoch",
    ylim=(-lim, lim),
)
error_setup = ng.PlotSetup(ylabel="delta eph [km]", xlabel="epoch")

with ng.Mosaic("a;b", figure_setup) as fig:

    with fig.subplot(residual_setup) as rfig:

        for station, (cdata, rdata) in results_per_station.items():

            rfig.line(
                cdata.cepochs, rdata.cres - cdata.cres, fmt=".", label=station
            )

    with fig.subplot(error_setup) as efig:

        for station, error in poserr_per_station.items():

            efig.line(
                results_per_station[station][0].cepochs,
                error,
                fmt=".",
                label=station,
            )

# with ng.SingleAxis(figure_setup) as fig:

#     for rfile in refdir.glob("*.npy"):

#         cfile = comparedir / rfile.name

#         cdata = pio.PrefitResults(cfile)
#         rdata = pio.PrefitResults(rfile)

#         fig.line(
#             cdata.cepochs, rdata.cres - cdata.cres, fmt=".", label=rfile.stem
#         )
