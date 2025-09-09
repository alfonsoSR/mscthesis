from nastro import graphics as ng
import numpy as np
from pathlib import Path
from prefit import paths as ppaths


source_file = ppaths.outdir / "prefit-closed-loop/NWNORCIA.npy"
cepochs, cobs, robs, repochs, rres, cres = np.load(source_file)
cdt = (cepochs - cepochs[0]) / 3600
rdt = (repochs - repochs[0]) / 3600

outdir = ppaths.figdir / "individual-stations-prefit-cl"
outdir.mkdir(exist_ok=True)
outpath = outdir / "NWNORCIA.png"

setup = ng.PlotSetup(
    save=True,
    dir=outpath.parent,
    name=outpath.name,
    canvas_title=f"Analysis: {outpath.stem}",
)

with ng.Mosaic("a;b", setup) as fig:

    with fig.subplot() as both:

        both.line(cdt, cres, fmt=".", label="Current")
        both.line(rdt, rres, fmt=".", label="Reference")

    with fig.subplot() as comp:

        comp.line(rdt, rres - cres, fmt=".")
