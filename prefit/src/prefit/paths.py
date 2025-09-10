from pathlib import Path

basedir: Path = Path("/Users/alfonso/home/thesis/phase-a")
codedir = basedir / "code"
datadir = basedir / "data/mex_phobos_flyby"
kerneldir = Path().home() / ".spice_kernels/mex/spk"
psadir = Path().home() / ".psa/mex/data"
outdir = basedir / "output"
figdir = basedir / "figures"
