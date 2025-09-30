# from nastro import graphics as ng, types as nt
from proptools import io as pio
from prefit import paths as ppaths
from pathlib import Path
from tudatpy.math import interpolators as tint
from tudatpy.astro import time_representation as ttime
import numpy as np
from matplotlib import pyplot as plt

base_dir: Path = (
    ppaths.outdir / "propagation-cl/integrator/benchmark_closest_fixed"
)
choice_dir = base_dir / "rkf78/4"

# Load results
choice = pio.PropagationOutput.from_config_file(choice_dir / "prop-config.yaml")
choice_epochs = np.array([ttime.Time(ti) for ti in choice.epochs])
choice_position_history = {
    ti: ci for ti, ci in zip(choice_epochs[::2], choice.cstate_j2000[::2, :3])
}
state_interpolator = (
    tint.create_one_dimensional_vector_interpolator_time_object(
        data_to_interpolate=choice_position_history,
        interpolator_settings=tint.lagrange_interpolation(8),
        data_first_derivatives=choice.cstate_j2000[::2, 3:],
    )
)
interpolated_position = np.array(
    [state_interpolator.interpolate(epoch) for epoch in choice_epochs]
)
error_interp = np.linalg.norm(
    choice.cstate_j2000[:, :3] - interpolated_position, axis=-1
)

fig, ax = plt.subplots()
ax.plot(choice.epochs[10:-10], error_interp[10:-10], ".")
plt.show()


# Interpolate with chosen step-size
# position_interpolator = interpolate.BarycentricInterpolator(
#     choice.epochs[::2], cstate.r_vec[:, ::2], axis=1
# )
# velocity_interpolator = interpolate.Akima1DInterpolator(
#     choice.epochs[::2], cstate.v_vec[:, ::2], axis=1
# )
# xinterp = position_interpolator(choice.epochs)
# vinterp = velocity_interpolator(choice.epochs)
# cstate_interp = nt.CartesianState[nt.Vector](*xinterp, *vinterp)

# # Difference between original and interpolated states
# cstate_err = cstate - cstate_interp

# # Interpolate choice
# fig_setup = ng.PlotSetup(yscale="log")
# with ng.SingleAxis(fig_setup) as fig:

#     fig.line(choice.epochs, cstate_err.r_mag, fmt=".")
