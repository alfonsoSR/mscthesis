from tudatpy.interface import spice
from tudatpy.dynamics.environment_setup import ephemeris
from tudatpy.astro import time_representation as ttime
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

metakernel = str(
    Path(
        "/Users/alfonso/home/thesis/phase-a/output/propagation-cl/integrator/benchmark_fixed/rkf78/8/metak.tm"
    )
)

try:
    spice.load_kernel(metakernel)

    tref = ttime.DateTime.from_iso_string(
        "2013-12-29 07:10:09.000000000000000"
    ).to_epoch_time_object()
    buffer = ttime.Time(60.0)
    t0 = tref - buffer
    tend = tref + buffer

    # t0 = ttime.DateTime.from_iso_string(
    #     "2013-12-27T00:00:00"
    # ).to_epoch_time_object()
    # tend = ttime.DateTime.from_iso_string(
    #     "2014-01-02T00:00:00"
    # ).to_epoch_time_object()
    step = ttime.Time(0.0001)

    epochs = [t0]
    while epochs[-1] <= tend:
        epochs.append(epochs[-1] + step)
    float_epochs = [ti.to_float() for ti in epochs]

    # eph_set = ephemeris.direct_spice(
    #     frame_orientation="J2000",
    #     frame_origin="Phobos",
    # )
    # eph = ephemeris.create_ephemeris(eph_set, "MEX")

    distances = np.linalg.norm(
        np.array(
            [
                # eph.cartesian_position(ti)
                spice.get_body_cartesian_position_at_epoch(
                    "MEX", "Phobos", "J2000", "NONE", ti
                )
                for ti in epochs
            ]
        ),
        axis=1,
    )

    min_distance_idx = np.argmin(distances)

    print(
        ttime.DateTime.from_epoch_time_object(
            epochs[min_distance_idx]
        ).to_iso_string()
    )

    print(distances[min_distance_idx])

    fig, ax = plt.subplots()
    ax.plot(float_epochs, distances)
    ax.plot(float_epochs[min_distance_idx], distances[min_distance_idx], "o")
    plt.show()

finally:
    spice.clear_kernels()
