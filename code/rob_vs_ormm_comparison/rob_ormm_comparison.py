from pathlib import Path
from tudatpy.interface import spice
from tudatpy.astro import (
    time_representation as ttime,
    frame_conversion as tframe,
)
from prefit import utils as putils, paths as ppaths
import numpy as np
from matplotlib import pyplot as plt


def get_epochs_in_et(
    initial_utc_epoch: ttime.DateTime, step: ttime.Time, duration_in_steps: int
) -> list[ttime.Time]:

    # Time conversion settings
    time_converter = ttime.default_time_scale_converter()
    j2000_tt = ttime.DateTime.from_iso_string(
        "2000-01-01T12:00:00"
    ).to_epoch_time_object()
    j2000_tdb = time_converter.convert_time_object(
        ttime.TimeScales.tt_scale, ttime.TimeScales.tdb_scale, j2000_tt
    )

    # Time interval (UTC)
    initial_epoch = initial_utc_epoch.to_epoch_time_object()
    epochs_utc = [
        initial_epoch + step * i for i in range(duration_in_steps + 1)
    ]

    # Time interval ET (TDB)
    epochs_tdb = [
        time_converter.convert_time_object(
            ttime.TimeScales.utc_scale, ttime.TimeScales.tdb_scale, epoch
        )
        for epoch in epochs_utc
    ]
    epochs_et = [epoch - j2000_tdb for epoch in epochs_tdb]

    return epochs_et


if __name__ == "__main__":

    # Define time interval
    epochs_et = get_epochs_in_et(
        initial_utc_epoch=ttime.DateTime.from_iso_string("2013-12-28T00:00:00"),
        step=ttime.Time(60.0),
        duration_in_steps=60 * 24 * 3,
    )

    # Get position from ROB ephemerides
    try:
        spice.load_kernel(
            str(ppaths.kerneldir / "MEX_ROB_130101_131231_001.BSP")
        )

        # Get position of MEX
        cstate_mex_mars_j2000_rob = np.array(
            [
                spice.get_body_cartesian_state_at_epoch(
                    target_body_name="MARS EXPRESS",
                    observer_body_name="MARS",
                    reference_frame_name="J2000",
                    aberration_corrections="NONE",
                    ephemeris_time=eti,
                )
                for eti in epochs_et
            ]
        )

        # Convert to RSW
        icrf2rsw_rob = np.array(
            [
                tframe.inertial_to_rsw_rotation_matrix(xi)
                for xi in cstate_mex_mars_j2000_rob
            ]
        )
        xpos_mex_mars_rsw_rob = np.array(
            [
                icrf2rsw_rob[i] @ cstatei[:3]
                for i, cstatei in enumerate(cstate_mex_mars_j2000_rob)
            ]
        )

    finally:
        spice.clear_kernels()

    # Get position from Luigi's ephemerides
    try:
        spice.load_kernel(str(ppaths.datadir / "metak_mex.tm"))
        # spice.load_kernel(str(kdir / "ORMM_T19_140101000000_01041.BSP"))

        # Get position of MEX
        cstate_mex_mars_j2000_ormm = np.array(
            [
                spice.get_body_cartesian_state_at_epoch(
                    target_body_name="MARS EXPRESS",
                    observer_body_name="MARS",
                    reference_frame_name="J2000",
                    aberration_corrections="NONE",
                    ephemeris_time=eti,
                )
                for eti in epochs_et
            ]
        )

        # Convert to RSW
        icrf2rsw_ormm = np.array(
            [
                tframe.inertial_to_rsw_rotation_matrix(xi)
                for xi in cstate_mex_mars_j2000_ormm
            ]
        )
        xpos_mex_mars_rsw_ormm = np.array(
            [
                icrf2rsw_ormm[i] @ cstatei[:3]
                for i, cstatei in enumerate(cstate_mex_mars_j2000_ormm)
            ]
        )

    finally:
        spice.clear_kernels()

    epochs_et_float = np.array([eti.to_float() for eti in epochs_et])

    outdir = ppaths.outdir / "rob_ormm_comparison"
    outdir.mkdir(exist_ok=True)
    np.save(outdir / "ormm_rsw", xpos_mex_mars_rsw_ormm)
    np.save(outdir / "ormm_icrf", cstate_mex_mars_j2000_ormm)
    np.save(outdir / "rob_rsw", xpos_mex_mars_rsw_rob)
    np.save(outdir / "rob_icrf", cstate_mex_mars_j2000_rob)
    np.save(outdir / "epochs", epochs_et_float)

    # diff = xpos_mex_mars_j2000_rob - xpos_mex_mars_j2000_ormm

    # fig, ax = plt.subplots(layout="tight")
    # ax.plot(epochs_et_float, np.linalg.norm(diff[:, :3], axis=-1))
    # plt.show()

    # with spiceypy.KernelPool(metak):

    #     print(spice.get_total_count_of_kernels_loaded())

    #     dsk_kernels = get_loaded_kernels("all")
    #     print([dski.name for dski in dsk_kernels])
