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
        initial_utc_epoch=ttime.DateTime.from_iso_string("2013-12-28T18:19:43"),
        step=ttime.Time(10.0),
        duration_in_steps=60 * 24 * 1,
    )

    # Path to metakernel
    metak = ppaths.datadir / "metak_mex.tm"

    # NAIF codes for relevant components
    hga_id: str = "-41020"  # High-gain antenna
    marsis_id: str = "-41300"  # MARSIS
    launch_vehicle_interface_id: str = "-41"  # This is the same as MEX

    # NAIF codes for relevant frames
    mex_sc_frameid: str = ""  # Spacecraft reference frame
    mex_spacecraft_frameid: str = "MEX_SPACECRAFT"  # Mechanical/structure frame
    mex_hga_frameid: str = "MEX_HGA"  # HGA frame

    marsis_frame_id: str = "MEX_MARSIS"  # Marsis frame wr

    try:
        spice.load_kernel(str(metak))

        kloaded = putils.get_loaded_kernels("all")
        for ki in kloaded:
            print(ki.name)

        # Get position of LVI with respect to Mars barycenter
        lvi_mars_j2000 = np.array(
            [
                spice.get_body_cartesian_position_at_epoch(
                    target_body_name=launch_vehicle_interface_id,
                    observer_body_name="MARS BARYCENTER",
                    reference_frame_name="J2000",
                    aberration_corrections="NONE",
                    ephemeris_time=eti,
                )
                for eti in epochs_et
            ]
        )

        # rotation = spice.compute_rotation_matrix_between_frames(
        #     original_frame=mex_spacecraft_frameid,
        #     new_frame=mex_hga_frameid,
        #     ephemeris_time=epochs_et[0],
        # )
        # print(rotation)

        # Get position of HGA wrt to LVI
        hga_lvi_mexsc = np.array(
            [
                spice.get_body_cartesian_position_at_epoch(
                    target_body_name=hga_id,
                    observer_body_name=launch_vehicle_interface_id,
                    reference_frame_name=mex_spacecraft_frameid,
                    aberration_corrections="NONE",
                    ephemeris_time=eti,
                )
                for eti in epochs_et
            ]
        )

        # Get position of MARSIS origin wrt to LVI
        marsis_lvi_mexsc = np.array(
            [
                spice.get_body_cartesian_position_at_epoch(
                    target_body_name=launch_vehicle_interface_id,
                    observer_body_name=marsis_id,
                    reference_frame_name=marsis_frame_id,
                    aberration_corrections="NONE",
                    ephemeris_time=eti,
                )
                for eti in epochs_et
            ]
        )

    finally:
        spice.clear_kernels()

    outdir = ppaths.outdir / "structure_analysis"
    outdir.mkdir(exist_ok=True)
    np.save(outdir / "lvi_mars_j2000.npy", lvi_mars_j2000)
    np.save(outdir / "hga_lvi_mexframe.npy", hga_lvi_mexsc)
    np.save(outdir / "marsis_lvi_mexframe.npy", marsis_lvi_mexsc)

    # diff = xpos_mex_mars_j2000_rob - xpos_mex_mars_j2000_ormm

    # fig, ax = plt.subplots(layout="tight")
    # ax.plot(epochs_et_float, np.linalg.norm(diff[:, :3], axis=-1))
    # plt.show()

    # with spiceypy.KernelPool(metak):

    #     print(spice.get_total_count_of_kernels_loaded())

    #     dsk_kernels = get_loaded_kernels("all")
    #     print([dski.name for dski in dsk_kernels])
