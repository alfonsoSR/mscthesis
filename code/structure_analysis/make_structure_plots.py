from nastro import graphics as ng, types as nt
import numpy as np
from pathlib import Path
import spiceypy as spice
from astropy import time
from prefit import paths as ppaths

source_dir = ppaths.outdir / "structure_analysis"
metak = ppaths.datadir / "metak_mex.tm"

if __name__ == "__main__":

    epoch = time.Time("2013-12-28T18:19:43", scale="tdb")
    j2000 = time.Time("2000-01-01T12:00:00", scale="tt").tdb
    et = (epoch - j2000).to("s").value  # type: ignore

    ets = np.arange(et, et + (2 * 86400) + 1, 86400)

    with spice.KernelPool(str(metak)):

        # Position of LVI wrt Mars in J2000
        lvi_mars_j2000 = (
            np.array(spice.spkezr("-41020", et, "J2000", "NONE", "499")[0])[:3]
            * 1e3
        )
        print(lvi_mars_j2000)

        # Position of the HGA wrt LVI in MEX_SPACECRAFT
        hga_lvi_mexframe = (
            np.array(spice.spkez(-41020, et, "MEX_SPACECRAFT", "NONE", -41)[0])[
                :3
            ]
            * 1e3
        )
        print(hga_lvi_mexframe)

        # Rotation from MEX_SPACECRAFT to J2000
        rotation = np.array(spice.pxform("MEX_SPACECRAFT", "J2000", et))
        hga_lvi_j2000_rot = rotation @ hga_lvi_mexframe
        print(hga_lvi_j2000_rot)

        hga_lvi_j2000 = (
            np.array(spice.spkez(-41020, et, "J2000", "NONE", -41)[0])[:3] * 1e3
        )
        print(hga_lvi_j2000 - hga_lvi_j2000_rot)

        # # Position of HGA wrt LVI in MEX_SPACECRAFT
        # mars_com_marsis = (
        #     np.array(
        #         [
        #             spice.spkez(499, eti, "J2000", "NONE", -41000)[0][:3]
        #             for eti in ets
        #         ]
        #     ).T
        #     * 1e3
        # )
        # print(mars_com_marsis.T[0])
        # exit(0)

        # # Position of typical center wrt LVI in MEX_SPACECRAFT
        # marsis_mars_marsis = (
        #     np.array(
        #         [
        #             spice.spkez(4, eti, "MEX_MARSIS", "NONE", -41300)[0][:3]
        #             for eti in ets
        #         ]
        #     ).T
        #     * 1e3
        # )
        # print(marsis_mars_marsis.T)
        # exit(0)

#         # LVI wrt MARSIS in MARSIS frame
#         marsis_lvi_marsis = (
#             np.array(
#                 [
#                     spice.spkez(-41020, eti, "MEX_SPACECRAFT", "NONE", -41)[0][
#                         :3
#                     ]
#                     for eti in ets
#                 ]
#             ).T
#             * 1e3
#         )
#         print(marsis_lvi_marsis)
#         exit(0)

#         corners = {}
#         corners_marsis = {}
#         for idx, corner in enumerate(
#             (-410012, -41003, -41004, -41005, -410016, -41007, -41008, -41009)
#         ):

#             corners[f"c{idx+1}"] = (
#                 np.array(
#                     [
#                         spice.spkez(corner, eti, "MEX_SPACECRAFT", "NONE", -41)[
#                             0
#                         ][:3]
#                         for eti in ets
#                     ]
#                 ).T
#                 * 1e3
#             )
#             corners_marsis[f"c{idx+1}"] = (
#                 np.array(
#                     [
#                         spice.spkez(corner, eti, "j2000", "NONE", -41)[0][:3]
#                         for eti in ets
#                     ]
#                 ).T
#                 * 1e3
#             )

#     boxes = np.array([corner for corner in corners.values()]).swapaxes(0, -1)
#     boxesm = np.array([corner for corner in corners_marsis.values()]).swapaxes(
#         0, -1
#     )

#     # lvi_mars_j2000 = np.load(source_dir / "lvi_mars_j2000.npy").T
#     # hga_lvi_mexframe = np.load(source_dir / "hga_lvi_mexframe.npy").T
#     # marsis_lvi_mexframe = np.load(source_dir / "marsis_lvi_mexframe.npy").T

#     figsetup = ng.PlotSetup(
#         xlabel="x_sc",
#         ylabel="y_sc",
#         zlabel="z_sc",
#     )

#     with ng.Plot3D(figsetup) as fig:

#         fig.line(0, 0, 0, fmt="o", label="LVI")
#         # fig.line(*hga_lvi_mexframe, fmt="o", label="HGA")
#         # fig.line(*marsis_lvi_mexframe, fmt="o", label="marsis")
#         # fig.line(*marsis_lvi_marsis, label="lvi")

#         for box, boxm in zip(boxes, boxesm):
#             # fig.line(*box, color="black")
#             fig.line(*boxm, fmt="o")

#         # fig.line(*marsis_lvi_mexframe, fmt="o", label="MARSIS")

#     # diff = np.linalg.norm(lvi_mars_j2000 - mex_mars_j2000, axis=0)

#     # with ng.SingleAxis() as fig:
#     #     fig.line(diff)

#     # with ng.Plot3D(ng.PlotSetup()) as fig:

#     #     fig.line(*diff)

#     #     # fig.line(0, 0, 0, fmt="o", label="Mars")
#     # fig.line(*lvi_mars_j2000, label="LVI")
#     # fig.line(*mex_mars_j2000, label="MEX")


# # antenna = np.load(src / "antenna.npy").T

# # with ng.Plot3D(ng.PlotSetup()) as fig:

# #     fig.line(*antenna, fmt="o")
