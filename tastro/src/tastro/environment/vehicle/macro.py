from tudatpy.dynamics.environment_setup import (
    rotation_model as trots,
    radiation_pressure as trad,
    vehicle_systems as tvs,
)
import numpy as np
from ...logging import log
import typing
import traceback

if typing.TYPE_CHECKING:
    from ...config import CaseSetup
    from ...config.environment.vehicles import VehicleSetup
    from ...config.environment.vehicles import (
        SinglePanelShapeSettings,
        VehicleShapeSetup,
        VehicleRadiationTargetSetup,
    )


def front_back_panel_geometry(
    area: float, frame: str, normal: typing.Literal["x", "y", "z"]
) -> tuple[tvs.BodyPanelGeometrySettings, tvs.BodyPanelGeometrySettings]:

    match normal:

        case "x":
            normal_vector = np.array([1.0, 0.0, 0.0])
        case "y":
            normal_vector = np.array([0.0, 1.0, 0.0])
        case "z":
            normal_vector = np.array([0.0, 0.0, 1.0])
        case _:
            raise ValueError(f"Invalid normal: {normal}")

    positive_panel = tvs.frame_fixed_panel_geometry(
        surface_normal=normal_vector,
        area=area,
        frame_orientation=frame,
    )
    negative_panel = tvs.frame_fixed_panel_geometry(
        surface_normal=-normal_vector,
        area=area,
        frame_orientation=frame,
    )

    return positive_panel, negative_panel


def create_single_panel_geometry_from_config(
    panel_config: "SinglePanelShapeSettings",
) -> "tvs.BodyPanelGeometrySettings":

    log.debug(f"Geometry {panel_config.panel_id} :: {panel_config.frame}")

    # Define settings for panel geometry
    return tvs.frame_fixed_panel_geometry(
        surface_normal=panel_config.normal_direction,
        area=panel_config.area,
        frame_orientation=panel_config.frame,
    )


def create_paneled_geometry_from_config(shape_config: "VehicleShapeSetup"):

    # Ensure chosen shape model is paneled
    if shape_config.model != "paneled":
        log.fatal(
            f"Requested generation of paneled geometry with "
            f"{shape_config.model} shape settings"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Ensure paneled settings are defined
    if len(shape_config.paneled_settings) == 0:
        log.fatal(
            f"Requested generation of paneled geometry without "
            f"specifying any panels"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Generate settings for individual panels and return
    panel_collection: dict[str, "tvs.BodyPanelGeometrySettings"] = {
        panel_config.panel_id: create_single_panel_geometry_from_config(
            panel_config
        )
        for panel_config in shape_config.paneled_settings
    }
    return panel_collection


def create_paneled_reflection_laws_from_config(
    radiation_config: "VehicleRadiationTargetSetup",
) -> dict[str, trad.BodyPanelReflectionLawSettings]:

    # Ensure chosen radiation model is paneled
    if radiation_config.model != "paneled":
        log.fatal(
            f"Requested definition of paneled reflection laws with "
            f"{radiation_config.model} radiation settings"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Ensure paneled settings are defined
    if len(radiation_config.paneled_settings) == 0:
        log.fatal(
            f"Requested generation of paneled radiation interface without "
            f"specifying any panels"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Define reflection laws per panel
    reflection_laws: dict[str, trad.BodyPanelReflectionLawSettings] = {}
    for settings in radiation_config.paneled_settings:

        # Define reflection law according to settings
        reflection_law = trad.specular_diffuse_body_panel_reflection(
            specular_reflectivity=settings.specular_reflectivity,
            diffuse_reflectivity=settings.diffuse_reflectivity,
            with_instantaneous_reradiation=settings.instantaneous_reradiation,
        )

        # Assign to all affected panels
        for panel_id in settings.panels:

            log.debug(
                f"Reflection {panel_id} :: {settings.instantaneous_reradiation}"
            )

            reflection_laws[panel_id] = reflection_law

    return reflection_laws


def create_paneled_rotation_models_from_config(
    shape_config: "VehicleShapeSetup", base_frame: str
) -> dict[str, "trots.RotationModelSettings"]:

    # Ensure chosen shape model is paneled
    if shape_config.model != "paneled":
        log.fatal(
            f"Attempted to define panel orientations with "
            f"{shape_config.model} shape settings"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Ensure paneled settings are defined
    if len(shape_config.paneled_settings) == 0:
        log.fatal(
            f"Requested generation of paneled geometry without "
            f"specifying any panels"
        )
        log.fatal(traceback.extract_stack()[-2])
        exit(1)

    # Generate settings for individual panels and return
    rotation_laws: dict[str, "trots.RotationModelSettings"] = {}
    for panel_config in shape_config.paneled_settings:

        # Skip if geometry is defined in base frame
        if panel_config.frame == base_frame:
            continue

        # Skip if frame is already present
        if panel_config.frame in rotation_laws:
            continue

        log.debug(f"Rotation {panel_config.frame} to {base_frame}")

        # Update rotation laws with new frame
        rotation_laws[panel_config.frame] = trots.spice(
            base_frame=base_frame, target_frame=panel_config.frame
        )

    return rotation_laws


def create_paneled_model(
    vehicle_config: "VehicleSetup",
) -> "tvs.FullPanelledBodySettings":

    # Define geometry
    panel_collection = create_paneled_geometry_from_config(vehicle_config.shape)

    # Define reflection laws
    reflection_laws = create_paneled_reflection_laws_from_config(
        vehicle_config.radiation
    )

    # Define panel settings
    panel_settings: list[tvs.BodyPanelSettings] = []
    for panel_id, geometry in panel_collection.items():

        # Ensure that reflection law is defined for current panel
        if panel_id not in reflection_laws:
            log.fatal(f"Reflection law not defined for {panel_id}")
            log.fatal(traceback.extract_stack()[-2])
            exit(1)

        # Define settings
        panel_settings.append(
            tvs.body_panel_settings(
                panel_geometry=geometry,
                panel_reflection_law=reflection_laws[panel_id],
                panel_type_id=panel_id,
            )
        )

    # Define rotation models for all panels
    rotation_laws = create_paneled_rotation_models_from_config(
        shape_config=vehicle_config.shape,
        base_frame=vehicle_config.rotation.target_frame,
    )

    # Return paneled model
    return tvs.full_panelled_body_settings(
        panel_settings=panel_settings,
        part_rotation_model_settings=rotation_laws,
    )


def paneled_mex_model(
    vehicle_config: "VehicleSetup",
) -> tvs.FullPanelledBodySettings:

    # Define panels for HGA
    hga_area = 2.270  # OnShape
    hga_frame = "MEX_HGA"
    hga_geometries = front_back_panel_geometry(hga_area, hga_frame, "z")
    # hga_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, 1.0]),
    #     area=hga_area,
    #     frame_orientation=hga_frame,
    # )
    # hga_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, -1.0]),
    #     area=hga_area,
    #     frame_orientation=hga_frame,
    # )

    # Define panels for solar array Y+
    sp_area = 6.351  # Onshape
    spp_frame = "MEX_SA+Y_ZERO"
    spn_frame = "MEX_SA-Y_ZERO"
    spp_geometry = front_back_panel_geometry(sp_area, spp_frame, "z")
    spn_geometry = front_back_panel_geometry(sp_area, spn_frame, "z")
    sp_geometries = list(spp_geometry) + list(spn_geometry)
    spp_zp_panel_geometry, spp_zn_panel_geometry = front_back_panel_geometry(
        sp_area, spp_frame, "z"
    )

    # Define panels for solar array Y-
    spn_zp_panel_geometry, spn_zn_panel_geometry = front_back_panel_geometry(
        sp_area, spn_frame, "z"
    )
    sp_geometries = (
        spp_zp_panel_geometry,
        spp_zn_panel_geometry,
        spn_zp_panel_geometry,
        spn_zn_panel_geometry,
    )
    # spn_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, 1.0]),
    #     area=sp_area,
    #     frame_orientation=spn_frame,
    # )
    # spn_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, -1.0]),
    #     area=sp_area,
    #     frame_orientation=spn_frame,
    # )

    # Define panels for bus
    bus_frame = "MEX_SPACECRAFT"
    bus_xy_area = 2.434
    bus_xz_area = 2.006
    bus_yz_area = 2.398
    bus_xy_geometry = front_back_panel_geometry(bus_xy_area, bus_frame, "z")
    bus_xz_geometry = front_back_panel_geometry(bus_xz_area, bus_frame, "y")
    bus_yz_geometry = front_back_panel_geometry(bus_yz_area, bus_frame, "x")
    bus_geometries = (
        list(bus_xy_geometry) + list(bus_xz_geometry) + list(bus_yz_geometry)
    )
    # bus_xy_zp_panel_geometry, bus_xy_zn_panel_geometry = (
    #     front_back_panel_geometry(bus_xy_area, bus_frame, "z")
    # )
    # bus_xz_yp_panel_geometry, bus_xz_yn_panel_geometry = (
    #     front_back_panel_geometry(bus_xz_area, bus_frame, "y")
    # )
    # bus_yz_xp_panel_geometry, bus_yz_xn_panel_geometry = (
    #     front_back_panel_geometry(bus_yz_area, bus_frame, "x")
    # )
    # bus_geometries = (
    #     bus_xy_zp_panel_geometry,
    #     bus_xy_zn_panel_geometry,
    #     bus_xz_yp_panel_geometry,
    #     bus_xz_yn_panel_geometry,
    #     bus_yz_xp_panel_geometry,
    #     bus_yz_xn_panel_geometry,
    # )

    # bus_xy_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, 1.0]),
    #     area=bus_xy_area,
    #     frame_orientation=bus_frame,
    # )
    # bus_xy_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 0.0, -1.0]),
    #     area=bus_xy_area,
    #     frame_orientation=bus_frame,
    # )

    # bus_xz_yp_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, 1.0, 0.0]),
    #     area=bus_xz_area,
    #     frame_orientation=bus_frame,
    # )
    # bus_xz_yn_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([0.0, -1.0, 0.0]),
    #     area=bus_xz_area,
    #     frame_orientation=bus_frame,
    # )

    # bus_yz_xp_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([1.0, 0.0, 0.0]),
    #     area=bus_yz_area,
    #     frame_orientation=bus_frame,
    # )
    # bus_yz_xn_panel_geometry = tvs.frame_fixed_panel_geometry(
    #     surface_normal=np.array([-1.0, 0.0, 0.0]),
    #     area=bus_yz_area,
    #     frame_orientation=bus_frame,
    # )

    # Define reflection laws (Check with Dominic)
    sp_reflection_law = trad.specular_diffuse_body_panel_reflection(
        specular_reflectivity=0.05,
        diffuse_reflectivity=0.05,
        with_instantaneous_reradiation=False,
    )
    bus_reflection_law = trad.specular_diffuse_body_panel_reflection(
        specular_reflectivity=0.1,
        diffuse_reflectivity=0.3,
        with_instantaneous_reradiation=False,
    )
    hga_reflection_law = trad.specular_diffuse_body_panel_reflection(
        specular_reflectivity=0.3,
        diffuse_reflectivity=0.2,
        with_instantaneous_reradiation=False,
    )

    # sp_reflection_law = trad.lambertian_body_panel_reflection(1 - 0.72)
    # bus_reflection_law = trad.lambertian_body_panel_reflection(0.9)
    # hga_reflection_law = trad.lambertian_body_panel_reflection(0.9)

    # Define panels for solar arrays
    sp_panels = [
        tvs.body_panel_settings(
            panel_geometry=geometry,
            panel_reflection_law=sp_reflection_law,
        )
        for geometry in sp_geometries
    ]

    # spp_zp_panel = tvs.body_panel_settings(
    #     panel_geometry=spp_zp_panel_geometry,
    #     panel_reflection_law=sp_reflection_law,
    # )
    # spp_zn_panel = tvs.body_panel_settings(
    #     panel_geometry=spp_zn_panel_geometry,
    #     panel_reflection_law=sp_reflection_law,
    # )
    # spn_zp_panel = tvs.body_panel_settings(
    #     panel_geometry=spn_zp_panel_geometry,
    #     panel_reflection_law=sp_reflection_law,
    # )
    # spn_zn_panel = tvs.body_panel_settings(
    #     panel_geometry=spn_zn_panel_geometry,
    #     panel_reflection_law=sp_reflection_law,
    # )

    # Define panels for bus
    bus_panels = [
        tvs.body_panel_settings(
            panel_geometry=geometry,
            panel_reflection_law=bus_reflection_law,
        )
        for geometry in bus_geometries
    ]
    # bus_xy_zp_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_xy_zp_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )
    # bus_xy_zn_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_xy_zn_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )
    # bus_xz_yp_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_xz_yp_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )
    # bus_xz_yn_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_xz_yn_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )
    # bus_yz_xp_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_yz_xp_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )
    # bus_yz_xn_panel = tvs.body_panel_settings(
    #     panel_geometry=bus_yz_xn_panel_geometry,
    #     panel_reflection_law=bus_reflection_law,
    # )

    # Define panels for HGA
    hga_panels = [
        tvs.body_panel_settings(
            panel_geometry=geometry,
            panel_reflection_law=hga_reflection_law,
        )
        for geometry in hga_geometries
    ]
    # hga_zp_panel = tvs.body_panel_settings(
    #     panel_geometry=hga_zp_panel_geometry,
    #     panel_reflection_law=hga_reflection_law,
    # )
    # hga_zn_panel = tvs.body_panel_settings(
    #     panel_geometry=hga_zn_panel_geometry,
    #     panel_reflection_law=hga_reflection_law,
    # )

    # Define panelled model
    panels = hga_panels + sp_panels + bus_panels
    # panels = [
    #     hga_zp_panel,
    #     hga_zn_panel,
    #     spp_zp_panel,
    #     spp_zn_panel,
    #     spn_zp_panel,
    #     spn_zn_panel,
    #     bus_xy_zp_panel,
    #     bus_xy_zn_panel,
    #     bus_xz_yp_panel,
    #     bus_xz_yn_panel,
    #     bus_yz_xp_panel,
    #     bus_yz_xn_panel,
    # ]

    # Define rotation models for HGA and solar arrays
    rotation_models = {
        frame: trots.spice(base_frame=bus_frame, target_frame=frame)
        for frame in (spp_frame, spn_frame, hga_frame)
    }
    mex_paneled_model = tvs.full_panelled_body_settings(
        panel_settings=panels, part_rotation_model_settings=rotation_models
    )

    return mex_paneled_model
