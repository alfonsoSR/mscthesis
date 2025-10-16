from ..config.environment.vehicles import VehicleSetup
from tudatpy.dynamics.environment_setup import (
    ephemeris as tephs,
    rotation_model as trots,
    radiation_pressure as trad,
    shape as tshape,
    vehicle_systems as tvs,
)
import numpy as np
from ..core import SettingsGenerator
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.interface import spice
from ..logging import log
from ..config import CaseSetup


def paneled_mex_model() -> tvs.FullPanelledBodySettings:

    # Define panels for HGA
    hga_area = 2.270  # OnShape
    hga_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, 1.0]),
        area=hga_area,
        frame_orientation="MEX_HGA",
    )
    hga_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, -1.0]),
        area=hga_area,
        frame_orientation="MEX_HGA",
    )

    # Define panels for solar array Y+
    sp_area = 6.351  # Onshape
    spp_frame = "MEX_SPACECRAFT"
    spp_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, 1.0]),
        area=sp_area,
        frame_orientation=spp_frame,
    )
    spp_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, -1.0]),
        area=sp_area,
        frame_orientation=spp_frame,
    )

    # Define panels for solar array Y-
    spn_frame = "MEX_SPACECRAFT"
    spn_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, 1.0]),
        area=sp_area,
        frame_orientation=spn_frame,
    )
    spn_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, -1.0]),
        area=sp_area,
        frame_orientation=spn_frame,
    )

    # Define panels for bus
    bus_frame = "MEX_SPACECRAFT"
    bus_xy_area = 2.434
    bus_xz_area = 2.006
    bus_yz_area = 2.398

    bus_xy_zp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, 1.0]),
        area=bus_xy_area,
        frame_orientation=bus_frame,
    )
    bus_xy_zn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 0.0, -1.0]),
        area=bus_xy_area,
        frame_orientation=bus_frame,
    )

    bus_xz_yp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, 1.0, 0.0]),
        area=bus_xz_area,
        frame_orientation=bus_frame,
    )
    bus_xz_yn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([0.0, -1.0, 0.0]),
        area=bus_xz_area,
        frame_orientation=bus_frame,
    )

    bus_yz_xp_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([1.0, 0.0, 0.0]),
        area=bus_yz_area,
        frame_orientation=bus_frame,
    )
    bus_yz_xn_panel_geometry = tvs.frame_fixed_panel_geometry(
        surface_normal=np.array([-1.0, 0.0, 0.0]),
        area=bus_yz_area,
        frame_orientation=bus_frame,
    )

    # Define reflection laws (Check with Dominic)
    sp_reflection_law = trad.lambertian_body_panel_reflection(1 - 0.72)
    bus_reflection_law = trad.lambertian_body_panel_reflection(0.9)
    hga_reflection_law = trad.lambertian_body_panel_reflection(0.9)

    # Define panels for solar arrays
    spp_zp_panel = tvs.body_panel_settings(
        panel_geometry=spp_zp_panel_geometry,
        panel_reflection_law=sp_reflection_law,
    )
    spp_zn_panel = tvs.body_panel_settings(
        panel_geometry=spp_zn_panel_geometry,
        panel_reflection_law=sp_reflection_law,
    )
    spn_zp_panel = tvs.body_panel_settings(
        panel_geometry=spn_zp_panel_geometry,
        panel_reflection_law=sp_reflection_law,
    )
    spn_zn_panel = tvs.body_panel_settings(
        panel_geometry=spn_zn_panel_geometry,
        panel_reflection_law=sp_reflection_law,
    )

    # Define panels for bus
    bus_xy_zp_panel = tvs.body_panel_settings(
        panel_geometry=bus_xy_zp_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )
    bus_xy_zn_panel = tvs.body_panel_settings(
        panel_geometry=bus_xy_zn_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )
    bus_xz_yp_panel = tvs.body_panel_settings(
        panel_geometry=bus_xz_yp_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )
    bus_xz_yn_panel = tvs.body_panel_settings(
        panel_geometry=bus_xz_yn_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )
    bus_yz_xp_panel = tvs.body_panel_settings(
        panel_geometry=bus_yz_xp_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )
    bus_yz_xn_panel = tvs.body_panel_settings(
        panel_geometry=bus_yz_xn_panel_geometry,
        panel_reflection_law=bus_reflection_law,
    )

    # Define panels for HGA
    hga_zp_panel = tvs.body_panel_settings(
        panel_geometry=hga_zp_panel_geometry,
        panel_reflection_law=hga_reflection_law,
    )
    hga_zn_panel = tvs.body_panel_settings(
        panel_geometry=hga_zn_panel_geometry,
        panel_reflection_law=hga_reflection_law,
    )

    # Define panelled model
    panels = [
        # hga_zp_panel,
        # hga_zn_panel,
        spp_zp_panel,
        spp_zn_panel,
        spn_zp_panel,
        spn_zn_panel,
        bus_xy_zp_panel,
        bus_xy_zn_panel,
        bus_xz_yp_panel,
        bus_xz_yn_panel,
        bus_yz_xp_panel,
        bus_yz_xn_panel,
    ]
    mex_paneled_model = tvs.full_panelled_body_settings(panels)

    return mex_paneled_model


class VehicleSettings(SettingsGenerator[VehicleSetup]):

    def ephemeris_settings(self) -> tephs.EphemerisSettings:

        # Define frame origin
        frame_origin = self.local.ephemerides.ephemeris_frame_origin
        if frame_origin == "global":
            frame_origin = self.config.environment.general.global_frame_origin

        # Define frame orientation
        frame_orientation = self.local.ephemerides.ephemeris_frame_orientation
        if frame_orientation == "global":
            frame_orientation = (
                self.config.environment.general.global_frame_orientation
            )

        match self.local.ephemerides.model:

            case "direct":

                log.debug("Direct ephemerides")
                return tephs.direct_spice(
                    frame_origin=frame_origin,
                    frame_orientation=frame_orientation,
                )

            case "interpolated":

                log.debug("Interpolated ephemerides")
                # Ensure interpolation step is defined
                interpolation_step = self.local.ephemerides.interpolation_step
                if interpolation_step is None:
                    raise ValueError("Interpolation step is not defined")

                # Ensure interpolation buffer is defined
                interpolation_buffer = (
                    self.local.ephemerides.interpolation_buffer
                )
                if interpolation_buffer is None:
                    raise ValueError("Interpolation buffer not defined")

                return tephs.interpolated_spice(
                    frame_origin=frame_origin,
                    frame_orientation=frame_orientation,
                    initial_time=self.config.time.initial_epoch
                    - interpolation_buffer,
                    final_time=self.config.time.final_epoch
                    + interpolation_buffer,
                    time_step=interpolation_step,
                )

            case _:

                raise ValueError(
                    f"Invalid ephemerides model: {self.local.ephemerides.model}"
                )

    def rotation_settings(self) -> trots.RotationModelSettings:

        # Define base frame
        base_frame = self.local.rotation.base_frame
        if base_frame == "global":
            base_frame = (
                self.config.environment.general.global_frame_orientation
            )

        # # Define target frame
        # target_frame = self.local.rotation.target_frame
        # if target_frame == "global":
        #     target_frame = (
        #         self.config.environment.general.global_frame_orientation
        #     )

        match self.local.rotation.model:

            case "spice":

                log.debug("Spice rotation settings")
                return trots.spice(
                    base_frame=base_frame,
                    target_frame=self.local.rotation.target_frame,
                )

            case "precise":

                log.debug("Precise rotation settings")
                raise NotImplementedError(
                    f"Precise rotation model not implemented for {self.name}"
                )

            case _:
                raise ValueError(
                    f"Invalid rotation model: {self.local.rotation.model}"
                )

    def radiation_target_settings(
        self,
    ) -> trad.RadiationPressureTargetModelSettings:

        # Get dictionary of occulting bodies per source
        acceleration_setup = self.config.propagation.accelerations[self.name]
        per_source_occulting_bodies: dict[str, list[str]] = {}
        for planet, planet_setup in acceleration_setup.external.items():

            # Skip if source does not exert radiation pressure
            if not planet_setup.radiation.present:
                continue

            # Skip if source does not have occulting bodies
            if planet_setup.radiation.occulting_bodies is None:
                continue

            # Check that source has shape settings
            if not self.config.environment.planets[planet].shape.present:
                log.error(
                    f"Shape settings required to use {planet} as radiation "
                    "source with occulting bodies"
                )
                exit(1)

            # Check that all the occulting bodies have shape settings
            for occulting_body in planet_setup.radiation.occulting_bodies:

                if not self.config.environment.planets[
                    occulting_body
                ].shape.present:
                    log.error(
                        f"Shape settings required to use {occulting_body}"
                        " as occulting body"
                    )
                    exit(1)

            # Set occulting bodies for source
            per_source_occulting_bodies[planet] = (
                planet_setup.radiation.occulting_bodies
            )

        match self.local.radiation.model:

            case "cannonball":

                log.debug("Cannonball radiation pressure")

                # Load setup for cannonball radiation model
                cannonball_setup = self.local.radiation.cannonball_settings
                reference_area = cannonball_setup.reference_area
                if reference_area is None:
                    raise ValueError(
                        "Missing reference area for cannonball radiation"
                    )
                radiation_coefficient = (
                    cannonball_setup.radiation_pressure_coefficient
                )
                if radiation_coefficient is None:
                    raise ValueError(
                        "Missing radiation coefficient for cannonball radiation"
                    )

                # Return radiation target settings
                return trad.cannonball_radiation_target(
                    reference_area=reference_area,
                    radiation_pressure_coefficient=radiation_coefficient,
                    per_source_occulting_bodies=per_source_occulting_bodies,
                )

            case "paneled":

                log.debug("Paneled radiation target settings")

                if not (
                    self.local.shape.present
                    and (self.local.shape.model == "paneled")
                ):
                    raise ValueError(
                        "Requested paneled radiation without specifying shape settings"
                    )

                # Add self-shadowing option to configuration
                return trad.panelled_radiation_target(
                    source_to_target_occulting_bodies=per_source_occulting_bodies,
                    maximum_number_of_pixels_per_source={"Sun": 0},
                )

            case _:
                raise NotImplementedError(
                    f"Invalid radiation target model: {self.local.radiation.model}"
                )

    def shape_settings(self) -> tvs.FullPanelledBodySettings:

        match self.local.shape.model:

            case "paneled":

                log.debug("Paneled shape settings")

                if self.name != "MEX":
                    raise ValueError(
                        "Paneled shape model only available for MEX"
                    )

                return paneled_mex_model()

            case _:
                raise NotImplementedError(
                    f"Invalid vehicle shape model: {self.local.shape.model}"
                )

    def doppler_tracking_settings(
        self, bodies: "SystemOfBodies"
    ) -> "SystemOfBodies":

        raise PendingDeprecationWarning(
            "Use ClosedLoopSettingsGenerator.update_vehicles instead"
        )

        for vehicle, vehicle_setup in self.config.environment.vehicles.items():

            # Skip if Doppler settings are not present
            if not vehicle_setup.doppler.present:
                continue

            if vehicle != "MEX":
                raise NotImplementedError("Only MEX supported")

            # Set turnaround ratio
            match vehicle_setup.doppler.turnaround_ratio:

                case "default":
                    bodies.get(
                        vehicle
                    ).system_models.set_default_transponder_turnaround_ratio_function()
                case _:
                    raise NotImplementedError("Invalid turnaround ratio setup")

            # Set reference point for tracking
            match vehicle_setup.doppler.reference_point:

                case "origin":
                    pass
                case "HGA":
                    # Get state of HGA wrt LVI in fixed frame
                    cstate_hga_lvi_fixed = (
                        spice.get_body_cartesian_state_at_epoch(
                            target_body_name="MEX_HGA",
                            observer_body_name="MEX_SPACECRAFT",
                            reference_frame_name="MEX_SPACECRAFT",
                            aberration_corrections="none",
                            ephemeris_time=self.config.time.initial_epoch,
                        )
                    )

                    # Define constant ephemerides
                    hga_ephemeris_settings = tephs.constant(
                        constant_state=cstate_hga_lvi_fixed,
                        frame_origin="MEX_SPACECRAFT",
                        frame_orientation="MEX_SPACECRAFT",
                    )
                    hga_ephemerides = tephs.create_ephemeris(
                        hga_ephemeris_settings, "HGA"
                    )

                    # Update MEX with reference point
                    bodies.get(vehicle).system_models.set_reference_point(
                        reference_point="HGA",
                        ephemeris=hga_ephemerides,
                    )

                case _:
                    raise NotImplementedError(
                        "Invalid reference point settings"
                    )

        return bodies
