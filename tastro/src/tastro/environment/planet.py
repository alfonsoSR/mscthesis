from ..config.environment.planets import PlanetSetup
from tudatpy.dynamics.environment_setup import (
    ephemeris as tephs,
    rotation_model as trots,
    shape as tshapes,
    gravity_field as tgravs,
    radiation_pressure as trad,
)
from tudatpy.math import interpolators as tint
from tudatpy.astro import time_representation as ttime
from tudatpy import data as tdata
from pathlib import Path
from ..core import SettingsGenerator
from ..logging import log


class PlanetSettings(SettingsGenerator[PlanetSetup]):

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

            case "horizons":

                log.debug("Horizons ephemerides")

                # Get horizons code from name
                if len(self.name.split(" ")) != 2:
                    raise ValueError(
                        "The name of the planet must include the MPC code to "
                        f"be able to set up horizons ephemerides: {self.name}"
                    )
                horizons_code = int(self.name.split(" ")[0])

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

                query = tephs.HorizonsQuery(
                    query_id=f"{horizons_code};",
                    location=f"@{frame_origin}",
                    epoch_start=(
                        self.config.time.initial_epoch - interpolation_buffer
                    ).to_float(),
                    epoch_end=(
                        self.config.time.final_epoch + interpolation_buffer
                    ).to_float(),
                    epoch_step=f"{int(interpolation_step.to_float()/60.)}m",
                    extended_query=True,
                )

                return query.create_ephemeris_tabulated(
                    frame_origin=frame_origin,
                    frame_orientation=frame_orientation,
                    aberations="geometric",
                )

            case _:

                raise ValueError(
                    f"Invalid ephemerides model: {self.local.ephemerides.model}"
                )

    def rotation_settings(self) -> trots.RotationModelSettings:

        # Define base and target frames
        base_frame = self.local.rotation.base_frame
        if base_frame == "global":
            base_frame = (
                self.config.environment.general.global_frame_orientation
            )
        target_frame = self.local.rotation.target_frame

        match self.local.rotation.model:

            case "spice":

                log.debug("Spice rotation settings")
                # Ensure target frame is set
                if target_frame is None:
                    raise ValueError(
                        f"Target frame not set with spice rotation model: {self.name}"
                    )

                return trots.spice(
                    base_frame=base_frame,
                    target_frame=target_frame,
                )

            case "precise":

                log.debug("Precise rotation settings")

                match self.name:

                    case "Mars":

                        return trots.mars_high_accuracy(
                            base_frame=base_frame,
                        )

                    case "Earth":

                        # TODO: Move to configuration?
                        cio_tdb_interpolation_step = ttime.Time(86400.0)
                        eop_interpolation_step = ttime.Time(60.0)

                        # Define buffers based on interpolation step
                        cio_tdb_buffer = cio_tdb_interpolation_step * 4
                        eop_buffer = eop_interpolation_step * 4

                        # Define settings for interpolators
                        cio_and_tdb_interp_settings = tint.InterpolatorGenerationSettings(
                            interpolator_settings=tint.cubic_spline_interpolation(),
                            initial_time=self.config.time.initial_epoch
                            - cio_tdb_buffer,
                            final_time=self.config.time.final_epoch
                            + cio_tdb_buffer,
                            time_step=cio_tdb_interpolation_step,
                        )
                        eop_interp_settings = tint.InterpolatorGenerationSettings(
                            interpolator_settings=tint.cubic_spline_interpolation(),
                            initial_time=self.config.time.initial_epoch
                            - eop_buffer,
                            final_time=self.config.time.final_epoch
                            + eop_buffer,
                            time_step=eop_interpolation_step,
                        )

                        # Return IAU 2000/2006 rotation model settings
                        return trots.gcrs_to_itrs(
                            precession_nutation_theory=trots.IAUConventions.iau_2006,
                            base_frame=base_frame,
                            cio_interpolation_settings=cio_and_tdb_interp_settings,
                            tdb_to_tt_interpolation_settings=cio_and_tdb_interp_settings,
                            short_term_eop_interpolation_settings=eop_interp_settings,
                        )

                    case _:

                        raise NotImplementedError(
                            f"Precise rotation model not set for {self.name}"
                        )

            case _:
                raise ValueError(
                    f"Invalid rotation model: {self.local.rotation.model}"
                )

    def shape_settings(self) -> tshapes.BodyShapeSettings:

        match self.local.shape.model:

            case "spherical_spice":
                log.debug("Spherical spice shape settings")
                return tshapes.spherical_spice()

            case "spherical":
                log.debug("Spherical shape settings")

                # Ensure equatorial radius is set
                if self.local.shape.equatorial_radius is None:
                    raise ValueError(
                        f"Requested spherical shape model for {self.name}: "
                        "Missing equatorial radius"
                    )

                return tshapes.spherical(
                    radius=self.local.shape.equatorial_radius
                )

            case "ellipsoid":
                log.debug("Ellipsoid shape settings")

                # Ensure equatorial radius is set
                if self.local.shape.equatorial_radius is None:
                    raise ValueError(
                        f"Requested ellipsoid shape model for {self.name}: "
                        "Missing equatorial radius"
                    )

                # Ensure flattening factor is set
                if self.local.shape.flattening_factor is None:
                    raise ValueError(
                        f"Requested ellipsoid shape model for {self.name}: "
                        "Missing flattening factor"
                    )

                return tshapes.oblate_spherical(
                    equatorial_radius=self.local.shape.equatorial_radius,
                    flattening=(1.0 / self.local.shape.flattening_factor),
                )

            case _:
                raise NotImplementedError(
                    f"Invalid shape model: {self.local.shape.model}"
                )

    def gravity_settings(self) -> tgravs.GravityFieldSettings:

        match self.local.gravity.model:

            case "spice_point_mass":
                log.debug("Point mass gravity settings")
                return tgravs.central_spice(self.name)

            case "spherical_harmonics":
                log.debug("Spherical harmonics gravity settings")

                # Extract information about gravity field model
                sh_model_file = self.local.gravity.spherical_harmonics_file
                if sh_model_file is None:
                    raise ValueError(
                        "Missing source file for spherical harmonics model "
                        f"of {self.name}"
                    )
                sh_degree = self.local.gravity.spherical_harmonics_degree
                if sh_degree is None:
                    raise ValueError(
                        "Missing max degree for spherical harmonics model "
                        f"of {self.name}"
                    )
                sh_order = self.local.gravity.spherical_harmonics_order
                if sh_order is None:
                    raise ValueError(
                        "Missing max order for spherical harmonics model "
                        f"of {self.name}"
                    )
                sh_frame = self.local.gravity.spherical_harmonics_frame
                if sh_frame is None:
                    raise ValueError(
                        "Missing frame for spherical harmonics model "
                        f"of {self.name}"
                    )

                # Get path to spherical harmonics gravity field model
                sh_source = Path(tdata.get_gravity_models_path()).resolve()
                sh_source /= f"{self.name}/{sh_model_file}"
                if not sh_source.exists():
                    raise FileNotFoundError(
                        f"Invalid spherical harmonics model: {sh_source}"
                    )

                return tgravs.from_file_spherical_harmonic(
                    file=str(sh_source),
                    maximum_degree=sh_degree,
                    maximum_order=sh_order,
                    associated_reference_frame=sh_frame,
                )

            case "central_sbdb":

                log.debug("SBDB point mass gravity settings")

                # Get horizons code from name
                if len(self.name.split(" ")) != 2:
                    raise ValueError(
                        "The name of the planet must include the MPC code to "
                        f"be able to set up horizons ephemerides: {self.name}"
                    )
                horizons_code = int(self.name.split(" ")[0])

                return tgravs.central_sbdb(horizons_code)

            case _:
                raise NotImplementedError(
                    f"Invalid gravity model: {self.local.gravity.model}"
                )

    def radiation_source_settings(self) -> trad.RadiationSourceModelSettings:

        match self.local.radiation.model:

            case "direct_isotropic":

                log.debug("Direct isotropic radiation source")

                # Define settings for luminosity model
                luminosity = self.local.radiation.direct_setup.luminosity
                if luminosity is None:
                    raise ValueError(
                        "Missing luminosity for direct radiation source"
                    )
                luminosity_settings = trad.constant_luminosity(luminosity)

                # Return isotropic source settings
                return trad.isotropic_radiation_source(luminosity_settings)

            case _:
                raise NotImplementedError(
                    f"Invalid radiation source model: {self.local.radiation.model}"
                )
