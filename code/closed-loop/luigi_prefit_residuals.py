######################### # IMPORTANT #############################################################################

# This example computes the residuals for the MEX IFMS, one at the time (uses observations_from_ifms_files function)
# In order to test this example, I am using a Phobos Flyby IFMS file missing the few last/first lines...
# The removed lines were classified as outliers, but they should be filtered with the proper tudat functionality,
# rather than manually (as done for now)
##################################################################################################################
import os
import csv

import numpy as np
from matplotlib import pyplot as plt
from astropy.time import Time
from pathlib import Path

from tudatpy.interface import spice
from tudatpy.astro import time_representation as ttime
from tudatpy.math import interpolators
from tudatpy.dynamics import (
    environment,
    environment_setup,
)
from tudatpy.estimation import (
    observations,
    observations_setup,
    observable_models_setup,
)
from tudatpy import data as tdata

from datetime import datetime
import matplotlib.dates as mdates
from prefit import paths as ppaths

# Set Folders Containing Relevant Files
mex_kernels_folder = str(ppaths.datadir / "kernels")
mex_fdets_folder = str(ppaths.datadir / "fdets/complete")
mex_ifms_folder = str(ppaths.datadir / "ifms/filtered")
mex_odf_folder = str(ppaths.datadir / "odf")
ifms_residuals_path = ppaths.outdir / "residuals_alfonso"
ifms_residuals_path.mkdir(exist_ok=True, parents=True)


def create_single_ifms_collection_from_file(
    ifms_file,
    bodies,
    spacecraft_name,
    transmitting_station_name,
    reception_band,
    transmission_band,
):

    # Loading IFMS file
    # print(f'IFMS file: {ifms_file}\n with transmitting station: {transmitting_station_name} will be loaded.')
    single_ifms_collection = (
        observations_setup.observations_wrapper.observations_from_ifms_files(
            [ifms_file],
            bodies,
            spacecraft_name,
            transmitting_station_name,
            reception_band,
            transmission_band,
        )
    )

    return single_ifms_collection


def save_residuals_to_csv(
    station_residuals: dict[str, list[tuple]],
) -> list[Path]:

    output_files: list[Path] = []

    for site_name, site_data in station_residuals.items():

        output_file = ifms_residuals_path / f"{site_name}_residuals.csv"

        with output_file.open("w", newline="") as buffer:

            writer = csv.writer(buffer)

            writer.writerow([f"# Station: {site_name}"])
            writer.writerow(["# Time | UTC Time | Residuals"])

            # Write the data rows
            for record in site_data:
                times, utc_times, concatenated_residuals, _, _ = record

                # Write each UTC time and residual in a separate row
                for time, utc_time, residual in zip(
                    times, utc_times, concatenated_residuals
                ):
                    writer.writerow(
                        [
                            time,
                            utc_time.strftime(
                                "%Y-%m-%d %H:%M:%S"
                            ),  # Convert datetime to string
                            residual,
                        ]
                    )

        output_files.append(output_file)

    return output_files


def main() -> dict[str, list[tuple]]:

    # General configuration: Time interval (UTC: Should we transform to TDB?)
    start = ttime.DateTime(2013, 12, 26)
    end = ttime.DateTime(2013, 12, 30)
    buffer_time = ttime.Time(86400.0)

    # General configuration: System of bodies
    bodies_to_create: list[str] = [
        "Earth",
        "Sun",
        "Mercury",
        "Venus",
        "Mars",
        "Jupiter",
        "Saturn",
        "Moon",
        "Phobos",
    ]
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"

    # General configuration: spacecraft
    spacecraft_name = "MEX"  # Set Spacecraft Name
    spacecraft_central_body = "Mars"  # Set Central Body (Mars)

    # Augment default time range with buffer time
    start_time = start.to_epoch_time_object() - buffer_time
    end_time = end.to_epoch_time_object() + buffer_time

    # Create default body settings for celestial bodies
    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies=bodies_to_create,
        initial_time=start_time,
        final_time=end_time,
        base_frame_origin=global_frame_origin,
        base_frame_orientation=global_frame_orientation,
    )

    # Update rotation settings for Earth: Use tudatpy IERS instead of SPICE
    # NOTE: How do I set this up?
    # NOTE: Why would I use J2000 as base frame instead of GCRF?
    # NOTE: What happens when the interpolation settings are none?
    cio_and_tdbtt_interp_settings = (
        interpolators.InterpolatorGenerationSettings(
            interpolator_settings=interpolators.cubic_spline_interpolation(),
            initial_time=start_time,
            final_time=end_time,
            time_step=ttime.Time(3600.0),
        )
    )
    eop_interp_settings = interpolators.InterpolatorGenerationSettings(
        interpolator_settings=interpolators.cubic_spline_interpolation(),
        initial_time=start_time,
        final_time=end_time,
        time_step=ttime.Time(60.0),
    )
    body_settings.get("Earth").rotation_model_settings = (
        environment_setup.rotation_model.gcrs_to_itrs(
            precession_nutation_theory=environment_setup.rotation_model.IAUConventions.iau_2006,
            base_frame=global_frame_orientation,
            cio_interpolation_settings=cio_and_tdbtt_interp_settings,
            tdb_to_tt_interpolation_settings=cio_and_tdbtt_interp_settings,
            short_term_eop_interpolation_settings=eop_interp_settings,
        )
    )

    # Update gravity field settings of Earth to use the same reference frame
    # as the rotation model
    earth_gfield_settings = body_settings.get("Earth").gravity_field_settings
    assert isinstance(
        earth_gfield_settings,
        environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings,
    )
    earth_gfield_settings.associated_reference_frame = body_settings.get(
        "Earth"
    ).rotation_model_settings.target_frame
    body_settings.get("Earth").gravity_field_settings = earth_gfield_settings

    # Vehicle properties
    body_settings.add_empty_settings(spacecraft_name)

    # Retrieve translational ephemeris from SPICE
    # body_settings.get(spacecraft_name).ephemeris_settings = environment_setup.ephemeris.direct_spice(
    #    'SSB', 'J2000', 'MEX')
    # body_settings.get(spacecraft_name).ephemeris_settings = (
    #     environment_setup.ephemeris.direct_spice(
    #         frame_origin="Mars",
    #         frame_orientation="J2000",
    #         body_name_to_use="MEX",
    #     )
    # )
    body_settings.get(spacecraft_name).ephemeris_settings = (
        environment_setup.ephemeris.interpolated_spice(
            start_time,
            end_time,
            10.0,
            spacecraft_central_body,
            global_frame_orientation,
        )
    )

    # Retrieve rotational ephemeris from SPICE
    body_settings.get(spacecraft_name).rotation_model_settings = (
        environment_setup.rotation_model.spice(
            global_frame_orientation, spacecraft_name + "_SPACECRAFT", ""
        )
    )
    body_settings.get("Earth").ground_station_settings = (
        environment_setup.ground_station.radio_telescope_stations()
    )

    # Create System of Bodies using the above-defined body_settings
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Set default DSN turnaround ratio for the transponder of MEX
    bodies.get_body(
        spacecraft_name
    ).system_models.set_default_transponder_turnaround_ratio_function()

    # Define frequency bands for transmission and reception
    reception_band = observations_setup.ancillary_settings.FrequencyBands.x_band
    transmission_band = (
        observations_setup.ancillary_settings.FrequencyBands.x_band
    )

    # # Get position of HGA with respect to LVI (Reference point for MEX)
    # __et = spice.convert_date_string_to_ephemeris_time(start.iso_string())
    # x_hga_lvi_mexframe = spice.get_body_cartesian_position_at_epoch(
    #     target_body_name="MEX_HGA",
    #     observer_body_name="MEX_SPACECRAFT",
    #     reference_frame_name="MEX_SPACECRAFT",
    #     aberration_corrections="NONE",
    #     ephemeris_time=__et,
    # ) + np.array([1.3, 0, 0])
    # bodies.get_body(spacecraft_name).system_models.set_reference_point(
    #     reference_point="HGA",
    #     location=x_hga_lvi_mexframe,
    #     frame_origin="MEX_SPACECRAFT",
    #     frame_orientation="MEX_SPACECRAFT",
    # )

    # Can probably be deleted
    single_ifms_collections_list = list()
    atmospheric_corrections_list = list()
    transmitting_stations_list = list()

    labels = set()
    ifms_station_residuals: dict[str, list[tuple]] = {}

    # Collect IFMS files per station
    # ifms_files_per_station: dict[str, list[Path]] = {}
    for ifms_file in Path(mex_ifms_folder).glob("*.TAB"):

        # Identify station for current file
        match ifms_file.name[1:3]:
            case "14":
                station_name = "DSS14"
            case "63":
                station_name = "DSS63"
            case "32":
                station_name = "NWNORCIA"
            case _:
                raise ValueError("Invalid station code")

        # if station_name not in ifms_files_per_station:
        #     ifms_files_per_station[station_name] = [ifms_file]
        # else:
        #     ifms_files_per_station[station_name].append(ifms_file)

        compressed_observations = observations_setup.observations_wrapper.observations_from_ifms_files(
            [str(ifms_file)],
            bodies,
            spacecraft_name,
            station_name,
            reception_band,
            transmission_band,
            apply_troposphere_correction=True,
        )

        # # Set the HGA as reference point for the observation collection
        # compressed_observations.set_reference_point(
        #     bodies,
        #     x_hga_lvi_mexframe,
        #     "HGA",
        #     "MEX",
        #     observable_models_setup.links.LinkEndType.retransmitter,
        # )

        # Compress Doppler observations from 1.0 s integration time to 60.0 s
        # compressed_observations = estimation_setup.observation.create_compressed_doppler_collection(
        #    ifms_collection, 60, 10)
        # compressed_observations = ifms_collection

        # Generate UTC datetime objects for observation epochs
        times = compressed_observations.get_concatenated_observation_times()
        mjd_times = [ttime.seconds_since_epoch_to_julian_day(t) for t in times]
        utc_times = np.array(
            [
                Time(mjd_time, format="jd", scale="utc").datetime
                for mjd_time in mjd_times
            ]
        )

        #  Create light-time corrections list
        light_time_correction_list: list[
            observable_models_setup.light_time_corrections.LightTimeCorrectionSettings
        ] = [
            observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction(
                ["Sun"]
            ),
        ]

        doppler_link_ends = compressed_observations.link_definitions_per_observable[
            observable_models_setup.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
        ]

        ########## IMPORTANT STEP #######################################################################
        # When working with IFMS, Add: subtract_doppler_signature = False, or it won't work
        observation_model_settings: list[
            observable_models_setup.model_settings.ObservationModelSettings
        ] = []

        for current_link_definition in doppler_link_ends:

            print("PRINTING LINK ENDS -------------------------")
            # print(
            #     current_link_definition.link_end_id(
            #         observable_models_setup.links.LinkEndType.retransmitter
            #     ).reference_point
            # )
            # print(
            #     current_link_definition.link_end_id(
            #         observable_models_setup.links.LinkEndType.transmitter
            #     ).reference_point
            # )
            # print(
            #     current_link_definition.link_end_id(
            #         observable_models_setup.links.LinkEndType.receiver
            #     ).reference_point
            # )

            observation_model_settings.append(
                observable_models_setup.model_settings.dsn_n_way_doppler_averaged(
                    link_ends=current_link_definition,
                    light_time_correction_settings=light_time_correction_list,
                    subtract_doppler_signature=False,
                )
            )
            print("DONE PRINTING LINK ENDS -------------------------")

        ###################################################################################################
        # Create observation simulators.
        observation_simulators = observations_setup.observations_simulation_settings.create_observation_simulators(
            observation_settings=observation_model_settings,
            bodies=bodies,
        )

        # Add elevation and SEP angles dependent variables to the IFMS observation collection
        elevation_angle_settings = observations_setup.observations_dependent_variables.elevation_angle_dependent_variable(
            link_end_type=observable_models_setup.links.LinkEndType.receiver
        )
        elevation_angle_parser = compressed_observations.add_dependent_variable(
            elevation_angle_settings, bodies
        )
        sep_angle_settings = observations_setup.observations_dependent_variables.avoidance_angle_dependent_variable(
            body_name="Sun",
            start_link_end_type=observable_models_setup.links.LinkEndType.retransmitter,
            end_link_end_type=observable_models_setup.links.LinkEndType.receiver,
        )
        sep_angle_parser = compressed_observations.add_dependent_variable(
            sep_angle_settings, bodies
        )

        # Compute and set residuals in the IFMS observation collection
        observations.compute_residuals_and_dependent_variables(
            observation_collection=compressed_observations,
            observation_simulators=observation_simulators,
            bodies=bodies,
        )
        concatenated_obs = (
            compressed_observations.get_concatenated_observations()
        )
        concatenated_computed_obs = (
            compressed_observations.get_concatenated_computed_observations()
        )
        # Retrieve RMS and mean of the residuals
        concatenated_residuals = (
            compressed_observations.get_concatenated_residuals()
        )
        rms_residuals = compressed_observations.get_rms_residuals()
        mean_residuals = compressed_observations.get_mean_residuals()

        # print(f'Residuals: {concatenated_residuals}')
        print(f"Mean Residuals: {mean_residuals}")
        print(f"RMS Residuals: {rms_residuals}")

        # Populate Station Residuals Dictionary
        site_name = station_name
        if site_name not in ifms_station_residuals.keys():
            ifms_station_residuals[site_name] = [
                (
                    times,
                    utc_times,
                    concatenated_residuals,
                    mean_residuals,
                    rms_residuals,
                )
            ]
        else:
            ifms_station_residuals[site_name].append(
                (
                    times,
                    utc_times,
                    concatenated_residuals,
                    mean_residuals,
                    rms_residuals,
                )
            )

        print("End of the first loop")

    return ifms_station_residuals


if __name__ == "__main__":

    metak_path = ppaths.datadir / "metak_mex.tm"

    try:
        # Load Required Spice Kernels
        spice.load_standard_kernels()
        for kernel in os.listdir(mex_kernels_folder):
            kernel_path = os.path.join(mex_kernels_folder, kernel)
            spice.load_kernel(kernel_path)

        # loaded_kernels = putils.get_loaded_kernels("all")
        # for k in loaded_kernels:
        #     print(k.name)
        # exit(0)

        # Run main function
        ifms_station_residuals = main()

        residuals_files = save_residuals_to_csv(ifms_station_residuals)

    finally:
        spice.clear_kernels()
