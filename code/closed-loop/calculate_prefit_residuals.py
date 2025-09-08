from tudatpy.astro import (
    time_representation as ttime,
    element_conversion as telc,
)
from tudatpy.dynamics import environment_setup as tenvs, environment as tenv
from prefit import utils as putils, io as pio, paths as ppaths
from tudatpy.interface import spice
from pathlib import Path
from tudatpy.estimation import (
    observations_setup as tobss,
    observable_models_setup as tomss,
    observations as tobs,
)
import numpy as np
from tudatpy.math import interpolators as tinterp
from tudatpy.dynamics.environment_setup.ground_station import (
    get_station_reference_state_itrf2000,
)


def define_system_of_bodies(spacecraft_name: str) -> tenv.SystemOfBodies:

    # General configuration: Time interval in UTC
    initial_epoch_utc = ttime.DateTime.from_iso_string("2013-12-26T00:00:00")
    final_epoch_utc = ttime.DateTime.from_iso_string("2013-12-30T00:00:00")
    buffer_time = ttime.Time(86400.0)

    # General configuration: System of bodies
    global_frame_origin: str = "SSB"
    global_frame_orientation: str = "J2000"
    central_body: str = "Mars"

    # Apply buffer to initial and final epochs (UTC)
    initial_epoch_buffer = (
        initial_epoch_utc.to_epoch_time_object() - buffer_time
    )
    final_epoch_buffer = final_epoch_utc.to_epoch_time_object() + buffer_time

    # Transform epochs to TDB
    time_converter = ttime.default_time_scale_converter()
    initial_epoch_buffer = time_converter.convert_time_object(
        ttime.TimeScales.utc_scale,
        ttime.TimeScales.tdb_scale,
        initial_epoch_buffer,
    )
    final_epoch_buffer = time_converter.convert_time_object(
        ttime.TimeScales.utc_scale,
        ttime.TimeScales.tdb_scale,
        final_epoch_buffer,
    )

    # Initialize environment_settings
    environment_settings = tenvs.BodyListSettings(
        frame_orientation=global_frame_orientation,
        frame_origin=global_frame_origin,
    )
    # environment_settings = tenvs.get_default_body_settings_time_limited(
    #     bodies=["Sun", "Moon"],
    #     initial_time=initial_epoch_buffer,
    #     final_time=final_epoch_buffer,
    #     base_frame_orientation=global_frame_orientation,
    #     base_frame_origin=global_frame_origin,
    # )

    # Define settings for spacecraft
    #############################################
    environment_settings.add_empty_settings(spacecraft_name)
    mex_settings = environment_settings.get(spacecraft_name)

    # Define settings for translational ephemerides
    mex_settings.ephemeris_settings = tenvs.ephemeris.direct_spice(
        frame_orientation=global_frame_orientation,
        frame_origin=central_body,
        body_name_to_use=spacecraft_name,
    )

    # # Missing orientation settings, potentially shape model
    # mex_settings.rotation_model_settings = tenvs.rotation_model.spice(
    #     base_frame=global_frame_orientation,
    #     target_frame="MEX_SPACECRAFT",
    # )

    # Define settings for Mars
    ############################################
    environment_settings.add_empty_settings("Mars")
    mars_settings = environment_settings.get("Mars")

    # Define settings for translational ephemerides
    mars_settings.ephemeris_settings = tenvs.ephemeris.interpolated_spice(
        initial_time=initial_epoch_buffer,
        final_time=final_epoch_buffer,
        time_step=ttime.Time(3600.0),
        frame_orientation=global_frame_orientation,
        frame_origin=global_frame_origin,
        body_name_to_use="Mars",
    )

    # Update settings for the Earth
    ############################################
    environment_settings.add_empty_settings("Earth")
    earth_settings = environment_settings.get("Earth")

    # Define settings for rotation model
    cio_and_tdbtt_interp_settings = tinterp.InterpolatorGenerationSettings(
        interpolator_settings=tinterp.cubic_spline_interpolation(),
        initial_time=initial_epoch_buffer,
        final_time=final_epoch_buffer,
        time_step=ttime.Time(3600.0),
    )
    eop_interp_settings = tinterp.InterpolatorGenerationSettings(
        interpolator_settings=tinterp.cubic_spline_interpolation(),
        initial_time=initial_epoch_buffer,
        final_time=final_epoch_buffer,
        time_step=ttime.Time(60.0),
    )
    earth_settings.rotation_model_settings = tenvs.rotation_model.gcrs_to_itrs(
        precession_nutation_theory=tenvs.rotation_model.IAUConventions.iau_2006,
        base_frame=global_frame_orientation,
        cio_interpolation_settings=cio_and_tdbtt_interp_settings,
        tdb_to_tt_interpolation_settings=cio_and_tdbtt_interp_settings,
        short_term_eop_interpolation_settings=eop_interp_settings,
    )

    # Define settings for translational ephemerides
    earth_settings.ephemeris_settings = tenvs.ephemeris.interpolated_spice(
        initial_time=initial_epoch_buffer,
        final_time=final_epoch_buffer,
        time_step=ttime.Time(3600.0),
        frame_origin=global_frame_origin,
        frame_orientation=global_frame_orientation,
        body_name_to_use="Earth",
    )

    # # Define shape deformation settings (Displacements of ground stations)
    # earth_settings.shape_deformation_settings = [
    #     tenvs.shape_deformation.basic_solid_body_tidal(
    #         tide_raising_bodies=["Sun", "Moon"],
    #         displacement_love_numbers={2: (0.6, 0.08)},
    #     )
    # ]

    # # Refer gravity field of the Earth to the same frame as its rotation model
    # __gfield_settings = environment_settings.get("Earth").gravity_field_settings
    # assert isinstance(
    #     __gfield_settings,
    #     tenvs.gravity_field.SphericalHarmonicsGravityFieldSettings,
    # )
    # __gfield_settings.associated_reference_frame = environment_settings.get(
    #     "Earth"
    # ).rotation_model_settings.target_frame
    # environment_settings.get("Earth").gravity_field_settings = __gfield_settings

    # Define settings for ground stations (TODO: Do this manually)
    stations_to_create = ["DSS63", "NWNORCIA"]
    ground_station_settings: list[
        tenvs.ground_station.GroundStationSettings
    ] = []
    for station_name in stations_to_create:

        # Reference epoch for ground station motion
        station_ref_epoch = ttime.DateTime.from_iso_string(
            "2000-01-01T00:00:00"
        )
        ref_state_itrf2000 = get_station_reference_state_itrf2000(station_name)

        # Transform reference state to ITRF2014 (Consistent with our EOPs)
        ref_state_itrf2014 = (
            tenvs.convert_ground_station_state_between_itrf_frames(
                ground_station_state=ref_state_itrf2000,
                epoch=station_ref_epoch.to_epoch_time_object(),
                base_frame="ITRF2000",
                target_frame="ITRF2014",
            )
        )

        # Define motion settings for ground station
        station_motion_settings = [
            tenvs.ground_station.LinearGroundStationMotionSettings(
                linear_velocity=ref_state_itrf2014[3:],
                reference_epoch=station_ref_epoch.to_epoch_time_object(),
            ),
            # tenvs.ground_station.BodyDeformationStationMotionSettings(),
        ]

        # Create object with station settings
        station_settings = tenvs.ground_station.basic_station(
            station_name=station_name,
            station_nominal_position=ref_state_itrf2014[:3],
            station_position_element_type=telc.PositionElementTypes.cartesian_position_type,
            station_motion_settings=station_motion_settings,
        )
        ground_station_settings.append(station_settings)

    earth_settings.ground_station_settings = ground_station_settings

    # shape_deformation_settings = [
    #     tenvs.shape_deformation.basic_solid_body_tidal(
    #         tide_raising_bodies=["Sun", "Moon"],
    #         displacement_love_numbers={2: (0.6, 0.08)},
    #     )
    # ]
    # environment_settings.get("Earth").shape_deformation_settings = (
    #     shape_deformation_settings
    # )

    # exit(0)

    # sta = environment_settings.get("Earth").ground_station_settings[0]
    # for item in sta.station_motion_settings:

    #     match item.model_type:
    #         case tenvs.ground_station.StationMotionModelTypes.linear:
    #             assert isinstance(
    #                 item, tenvs.ground_station.LinearGroundStationMotionSettings
    #             )
    #             print("Has linear motion settings")
    #             pass
    #         case tenvs.ground_station.StationMotionModelTypes.body_deformation:
    #             print("Has body deformation motion settings")
    #             pass
    #         case _:
    #             print(f"Type not considered: {item.model_type}")

    # exit(0)

    # Create system of bodies
    bodies = tenvs.create_system_of_bodies(environment_settings)

    # Set up spacecraft to use default turnaround ratios
    bodies.get_body(
        spacecraft_name
    ).system_models.set_default_transponder_turnaround_ratio_function()

    return bodies


def group_ifms_files_per_station(source_dir: Path) -> dict[str, list[Path]]:

    # Get IFMS files in directory and group them per station
    ifms_per_station: dict[str, list[Path]] = {}
    for file in source_dir.glob("*.TAB"):

        # Get metadata from file name
        metadata = pio.get_metadata_from_ifms_filename(file)

        # Get station name from metadata
        match metadata["station_id"]:
            case "32":
                station_name = "NWNORCIA"
            case "63":
                station_name = "DSS63"
            case "14":
                station_name = "DSS14"
            case _:
                raise ValueError("Invalid station id")

        # Update collection with file
        if station_name not in ifms_per_station:
            ifms_per_station[station_name] = [file]
        else:
            ifms_per_station[station_name].append(file)

    return ifms_per_station


def define_observation_collection_for_station(
    bodies: tenv.SystemOfBodies, station: str, station_ifms: list[Path]
) -> tuple[pio.TwoWayDopplerObservations, tobs.ObservationCollection]:

    # Load data from IFMS files into Python object
    content = pio.load_doppler_observations_from_ifms_files(station_ifms)

    # Define interpolator for transmitter frequency
    interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
        start_times=content.ramping_tdb0[:-1],
        end_times=content.ramping_tdb0[1:],
        ramp_rates=content.ramping_df[:-1].tolist(),
        start_frequency=content.ramping_f0[:-1].tolist(),
    )
    station_obj = bodies.get("Earth").get_ground_station(station)
    station_obj.set_transmitting_frequency_calculator(interpolator)

    # Define ancilliary settings for observation collection
    ancilliary = tobss.ancillary_settings.dsn_n_way_doppler_ancilliary_settings(
        frequency_bands=[
            tobss.ancillary_settings.FrequencyBands.x_band,
            tobss.ancillary_settings.FrequencyBands.x_band,
        ],
        reference_frequency_band=tobss.ancillary_settings.FrequencyBands.x_band,
        reference_frequency=0.0,
        integration_time=ttime.Time(1.0),
    )

    # Define link ends
    transponder = tomss.links.body_origin_link_end_id(spacecraft_name)
    ground_station = tomss.links.body_reference_point_link_end_id(
        body_name="Earth",
        reference_point_id=station,
    )
    link_ends = {
        tomss.links.LinkEndType.transmitter: ground_station,
        tomss.links.LinkEndType.retransmitter: transponder,
        tomss.links.LinkEndType.receiver: ground_station,
    }
    link_definition = tomss.links.LinkDefinition(link_ends)

    # Create observations set
    obs_format = [np.array([x]) for x in content.observation_values]
    observation_set = tobs.single_observation_set(
        observable_type=tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type,
        link_definition=link_definition,
        observations=obs_format,
        observation_times=content.observation_epochs_et.tolist(),
        reference_link_end=tomss.links.LinkEndType.receiver,
        ancilliary_settings=ancilliary,
    )
    observations = tobs.ObservationCollection([observation_set])

    return content, observations


# def compute_prefit_residuals_for_ifms_file(
#     bodies: tenv.SystemOfBodies, ifms: Path, target_name: str
# ) -> tuple[str, tobs.ObservationCollection] | None:

#     # Extract metadata from file name
#     metadata = putils.get_metadata_from_ifms_filename(ifms)

#     # Identify station
#     subtract_signature: bool = False
#     match metadata["station_id"]:
#         case "14":
#             station_name = "DSS14"
#             return None
#         case "63":
#             station_name = "DSS63"
#             return None
#         case "32":
#             station_name = "NWNORCIA"
#         case "62":
#             station_name = "DSS-62"
#             return None
#         case _:
#             return None
#             raise ValueError("Invalid station code")

#     # Load raw contents from file and update with metadata
#     ifms_content = tdata.read_ifms_file(
#         str(ifms), apply_tropospheric_correction=False
#     )

#     # Define link ends
#     uplink_station = tomss.links.body_reference_point_link_end_id(
#         body_name="Earth",
#         reference_point_id=station_name,
#     )
#     downlink_station = tomss.links.body_reference_point_link_end_id(
#         body_name="Earth",
#         reference_point_id=station_name,
#     )
#     transponder = tomss.links.body_reference_point_link_end_id(
#         body_name=spacecraft_name,
#         reference_point_id="HGA",
#     )
#     link_ends = {
#         tomss.links.LinkEndType.transmitter: uplink_station,
#         tomss.links.LinkEndType.retransmitter: transponder,
#         tomss.links.LinkEndType.receiver: downlink_station,
#     }

#     # Get ET observation epochs from IFMS
#     if (
#         tdata.TrackingDataType.tdb_reception_time_j2000
#         not in ifms_content.get_available_datatypes()
#     ):
#         raise ValueError(f"IFMS file does not contain observation epochs in ET")

#     observation_epochs_et = np.array(
#         [
#             ttime.Time(float(ti))
#             for ti in ifms_content.raw_datamap["tdb_seconds_since_j2000"]
#         ]
#     )

#     # Correct observation epochs with delays (I think irrelevant)
#     if (
#         tdata.TrackingDataType.time_tag_delay
#         in ifms_content.get_available_datatypes()
#     ):
#         raise NotImplementedError("Handling of delays in IFMS not implemented")

#     # Get observations from IFMS
#     if (
#         tdata.TrackingDataType.doppler_averaged_frequency
#         not in ifms_content.get_available_datatypes()
#     ):
#         raise ValueError(f"IFMS file does not contain observations")

#     observation_values = np.array(
#         ifms_content.raw_datamap["doppler_averaged_frequency_hz"], dtype=float
#     )

#     # Define transmission and reception bands, and observation model
#     match metadata["data_type"]:
#         case "DPX":
#             tx_band = tobss.ancillary_settings.FrequencyBands.x_band
#             rx_band = tobss.ancillary_settings.FrequencyBands.x_band
#             obs_type = (
#                 tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
#             )
#         case "D2X":
#             tx_band = tobss.ancillary_settings.FrequencyBands.x_band
#             rx_band = tobss.ancillary_settings.FrequencyBands.x_band
#             obs_type = (
#                 tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
#             )
#         case _:
#             raise ValueError(f"Data type {metadata['data_type']} not supported")

#     # Load observations from file
#     observations = tobss.observations_wrapper.observations_from_ifms_files(
#         ifms_file_names=[str(ifms)],
#         bodies=bodies,
#         ground_station_name=station_name,
#         target_name=target_name,
#         reception_band=rx_band,
#         transmission_band=tx_band,
#         apply_troposphere_correction=False,
#         earth_fixed_station_positions={
#             station_name: bodies.get_body("Earth")
#             .get_ground_station(station_name)
#             .station_state.cartesian_positon_at_reference_epoch
#         },
#     )

#     # Retrieve link ends from observation collection
#     __link_ends = observations.link_definitions_per_observable[obs_type]
#     assert len(__link_ends) == 1
#     link_ends = __link_ends[0]

#     # Define light-time corrections
#     # light_time_settings = [
#     #     tomss.light_time_corrections.first_order_relativistic_light_time_correction(
#     #         ["Sun"]
#     #     )
#     # ]

#     # Define observation settings for current link end
#     observation_model_settings: list[
#         tomss.model_settings.ObservationModelSettings
#     ] = [
#         tomss.model_settings.dsn_n_way_doppler_averaged(
#             link_ends=link_ends,
#             subtract_doppler_signature=False,
#             # light_time_correction_settings=light_time_settings,
#         )
#     ]

#     # Create observation simulator
#     observation_simulator = (
#         tobss.observations_simulation_settings.create_observation_simulators(
#             observation_settings=observation_model_settings, bodies=bodies
#         )
#     )

#     # Compute pre-fit residuals
#     tobs.compute_residuals_and_dependent_variables(
#         observation_collection=observations,
#         observation_simulators=observation_simulator,
#         bodies=bodies,
#     )

#     # Return updated observation collection and station name
#     return station_name, observations


if __name__ == "__main__":

    # Paths
    metakernel = ppaths.datadir / "metak_mex.tm"
    ifms_dir = ppaths.datadir / "ifms/filtered"
    output_dir = ppaths.outdir / "prefit-closed-loop"
    output_dir.mkdir(exist_ok=True)

    try:
        spice.load_kernel(str(metakernel))

        # Define system of bodies
        spacecraft_name: str = "MEX"
        bodies = define_system_of_bodies(spacecraft_name)

        # Group IFMS files per station
        ifms_per_station = group_ifms_files_per_station(ifms_dir)

        # Process observations for each station
        output_per_station: dict[str, np.ndarray] = {}
        for station_name, station_files in ifms_per_station.items():

            if station_name == "DSS14":
                continue

            print(f"Processing observations for {station_name}")

            # Define observation collection
            contents, observations = define_observation_collection_for_station(
                bodies=bodies, station=station_name, station_ifms=station_files
            )

            # Define observation model for station
            link_definitions = observations.get_link_definitions_for_observables(
                tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
            )
            observation_model_settings = [
                tomss.model_settings.dsn_n_way_doppler_averaged(
                    link_ends=link_definitions[0],
                    subtract_doppler_signature=False,
                )
            ]

            # Create observation simulator
            simulator = tobss.observations_simulation_settings.create_observation_simulators(
                observation_settings=observation_model_settings, bodies=bodies
            )

            # Compute pre-fit residuals
            tobs.compute_residuals_and_dependent_variables(
                observation_collection=observations,
                observation_simulators=simulator,
                bodies=bodies,
            )

            # Save output
            obs_epochs = [
                eti.to_float() for eti in contents.observation_epochs_et
            ]
            results = np.array(
                [
                    observations.get_concatenated_observation_times(),
                    observations.get_concatenated_computed_observations(),
                    observations.get_concatenated_observations(),
                    obs_epochs,
                    contents.residuals,
                    observations.get_concatenated_residuals(),
                ]
            )
            np.save(output_dir / f"{station_name}", results)

            # observation_model_settings = [
            #     tomss.model_settings.dsn_n_way_doppler_averaged(
            #         link_ends=observations.link
            #     )
            # ]

            # # Extract raw data from IFMS files
            # observation_data = pio.load_doppler_observations_from_ifms_files(
            #     station_files
            # )

            # # Transform ramping reference epochs to TDB
            # __converter = ttime.default_time_scale_converter()
            # ramping_tdb0 = np.array(
            #     [
            #         __converter.convert_time_object(
            #             input_scale=ttime.TimeScales.utc_scale,
            #             output_scale=ttime.TimeScales.tdb_scale,
            #             input_value=ti,
            #             earth_fixed_position=bodies.get_body("Earth")
            #             .get_ground_station(station_name)
            #             .station_state.cartesian_positon_at_reference_epoch,
            #         )
            #         for ti in observation_data.ramping_utc0
            #     ]
            # )

            # # Define interpolator for transmitted frequency
            # interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
            #     start_times=ramping_tdb0[:-1].tolist(),
            #     end_times=ramping_tdb0[1:].tolist(),
            #     ramp_rates=observation_data.ramping_df[:-1].tolist(),
            #     start_frequency=observation_data.ramping_f0[:-1].tolist(),
            # )
            # station = bodies.get_body("Earth").get_ground_station(station_name)
            # station.set_transmitting_frequency_calculator(interpolator)

            # # Define ancilliary settings
            # ancilliary_settings = tobss.ancillary_settings.dsn_n_way_doppler_ancilliary_settings(
            #     frequency_bands=[
            #         tobss.ancillary_settings.FrequencyBands.x_band,
            #         tobss.ancillary_settings.FrequencyBands.x_band,
            #     ],
            #     reference_frequency_band=tobss.ancillary_settings.FrequencyBands.x_band,
            #     reference_frequency=0.0,
            #     integration_time=1.0,
            # )

            # # Define links ends
            # uplink_station = tomss.links.body_reference_point_link_end_id(
            #     body_name="Earth",
            #     reference_point_id=station_name,
            # )
            # transponder = tomss.links.body_origin_link_end_id(spacecraft_name)
            # # transponder = tomss.links.body_reference_point_link_end_id(
            # #     body_name=spacecraft_name,
            # #     reference_point_id="HGA",
            # # )
            # downlink_station = tomss.links.body_reference_point_link_end_id(
            #     body_name="Earth", reference_point_id=station_name
            # )
            # link_ends = {
            #     tomss.links.LinkEndType.transmitter: uplink_station,
            #     tomss.links.LinkEndType.retransmitter: transponder,
            #     tomss.links.LinkEndType.receiver: downlink_station,
            # }

            # # Define observation set
            # observable_type = (
            #     tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
            # )
            # station_observations = tobs.create_single_observation_set(
            #     observable_type=observable_type,
            #     link_ends=link_ends,
            #     observations=[
            #         np.array([x]) for x in observation_data.observation_values
            #     ],
            #     observation_times=observation_data.observation_epochs_et.tolist(),
            #     reference_link_end=tomss.links.LinkEndType.receiver,
            #     ancillary_settings=ancilliary_settings,
            # )
            # observations = tobs.ObservationCollection([station_observations])

            # # Define observation model settings
            # observation_model_settings: list[
            #     tomss.model_settings.ObservationModelSettings
            # ] = [
            #     tomss.model_settings.dsn_n_way_doppler_averaged(
            #         link_ends=observations.link_definitions_per_observable[
            #             observable_type
            #         ][0],
            #         subtract_doppler_signature=False,
            #         # light_time_correction_settings=light_time_settings,
            #     )
            # ]

            # # Create observation simulator
            # observation_simulator = tobss.observations_simulation_settings.create_observation_simulators(
            #     observation_settings=observation_model_settings, bodies=bodies
            # )

            # # Compute pre-fit residuals
            # tobs.compute_residuals_and_dependent_variables(
            #     observation_collection=observations,
            #     observation_simulators=observation_simulator,
            #     bodies=bodies,
            # )

            # # Compute pre-fit residuals and save output
            # output_per_station: dict[str, np.ndarray] = {}
            # for ifms_file in ifms_dir.glob("*.TAB"):

            #     print(f"Processing observations in {ifms_file.name}")

            #     # Load observations and compute pre-fit residuals
            #     _output = compute_prefit_residuals_for_ifms_file(
            #         bodies, ifms_file, spacecraft_name
            #     )
            #     if _output is None:
            #         print(f"Skipping {ifms_file.name}")
            #         continue
            #     station_name, observations = _output

            # Get residual distributed with observations
            # contents = tdata.read_ifms_file(str(ifms_file), False)
            # ifms_residuals = np.array(
            #     contents.raw_datamap["doppler_noise_hz"], dtype=float
            # )

            # # Get observation epochs and residuals as array
            # results = np.array(
            #     [
            #         observations.get_concatenated_observation_times(),
            #         observations.get_concatenated_computed_observations(),
            #         observations.get_concatenated_observations(),
            #         observation_data.residuals,
            #         observations.get_concatenated_residuals(),
            #     ]
            # )

            # np.save(output_dir / f"{station_name}.npy", results)

        #     # Add results to output container
        #     if station_name not in output_per_station:
        #         output_per_station[station_name] = results
        #     else:
        #         output_per_station[station_name] = np.append(
        #             output_per_station[station_name], results, axis=1
        #         )

        # for station, results in output_per_station.items():

        #     # Sort results by observation epoch
        #     sorted_results = results[:, np.argsort(results[0])]

        #     # Save output to file
        #     np.save(output_dir / f"{station}.npy", sorted_results)

    finally:
        spice.clear_kernels()
