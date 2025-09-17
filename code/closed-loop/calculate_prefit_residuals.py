from tudatpy.astro import (
    time_representation as ttime,
    element_conversion as telc,
)
from tudatpy.dynamics import environment_setup as tenvs, environment as tenv
from prefit import (
    utils as putils,
    io as pio,
    paths as ppaths,
    observations as pobs,
)
from tudatpy.interface import spice
from pathlib import Path
from tudatpy.estimation import (
    observations_setup as tobss,
    observable_models_setup as tomss,
    observations as tobs,
)
import numpy as np
from tudatpy.math import interpolators as tinterp
import shutil
import argparse
from matplotlib import pyplot as plt

Parser = argparse.ArgumentParser()
Parser.add_argument("config_file", help="Path to configuration file")


def define_system_of_bodies_from_raw_observations(
    configuration: pio.PrefitSettings,
    raw_observations_per_station: dict[str, pio.TwoWayDopplerObservations],
) -> tenv.SystemOfBodies:

    # Define initial and final epochs based on observation span
    initial_epoch_et = np.infty
    final_epoch_et = 0
    for raw_observations in raw_observations_per_station.values():
        if raw_observations.observation_epochs_et[0] < initial_epoch_et:
            initial_epoch_et = raw_observations.observation_epochs_et[0]
        if raw_observations.observation_epochs_et[-1] > final_epoch_et:
            final_epoch_et = raw_observations.observation_epochs_et[-1]
    assert isinstance(initial_epoch_et, ttime.Time)
    assert isinstance(final_epoch_et, ttime.Time)

    # General configuration: Time interval in UTC
    buffer_time = ttime.Time(configuration.time["buffer"])

    # General configuration: System of bodies
    global_frame_origin: str = "SSB"
    global_frame_orientation: str = "J2000"
    central_body: str = "Mars"
    default_step: ttime.Time = ttime.Time(3600.0)
    correction_bodies: list[str] = configuration.light_time["massive_bodies"]
    spacecraft_name: str = configuration.bodies["Spacecraft"]["name"]

    # Transform ET epochs to TDB
    initial_epoch_tdb = ttime.DateTime.from_julian_day(
        ttime.seconds_since_epoch_to_julian_day(initial_epoch_et.to_float())
    ).to_epoch_time_object()
    final_epoch_tdb = ttime.DateTime.from_julian_day(
        ttime.seconds_since_epoch_to_julian_day(final_epoch_et.to_float())
    ).to_epoch_time_object()

    # Apply buffer to initial and final epochs (UTC)
    initial_epoch_buffer = initial_epoch_tdb - buffer_time
    final_epoch_buffer = final_epoch_tdb + buffer_time

    # Initialize environment_settings
    environment_settings = tenvs.BodyListSettings(
        frame_orientation=global_frame_orientation,
        frame_origin=global_frame_origin,
    )

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

    # Position of HGA with respect to reference point of spacecraft
    mex_settings.rotation_model_settings = tenvs.rotation_model.spice(
        base_frame=global_frame_orientation,
        target_frame=configuration.bodies["Spacecraft"]["body_fixed_frame"],
    )

    # Define settings for Mars
    ############################################
    environment_settings.add_empty_settings("Mars")
    mars_settings = environment_settings.get("Mars")
    mars_config = configuration.bodies["Mars"]

    # Define settings for translational ephemerides
    match mars_config["translational_ephemerides"]:
        case "spice":
            mars_settings.ephemeris_settings = tenvs.ephemeris.direct_spice(
                frame_origin=global_frame_origin,
                frame_orientation=global_frame_orientation,
                body_name_to_use="Mars",
            )
        case "interpolated_spice":
            mars_settings.ephemeris_settings = (
                tenvs.ephemeris.interpolated_spice(
                    initial_time=initial_epoch_buffer,
                    final_time=final_epoch_buffer,
                    time_step=default_step,
                    frame_orientation=global_frame_orientation,
                    frame_origin=global_frame_origin,
                    body_name_to_use="Mars",
                )
            )
        case _:
            raise ValueError("Invalid ephemeris type for Mars")

    # Define basic mass settings
    mars_settings.gravity_field_settings = tenvs.gravity_field.central_spice(
        "Mars"
    )

    # Update settings for the Earth
    ############################################
    environment_settings.add_empty_settings("Earth")
    earth_settings = environment_settings.get("Earth")
    earth_config = configuration.bodies["Earth"]

    # Basic gravity field settings for light-time correction
    earth_settings.gravity_field_settings = tenvs.gravity_field.central_spice(
        "Earth"
    )

    # Define settings for rotation model
    match earth_config["rotation_model"]:
        case "spice":
            earth_settings.rotation_model_settings = tenvs.rotation_model.spice(
                base_frame=global_frame_orientation,
                target_frame="IAU_EARTH",
            )
        case "iers":
            cio_and_tdbtt_interp_settings = (
                tinterp.InterpolatorGenerationSettings(
                    interpolator_settings=tinterp.cubic_spline_interpolation(),
                    initial_time=initial_epoch_buffer,
                    final_time=final_epoch_buffer,
                    time_step=default_step,
                )
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
        case _:
            raise ValueError("Invalid rotation model for Earth")

    # Define settings for translational ephemerides
    match earth_config["translational_ephemerides"]:
        case "spice":
            earth_settings.ephemeris_settings = tenvs.ephemeris.direct_spice(
                frame_origin=global_frame_origin,
                frame_orientation=global_frame_orientation,
                body_name_to_use="Earth",
            )
        case "interpolated_spice":
            earth_settings.ephemeris_settings = (
                tenvs.ephemeris.interpolated_spice(
                    initial_time=initial_epoch_buffer,
                    final_time=final_epoch_buffer,
                    time_step=default_step,
                    frame_orientation=global_frame_orientation,
                    frame_origin=global_frame_origin,
                    body_name_to_use="Earth",
                )
            )
        case _:
            raise ValueError("Invalid ephemeris type for Earth")

    # Define shape settings (GRS80 - Recommended for VMF3)
    match earth_config["shape_model"]:
        case "grs80":
            earth_settings.shape_settings = tenvs.shape.oblate_spherical(
                equatorial_radius=6378136.6,
                flattening=(1.0 / 298.25642),
            )
        case "spherical_spice":
            earth_settings.shape_settings = tenvs.shape.spherical_spice()
        case _:
            raise ValueError("Invalid shape model for Earth")

    # # Define atmospheric settings
    # earth_settings.atmosphere_settings = tenvs.atmosphere.nrlmsise00()

    # # Define shape deformation settings (Displacements of ground stations)
    # earth_settings.shape_deformation_settings = [
    #     tenvs.shape_deformation.basic_solid_body_tidal(
    #         tide_raising_bodies=["Sun", "Moon"],
    #         displacement_love_numbers={2: (0.6, 0.08)},
    #     )
    # ]

    # Define settings for ground stations (TODO: Do this manually)
    stations_to_create = list(raw_observations_per_station.keys())
    ground_station_settings: list[
        tenvs.ground_station.GroundStationSettings
    ] = []
    for station_name in stations_to_create:

        # Get configuration for station
        if station_name not in configuration.stations:
            raise ValueError(f"Missing configuration for {station_name}")
        station_config = configuration.stations[station_name]

        # TODO: Control what to do based on station configuration
        # Reference epoch for ground station motion
        station_ref_epoch = ttime.DateTime.from_iso_string(
            "2000-01-01T00:00:00"
        )
        try:
            ref_state_itrf2000 = putils.get_station_reference_state_itrf2000(
                station_name
            )
        except KeyError:
            print(f"Using approximated position for station: {station_name}")
            _station_name = f"DSS-{station_name.split('DSS')[-1]}"
            ref_pos_itrf2000 = tenvs.ground_station.get_approximate_dsn_ground_station_positions()[
                _station_name
            ]
            ref_state_itrf2000 = np.array([*ref_pos_itrf2000, 0, 0, 0])

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

    # Add empty bodies for corrections
    existing_bodies = list(configuration.bodies.keys())
    empty_bodies = [
        body for body in correction_bodies if body not in existing_bodies
    ]
    for cbody in empty_bodies:

        # Add empty body
        environment_settings.add_empty_settings(cbody)
        cbody_settings = environment_settings.get(cbody)

        # Define interpolated ephemerides
        cbody_settings.ephemeris_settings = tenvs.ephemeris.interpolated_spice(
            initial_time=initial_epoch_buffer,
            final_time=final_epoch_buffer,
            time_step=default_step,
            frame_orientation=global_frame_orientation,
            frame_origin=global_frame_origin,
        )

        # Define point-mass gravity field settings from spice
        cbody_settings.gravity_field_settings = (
            tenvs.gravity_field.central_spice(cbody)
        )

    # Create system of bodies
    bodies = tenvs.create_system_of_bodies(environment_settings)

    # Set up spacecraft to use default turnaround ratios
    bodies.get_body(
        spacecraft_name
    ).system_models.set_default_transponder_turnaround_ratio_function()

    # Define positon of reference point within the spacecraft
    spacecraft_config = configuration.bodies["Spacecraft"]
    ref_point_name: str = spacecraft_config["reference_point"]
    match ref_point_name:
        case "MEX":
            pass
        case "HGA":

            # Get offset of COM wrt LVI if defined
            # Using this does not make sense, so I always set it to zero
            if "x_com_offset" not in spacecraft_config:
                print("Not using offsets")
                com_pos = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            else:
                com_pos = np.array(
                    [
                        float(spacecraft_config["x_com_offset"]),
                        float(spacecraft_config["y_com_offset"]),
                        float(spacecraft_config["z_com_offset"]),
                        0.0,
                        0.0,
                        0.0,
                    ]
                )

            cstate_hga_lvi_mexframe = (
                spice.get_body_cartesian_state_at_epoch(
                    target_body_name="MEX_HGA",
                    observer_body_name="MEX_SPACECRAFT",
                    reference_frame_name="MEX_SPACECRAFT",
                    aberration_corrections="NONE",
                    ephemeris_time=initial_epoch_et,
                )
                - com_pos
            )

            hga_ephemeris_settings = tenvs.ephemeris.constant(
                constant_state=cstate_hga_lvi_mexframe,
                frame_origin="MEX_SPACECRAFT",
                frame_orientation="MEX_SPACECRAFT",
            )
            hga_ephemeris = tenvs.ephemeris.create_ephemeris(
                ephemeris_settings=hga_ephemeris_settings, body_name="HGA"
            )
            bodies.get(spacecraft_name).system_models.set_reference_point(
                reference_point="HGA",
                ephemeris=hga_ephemeris,
            )
        case _:
            raise ValueError("Invalid reference point for spacecraft")

    # Add tropospheric data to the ground stations
    vmf3_path = Path().home() / ".pride/data/tropospheric"
    vmf3_files = [str(xi) for xi in vmf3_path.glob("*.v3gr_r")]
    tomss.light_time_corrections.set_vmf_troposphere_data(
        data_files=vmf3_files,
        file_has_meteo=True,
        file_has_gradient=True,
        bodies=bodies,
        set_troposphere_data=True,
        set_meteo_data=True,
    )

    # Add ionospheric data to the ground stations
    ionex_path = Path().home() / ".pride/data/ionospheric"
    ionex_files = [str(xi) for xi in ionex_path.glob("*.13i")]
    tomss.light_time_corrections.set_ionosphere_model_from_ionex(
        data_files=ionex_files,
        bodies=bodies,
    )

    # earth_body = bodies.get("Earth")
    # for ground_station in earth_body.ground_station_list:

    #     gs = earth_body.get_ground_station(ground_station)
    #     some = ttime.DateTime.from_modified_julian_day(56654.0)
    #     temp = gs.temperature_function(some.to_epoch()) - 273.15
    #     print(
    #         f"Station: {ground_station} - Epoch: {some.to_iso_string()} "
    #         f"- Temperature: {temp}"
    #     )

    return bodies


def group_ifms_data_per_station(
    source_dir: Path, ignore: list[str]
) -> dict[str, list[Path]]:

    # Get IFMS files in directory and group them per station
    ifms_per_station: dict[str, list[Path]] = {}
    for file in source_dir.glob("*.TAB"):

        ignore_file: bool = False
        for item in ignore:
            if item in file.name:
                ignore_file = True
                break
        if ignore_file:
            continue

        # Get metadata from file name
        metadata = pio.get_metadata_from_ifms_filename(file)

        # Get station name from metadata
        station_name = pobs.identify_station_from_id(
            int(metadata["station_id"])
        )

        # Update collection with file
        if station_name not in ifms_per_station:
            ifms_per_station[station_name] = [file]
        else:
            ifms_per_station[station_name].append(file)

    return ifms_per_station


def define_observation_collection_for_station(
    bodies: tenv.SystemOfBodies,
    station: str,
    station_raw_data: pio.TwoWayDopplerObservations,
    configuration: pio.PrefitSettings,
) -> tobs.ObservationCollection:

    # Load data from IFMS files into Python object
    content = station_raw_data

    # Define interpolator for uplink frequency
    interpolator = tenv.PiecewiseLinearFrequencyInterpolator(
        start_times=content.ramping_start_tdb.tolist(),
        end_times=content.ramping_stop_tdb.tolist(),
        ramp_rates=content.ramping_df.tolist(),
        start_frequency=content.ramping_f0.tolist(),
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
    match configuration.bodies["Spacecraft"]["reference_point"]:
        case "MEX":
            transponder = tomss.links.body_origin_link_end_id(spacecraft_name)
        case "HGA":
            transponder = tomss.links.body_reference_point_link_end_id(
                body_name=spacecraft_name,
                reference_point_id="HGA",
            )
        case _:
            raise ValueError("Invalid transponder setup from configuration")
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

    if configuration.observations["compress"]:
        raise NotImplementedError("Still incompatible with compression")

    # print(f"{station}", end="\n-----------------------------\n")
    # print(len(original_observations.get_concatenated_observations()))

    # # Compress observations if requested
    # if configuration.observations["compress"]:
    #     observations = (
    #         tobss.observations_wrapper.create_compressed_doppler_collection(
    #             original_observation_collection=original_observations,
    #             compression_ratio=int(
    #                 configuration.observations["integration_time"]
    #             ),
    #         )
    #     )
    # else:
    #     observations = original_observations

    # print(len(observations.get_concatenated_observations()))
    # print(f"{station}", end="\n-----------------------------\n")

    return observations


def load_spice_kernels(configuration: pio.PrefitSettings) -> None:

    # Define path to metakernel
    metakernel = ppaths.datadir / "metak.tm"

    # Get path to specific kernels
    match configuration.ephemerides["mex_translational"]:
        case "ORMM":
            extra_kernels = [
                "ORMM_T19_131201000000_01033.BSP",
                "ORMM_T19_140101000000_01041.BSP",
            ]
        case "ROB":
            extra_kernels = ["MEX_ROB_130101_131231_001.BSP"]
        case _:
            raise ValueError(f"Invalid option for mex ephemerides")

    # # Choose version of DE ephemerides
    # match configuration.ephemerides["de_version"]:
    #     case "DE405":
    #         print("Using DE405")
    #         extra_kernels.append("DE405.BSP")
    #         extra_kernels.append("../pck/DE403-MASSES.TPC")
    #     case "DE440":
    #         print("Using DE440")
    #         extra_kernels.append("DE440.BSP")
    #         extra_kernels.append("../pck/gm_de440.tpc")
    #     case _:
    #         raise ValueError("Invalid kernels")

    # Load standard kernels
    spice.load_kernel(str(metakernel))

    # Load extra kernels
    for extra_kernel in extra_kernels:
        spice.load_kernel(str(ppaths.kerneldir / extra_kernel))

    return None


if __name__ == "__main__":

    __config_path = str(Parser.parse_args().config_file)
    config_path = Path(__config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Not found: {config_path}")
    config = pio.load_configuration(config_path)

    # Paths
    metakernel = ppaths.datadir / "metak.tm"
    ifms_dir = ppaths.psadir
    output_dir = config_path.parent
    output_dir.mkdir(exist_ok=True, parents=True)

    # Copy configuration file to output directory
    dst_config_path = output_dir / config_path.name
    if config_path.resolve() != dst_config_path.resolve():
        shutil.copy(config_path, output_dir / config_path.name)

    try:

        # Load spice kernels based on configuration
        load_spice_kernels(config)

        # Load data from IFMS files and group per station
        ifms_paths_per_station = group_ifms_data_per_station(
            ifms_dir,
            ["133612305", "133642046", "M84", "M62", "M14", "M63"],
        )
        data_per_station: dict[str, pio.TwoWayDopplerObservations] = {}
        for _station, _station_ifms in ifms_paths_per_station.items():
            data_per_station[_station] = (
                pobs.load_doppler_observations_from_list_of_ifms_files(
                    _station_ifms, config.stations
                )
            )

        # Load data from ODF files and group per station
        odf_paths = [file for file in ppaths.psadir.glob("*.DAT")]
        data_per_station = pobs.load_odf_data_per_station(
            odf_paths, data_per_station
        )

        fig, ax = plt.subplots()
        nnorcia = data_per_station["NWNORCIA"]
        nnorcia_epochs = np.array(
            [xi.to_float() for xi in nnorcia.observation_epochs_et]
        )
        ax.plot(nnorcia_epochs, nnorcia.observation_values, ".")
        plt.show()
        exit(0)

        # Define system of bodies
        spacecraft_name: str = config.bodies["Spacecraft"]["name"]
        bodies = define_system_of_bodies_from_raw_observations(
            config,
            data_per_station,
        )

        # # Update system of bodies with data for light-time corrections
        # bodies, lt_correction_settings_dict = (
        #     update_system_of_bodies_for_light_time_corrections(
        #         system_of_bodies=bodies,
        #         troposphere=correct_troposphere,
        #         ionosphere=correct_ionosphere,
        #         relativistic=correct_relativistic,
        #     )
        # )

        # Process observations for each station
        output_per_station: dict[str, np.ndarray] = {}
        for station_name, contents in data_per_station.items():

            print(f"Processing observations for {station_name}")

            # Define an observation collection for the station
            observations = define_observation_collection_for_station(
                bodies=bodies,
                station=station_name,
                station_raw_data=contents,
                configuration=config,
            )

            # Define settings for light-time corrections
            lt_correction_settings: list[
                tomss.light_time_corrections.LightTimeCorrectionSettings
            ] = []
            station_config = config.stations[station_name]

            # Troposphere
            if (
                station_config["troposphere"]
                and not station_config["tropo_from_file"]
            ):

                lt_correction_settings.append(
                    tomss.light_time_corrections.vmf3_tropospheric_light_time_correction(
                        body_with_atmosphere_name="Earth",
                        use_gradient_correction=True,
                    )
                )

            # Ionosphere
            if station_config["ionosphere"]:

                lt_correction_settings.append(
                    tomss.light_time_corrections.ionex_ionospheric_light_time_correction(
                        body_with_ionosphere_name="Earth",
                        ionosphere_height=10.0,  # Irrelevant
                    )
                )

            # Relativistic
            if station_config["relativistic"]:

                lt_correction_settings.append(
                    tomss.light_time_corrections.first_order_relativistic_light_time_correction(
                        config.light_time["massive_bodies"]
                    )
                )

            # Define observation model for station
            link_definitions = observations.get_link_definitions_for_observables(
                tomss.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
            )
            observation_model_settings = [
                tomss.model_settings.dsn_n_way_doppler_averaged(
                    link_ends=link_definitions[0],
                    subtract_doppler_signature=False,
                    light_time_correction_settings=lt_correction_settings,
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
            _concatenated_observations = (
                observations.get_concatenated_observations()
            )
            results = np.array(
                [
                    observations.get_concatenated_observation_times(),
                    observations.get_concatenated_computed_observations(),
                    _concatenated_observations,
                    np.zeros_like(_concatenated_observations),
                    np.zeros_like(_concatenated_observations),
                    observations.get_concatenated_residuals(),
                ]
            )
            np.save(output_dir / f"{station_name}", results)

            print(f"Saved results for {station_name}")

    finally:
        spice.clear_kernels()
